import math
import numpy as np
from pathlib import Path
import numba
from numba import cuda as numba_cuda
from .clib import mds_clib
from typing import Optional

try:
    import cupy as cp
    _HAS_CUPY = True
except ImportError:
    _HAS_CUPY = False


nm = 1e-9
µm = 1e-6
scale = {}
scale['micrometer'] = µm
scale['nanometer'] = nm


# Thread block geometry for the numba.cuda kernel: 16x16 = 256 threads/block,
# one thread per sensor in the (Sy_range, Sx_range) grid.
NUMBA_CUDA_BLOCK = (16, 16)
# Dipoles are streamed through shared memory in tiles of this many. We size the
# tile to the block thread-count so the cooperative load is exactly one dipole
# per thread per tile (256 here). Shared use = 2 * TILE * 3 * 8 B = 12 KiB/block.
NUMBA_CUDA_TILE = NUMBA_CUDA_BLOCK[0] * NUMBA_CUDA_BLOCK[1]


# @numba.jit(nopython=True)
@numba.njit(parallel=True)
def dipole_B(dip_r, dip_m, Sx_range, Sy_range, Sheight, B_grid, dir_vector):
    """
    Compute the z-component of the dipole field at the node position(s) of a
    scan grid, from a group of particles located at the dip_r dipole_positions,
    and which have magnetic dipole moments given in the dip_m array. The scan
    grid is defined by rectangular sensors whose center are defined by the
    Sx_range and Sy_range points, so the scan surface is defined by a `Nx x Ny`
    grid of points.

    For these arrays, N > 1

    Parameters
    ----------
    dip_r
        N x 3 array with dipole positions (m)
    dip_m
        N x 3 array with dipole moments (Am^2)
    Sx_range
        M-array with coordinates of measurement point (m) in the x-direction
    Sy_range
        P-array with coordinates of measurement point (m) in the y-direction
    Sheight
        Height or z-position of the scan surface
    B_grid
        M x P array to be populated with the dipole field values
    dir_vector
        Normalized direction vector
    """
    if B_grid.shape[0] * B_grid.shape[1] != Sx_range.shape[0] * Sy_range.shape[0]:
        raise Exception("Bz grid array with wrong dimensions")

    pos_r = np.zeros(3)
    pos_r[2] = Sheight
    for j, sy in enumerate(Sy_range):
        for i, sx in enumerate(Sx_range):
            pos_r[0], pos_r[1] = sx, sy
            # Subtract every row of dip_r (Nx3 array) from pos_r (1x3 array):
            r = pos_r - dip_r

            x, y, z = r[:, 0], r[:, 1], r[:, 2]

            rho2 = np.sum(r**2, axis=1)
            rho = np.sqrt(rho2)
            rho5 = rho2 * rho2 * rho

            mx, my, mz = dip_m[:, 0], dip_m[:, 1], dip_m[:, 2]
            m_dot_r = mx * x + my * y + mz * z

            Bx = (3e-7 * x * m_dot_r / rho5) - (1e-7 * mx / (rho2 * rho))
            By = (3e-7 * y * m_dot_r / rho5) - (1e-7 * my / (rho2 * rho))
            Bz = (3e-7 * z * m_dot_r / rho5) - (1e-7 * mz / (rho2 * rho))

            B_grid[j, i] = np.sum(
                Bx * dir_vector[0] + By * dir_vector[1] + Bz * dir_vector[2]
            )

    return None


@numba_cuda.jit
def dipole_B_numba_cuda(dip_r, dip_m, Sx_range, Sy_range, Sheight,
                        B_grid, dir_vector):
    """
    GPU port of `dipole_B`: one thread per sensor in the (Sx, Sy) scan grid.
    The dipole list is streamed through shared memory in tiles of
    NUMBA_CUDA_TILE, so each dipole is read from global memory once per *block*
    instead of once per *thread* (the classic N-body staging). This is what
    lets the kernel keep scaling once the dipole arrays no longer fit in L2.

    All arithmetic is float64

    Optimizations:
      * Project onto `dir_vector` analytically instead of forming (Bx,By,Bz):
            B.d = 3e-7 (m.r)(r.d)/rho^5 - 1e-7 (m.d)/rho^3
      * Build all radial powers from a single reciprocal (1 sqrt + 1 divide per
        dipole) instead of six separate divisions
    """
    sh_r = numba_cuda.shared.array((NUMBA_CUDA_TILE, 3), numba.float64)
    sh_m = numba_cuda.shared.array((NUMBA_CUDA_TILE, 3), numba.float64)

    # get the thread position in the 2D gpu grid (which we set as sensor pos)
    i, j = numba_cuda.grid(2)
    # thread index
    tid = (numba_cuda.threadIdx.y * numba_cuda.blockDim.x + numba_cuda.threadIdx.x)

    # Out-of-range threads must still reach every syncthreads() below (an early
    # return here would deadlock the staged load), so we flag them and simply
    # skip the final store.
    inb = (i < Sx_range.shape[0]) and (j < Sy_range.shape[0])
    sx = Sx_range[i] if inb else 0.0
    sy = Sy_range[j] if inb else 0.0
    sz = Sheight

    dx = dir_vector[0]
    dy = dir_vector[1]
    dz = dir_vector[2]

    acc = 0.0
    Ndip = dip_r.shape[0]

    for base in range(0, Ndip, NUMBA_CUDA_TILE):
        # Here, each thread loads a dipole r and m, into the tile's shared mem
        # and wait for all threads until data is loaded
        k = base + tid
        if k < Ndip:  # skip non existent dipoles beyond Ndip in last block
            sh_r[tid, 0] = dip_r[k, 0]
            sh_r[tid, 1] = dip_r[k, 1]
            sh_r[tid, 2] = dip_r[k, 2]
            sh_m[tid, 0] = dip_m[k, 0]
            sh_m[tid, 1] = dip_m[k, 1]
            sh_m[tid, 2] = dip_m[k, 2]
        numba_cuda.syncthreads()

        upper = NUMBA_CUDA_TILE
        if Ndip - base < NUMBA_CUDA_TILE:
            upper = Ndip - base

        # Here, each thread computes the B field from the 256 dipole data
        # in shared memory
        for t in range(upper):
            rx = sx - sh_r[t, 0]
            ry = sy - sh_r[t, 1]
            rz = sz - sh_r[t, 2]
            mx = sh_m[t, 0]
            my = sh_m[t, 1]
            mz = sh_m[t, 2]

            rho2 = rx*rx + ry*ry + rz*rz
            inv_rho = 1.0 / math.sqrt(rho2)
            inv_rho2 = inv_rho * inv_rho
            inv_rho3 = inv_rho2 * inv_rho
            inv_rho5 = inv_rho3 * inv_rho2

            mr = mx*rx + my*ry + mz*rz
            r_d = rx*dx + ry*dy + rz*dz
            m_d = mx*dx + my*dy + mz*dz

            acc += 3e-7 * mr * r_d * inv_rho5 - 1e-7 * m_d * inv_rho3

        # Guard the shared buffers before the next tile overwrites them.
        numba_cuda.syncthreads()

    if inb:
        B_grid[j, i] = acc


if _HAS_CUPY:
    # cp.ElementwiseKernel: one thread per output sensor; the body is CUDA C
    # mirroring the CPU kernel one-to-one. `dip_r` / `dip_m` come in flat
    # (size 3*Ndipoles), so indices are 3*k + {0,1,2}. `sx`, `sy` are
    # per-element scalars; cupy handles the broadcast against Sx_range and
    # Sy_range at call time.
    dipole_B_cupy_kernel = cp.ElementwiseKernel(
        in_params=('raw float64 dip_r, raw float64 dip_m, '
                   'raw float64 dir_vec, '
                   'float64 sx, float64 sy, float64 sz, int32 Ndip'),
        out_params='float64 B',
        operation=r'''
            double dx = dir_vec[0];
            double dy = dir_vec[1];
            double dz = dir_vec[2];
            double acc = 0.0;
            for (int k = 0; k < Ndip; k++) {
                double rx = sx - dip_r[3*k + 0];
                double ry = sy - dip_r[3*k + 1];
                double rz = sz - dip_r[3*k + 2];

                double rho2 = rx*rx + ry*ry + rz*rz;
                double inv_rho  = 1.0 / sqrt(rho2);
                double inv_rho2 = inv_rho * inv_rho;
                double inv_rho3 = inv_rho2 * inv_rho;
                double inv_rho5 = inv_rho3 * inv_rho2;

                double mx = dip_m[3*k + 0];
                double my = dip_m[3*k + 1];
                double mz = dip_m[3*k + 2];

                double mr  = mx*rx + my*ry + mz*rz;
                double r_d = rx*dx + ry*dy + rz*dz;
                double m_d = mx*dx + my*dy + mz*dz;

                acc += 3e-7 * mr * r_d * inv_rho5 - 1e-7 * m_d * inv_rho3;
            }
            B = acc;
        ''',
        name='dipole_B_cupy',
    )

    # Pre-warm at import time: triggers NVRTC compilation up front so the
    # first compute_scan_signal(method="cupy") call doesn't pay the cost.
    # Subsequent Python processes on the same machine hit cupy's on-disk
    # kernel cache (~/.cupy/kernel_cache) and skip NVRTC entirely. If no GPU
    # is reachable at import time the warm-up silently skips; the cupy
    # method branch raises a clear error when actually invoked.
    try:
        _wu_dip = cp.zeros(3, dtype=np.float64)
        _wu_dir = cp.zeros(3, dtype=np.float64)
        print("Compiling CUDA kernel (cupy)...", flush=True)
        dipole_B_cupy_kernel(
            _wu_dip, _wu_dip, _wu_dir,
            np.float64(0.0), np.float64(0.0), np.float64(0.0),
            np.int32(0),  # Ndip=0 -> inner loop is a no-op, no math executed
        )
        cp.cuda.Stream.null.synchronize()
        del _wu_dip, _wu_dir
    except Exception:
        pass


class MicroDemagSignature(object):
    def __init__(
        self,
        scan_limits: list[list[2], list[2]],
        scan_spacing: list[2] | tuple[2],
        scan_height: float,
    ):
        """
        This class allows to calculate the dipolar signal from the
        magnetization vectors of a finite element micromagnetic system
        simulated with a micromagnetic code. The dipolar signal is recorded in a
        rectangular scan grid that is defined when instatiating this class.
        Information of the micromagnetic model is defined in a vbox file
        produced by MERRILL.
        The readers for different simulation code outputs, generate the common attributes:
            self.r
            self.dip_moments
        Other attributes are specific to the reader used.

        Parameters
        ----------
        scan_limits
            A sequence containing 2 pairs of floating values, defining the scan
            limits, e.g. [[0, 0], [1e-6, 1e-6]]. These points are the centers
            of the sensors in the west-south corner and the north-east corner
            respectively (see Notes). Specified in meters
        scan_spacing
            A sequence with 2 values for the scan grid spacing in the x and y
            directions (see Notes). Specified in meters
        scan_height
            Position of the scan grid in the z-direction. Specified in meters

        Notes
        -----
        For a scan grid the scan limit pairs define its corners as:

                -----------------
                |   |   |   | o | -- scan_limits[1]
                -----------------
                |   |   |   |   |
                -----------------
                | o |   |   |   |
                -----------------
                  |
              scan_limits[0]

        By specifying a scan_spacing these points are included. For instance,
        using `scan_spacing=[[0, 0], [1e-6, 1e-6]]` and
        `scan_spacing=[0.05e-6, 0.05e-6]` the grid has 21 points in each
        direction. If `scan_spacing` is not dividing the limits exactly by an
        integer number, then the scan grid limits are approximated. Check the
        `Sx_range` and `Sy_range` variables to check this.

        """

        self.scan_limits = scan_limits
        self.scan_spacing = scan_spacing
        self.scan_height = scan_height

        self.Nx = round((scan_limits[1][0] - scan_limits[0][0]) / scan_spacing[0]) + 1
        self.Ny = round((scan_limits[1][1] - scan_limits[0][1]) / scan_spacing[1]) + 1

        self.Sx = scan_limits[0][0] + np.arange(self.Nx) * scan_spacing[0]
        self.Sy = scan_limits[0][1] + np.arange(self.Ny) * scan_spacing[1]
        # self.Sx, self.Sy = np.meshgrid(self.Sx, self.Sy)
        # Sgrid = np.stack(
        #     (self.Sx, self.Sy, np.ones_like(self.Sx) * scan_height), axis=2)

        self.Ms = None
        self.B_grid = np.zeros((self.Ny, self.Nx), dtype=np.float64)

        # Lazy caches of the static device-side dipole arrays + output buffer
        # for the GPU backends. Populated on the first compute_scan_signal
        # call that uses the corresponding method. Call reset_gpu_cache() to
        # invalidate them after rerunning a reader (which replaces self.r
        # and self.dip_moments).
        self._numba_cuda_cache = None
        self._cupy_cache = None

    def reset_gpu_cache(self):
        """
        Drop any cached GPU-side copies of self.r / self.dip_moments / Sx /
        Sy / output buffer so the next compute_scan_signal call with a GPU
        backend re-uploads them. Call this after any operation that mutates
        the dipole arrays (e.g. re-running a reader on the same instance).
        """
        self._numba_cuda_cache = None
        self._cupy_cache = None

    def _read_magnetic_params(self, log_file):
        """
        Reads the saturation magnetization value from the log file of a MERRILL simulation.
        TODO: Might require to check MERRILL keeps the log file format intact
        """
        line = ""
        with open(log_file, "r") as f:
            while not line.startswith("Material"):
                line = f.readline()
            self.headers = f.readline().strip().split()
            mat_params = np.array(f.readline().strip().split(), dtype=np.float64)
            self.material_parameters = {}
        # for i, mp in enumerate(headers):
        #     self.material_parameters[mp] = mat_params[i]
        # print('Mat params:', self.material_parameters)
        self.Ms = mat_params[1]

    def reader_merrill_vbox(self,
                            mm_sim_file: str | Path,
                            energy_log_file: Optional[str] | Path = None,
                            Ms: Optional[float] = None,
                            origin_to_geom_center: bool = False,
                            vbox_file_delimiter=None
                            ) -> None:
        """
        Reads a MERRILL vbox file with 7 columns: x y z mx my mz volumes

        Parameters
        ----------
        mm_sim_file
            Path to a magnetization output file from a MERRILL vbox file:
                x y z mx my mz volume
            where the spatial data is scaled in µm
        merrill_energy_log_file (optional)
            Path to the MERRILL energy file from which magnetic paramters are
            read in the header
        Ms
            If None, the magnetization is computed using the energy log file
            specified in the main class. If Ms is passed, or self.Ms is
            specified, use the corresponding value
        origin_to_geom_center
            If True, all coordinates of the vbox file are shifted with respect
            to the geometric center of the system, which is computed using all
            coordinates and volumes from the file
        vbox_file_delimiter
            The delimiter for the vbox file from MERRILL passed to the
            `loadtxt` function. Prior to commit starting with hash 3aaf7f7
            (30/07/2021) the delimiter was whitespace, so the default option is
            `None`. Newer MERRILL versions use CSV format, so a comma `,` is
            required
        """
        mm_sim_file = Path(mm_sim_file)

        if energy_log_file:
            self._read_magnetic_params(energy_log_file)
        elif Ms:
            self.Ms = Ms

        if self.Ms is None:
            raise ValueError(
                "Specify a value for the saturation by either setting the Ms argument or via"
                " a log file in the constructor")

        # Read the vbox file: TODO: we can use a faster method for large files
        self.mag_data = np.loadtxt(
            mm_sim_file, skiprows=1, ndmin=2, delimiter=vbox_file_delimiter)
        # Scale spatial data:
        self.mag_data[:, :3] *= µm
        self.mag_data[:, 6] *= µm**3

        self.dip_volumes = self.mag_data[:, 6]
        # "Unpacking" occurs in the 1st dimension (row) so we transpose
        # The unpacking should generate mem views of the arrays
        self.x, self.y, self.z = self.mag_data[:, :3].T
        self.r = self.mag_data[:, :3]
        self.mx, self.my, self.mz = self.mag_data[:, 3:6].T

        # Shift positions wrt to the geometric centre if True
        if origin_to_geom_center:
            geom_center = self.r * self.dip_volumes[:, np.newaxis]
            geom_center = geom_center.sum(axis=0)
            geom_center = geom_center / self.dip_volumes.sum()

            np.subtract(self.r, geom_center, out=self.r)

        # self.dip_moments = self.Ms * self.dip_volumes
        # # This requires self.dip_moments to have 3 columns :/ so it fails:
        # np.multiply(self.dip_moments[:, np.newaxis], self.mag_data[:, 3:6],
        #             out=self.dip_moments)
        self.dip_moments = (
            self.Ms * self.dip_volumes[:, np.newaxis] * self.mag_data[:, 3:6]
        )
        # Dipole arrays were just replaced; drop any stale GPU buffers.
        self.reset_gpu_cache()

    def reader_finmag(self,
                      mm_sim_file: str | Path,
                      Ms: float,
                      origin_to_geom_center: bool,
                      delimiter: Optional[str] = None,
                      units: str = 'nanometer',
                      traslation_vector: Optional[list[3]] = None
                      ) -> None:
        """
        Reads a finite element micromagnetic file produced by the code FinMag
        The file must contain 6 columns with spin data per mesh node:

            x y z mx my mz cell_volume

        where the cell volume is the effective volume per node averaged using
        neighboring tetrahedra
        This reader assumes the magnetization is non-zero in all sites, and that
        there is a single material

        Parameters
        ----------
        Ms
            From the micromagnetic simulation (assuming single material)
        origin_to_geom_center
            If True, all coordinates of the vbox file are shifted with respect
            to the geometric center of the system, which is computed using all
            coordinates from the file, and the cell volumes based on dx, dy, dz
        delimiter
            For the text file containing the cell positions and magnetizations
        units
            Scale units of positions. Usually in nm in micromagnetic codes
        traslation_vector
            3-element list with coordinates (same scale units than in the file)
            to translate center after (if specified) shifting the origin to
            the geometric center
        """
        # TODO: implement for different Ms in mesh, maybe accepting an additional column

        mm_sim_file = Path(mm_sim_file)
        if mm_sim_file.suffix in ['.txt', '.dat']:
            self.mag_data = np.loadtxt(mm_sim_file, ndmin=2, delimiter=delimiter)
        else:
            self.mag_data = np.load(mm_sim_file)  # numpy handles io errors

        # Check and discard points with zero magnetization
        # eps = 1e-4
        # materialCells = np.linalg.norm(self.mag_data[:, 3:6], axis=1) > eps
        # self.mag_data = self.mag_data[materialCells]

        # "Unpacking" occurs in the 1st dimension (row) so we transpose
        # The unpacking should generate mem views of the arrays
        self.x, self.y, self.z = self.mag_data[:, :3].T
        self.r = self.mag_data[:, :3]
        self.mx, self.my, self.mz = self.mag_data[:, 3:6].T

        self.fe_ntets = self.r.shape[0]
        # Compute cell and particle volume
        self.fe_tet_volumes = self.mag_data[:, 6]
        self.mesh_volume = np.sum(self.fe_tet_volumes)

        # Current geometric center
        self.geom_center = self.r * self.fe_tet_volumes[:, np.newaxis]
        self.geom_center = self.geom_center.sum(axis=0) / self.mesh_volume

        # Shift positions wrt to the geometric centre if True
        if origin_to_geom_center:
            np.subtract(self.r, self.geom_center, out=self.r)

            # Recompute gc: (should be zero)
            self.geom_center = self.r * self.fe_tet_volumes[:, np.newaxis]
            self.geom_center = self.geom_center.sum(axis=0) / self.mesh_volume

        # Translate positions if specified
        if traslation_vector:
            traslation_vector = np.array(traslation_vector)
            self.mag_data[:, :3] += traslation_vector
            # Shift geom center
            np.add(self.geom_center, traslation_vector, out=self.geom_center)

        # Scale spatial data:
        self.mag_data[:, :3] *= scale[units]
        self.fe_tet_volumes *= scale[units]**3
        self.mesh_volume *= scale[units]**3
        self.geom_center *= scale[units]

        self.dip_moments = Ms * self.fe_tet_volumes[:, np.newaxis] * self.mag_data[:, 3:6]
        # Dipole arrays were just replaced; drop any stale GPU buffers.
        self.reset_gpu_cache()

    def reader_fd_micromagnetic(self,
                                mm_sim_file: str | Path,
                                Ms: float,
                                origin_to_geom_center: bool,
                                dV: list[3],
                                delimiter: Optional[str] = None,
                                units: str = 'nanometer',
                                traslation_vector: Optional[list[3]] = None
                                ) -> None:
        """
        Reads a finite difference micromagnetic file with 6 columns: x y z mx my mz
        Sites with zero magnetization are discarded
        The number of cells is computed from sites with nonvanishing magnetization

        Parameters
        ----------
        Ms
            From the micromagnetic simulation (assuming single material)
        origin_to_geom_center
            If True, all coordinates of the vbox file are shifted with respect
            to the geometric center of the system, which is computed using all
            coordinates from the file, and the cell volumes based on dx, dy, dz
        dV
            3-element list with the cell dimensions: dx, dy, dz
            Used to compute volume
        delimiter
            For the text file containing the cell positions and magnetizations
        units
            Scale units of positions. Usually in nm in micromagnetic codes
        traslation_vector
            3-element list with coordinates (same scale units than in the file)
            to translate center after (if specified) shifting the origin to
            the geometric center
        """
        # if str(self.mm_sim_file).endswith('npy'):  # or use Path's suffix

        mm_sim_file = Path(mm_sim_file)
        if mm_sim_file.suffix in ['.txt', '.dat']:
            self.mag_data = np.loadtxt(mm_sim_file, ndmin=2, delimiter=delimiter)
        else:
            self.mag_data = np.load(mm_sim_file)  # numpy handles io errors

        # Check and discard points with zero magnetization
        eps = 1e-4
        materialCells = np.linalg.norm(self.mag_data[:, 3:6], axis=1) > eps
        self.mag_data = self.mag_data[materialCells]

        # "Unpacking" occurs in the 1st dimension (row) so we transpose
        # The unpacking should generate mem views of the arrays
        self.x, self.y, self.z = self.mag_data[:, :3].T
        self.r = self.mag_data[:, :3]
        self.mx, self.my, self.mz = self.mag_data[:, 3:6].T

        self.fd_ncells = self.r.shape[0]
        # Compute cell and particle volume (using sites with |m| > 0)
        self.fd_cell_volume = dV[0] * dV[1] * dV[2]
        self.mesh_volume = self.fd_cell_volume * self.fd_ncells

        # Current geometric center
        self.geom_center = self.r.sum(axis=0)
        self.geom_center = self.geom_center * self.fd_cell_volume / self.mesh_volume

        # Shift positions wrt to the geometric centre if True
        if origin_to_geom_center:
            np.subtract(self.r, self.geom_center, out=self.r)

            # Recompute gc: (should be zero)
            self.geom_center = self.r.sum(axis=0)
            self.geom_center = self.geom_center * self.fd_cell_volume / self.mesh_volume

        # Translate positions if specified
        if traslation_vector:
            traslation_vector = np.array(traslation_vector)
            self.mag_data[:, :3] += traslation_vector
            # Shift geom center
            np.add(self.geom_center, traslation_vector, out=self.geom_center)

        # Scale spatial data:
        self.mag_data[:, :3] *= scale[units]
        self.fd_cell_volume *= scale[units]**3
        self.mesh_volume *= scale[units]**3
        self.geom_center *= scale[units]

        self.dip_moments = Ms * self.fd_cell_volume * self.mag_data[:, 3:6]
        # Dipole arrays were just replaced; drop any stale GPU buffers.
        self.reset_gpu_cache()

    def compute_scan_signal(self, method="numba", direction_vector=(0.0, 0.0, 1.0)):
        """
        Computes the dipolar signal at the scan surface

        Parameters
        ----------
        method
            Backend used to evaluate the dipole-field sum. One of:
              * `numba`       — CPU, numba @njit(parallel=True)
              * `cython`      — CPU, C/OpenMP via the bundled C library
              * `numba_cuda`  — GPU, numba.cuda kernel (one thread per sensor)
              * `cupy`        — GPU, cupy.ElementwiseKernel (one thread per sensor)
            `numba_cuda` requires a CUDA-capable GPU and numba's CUDA runtime;
            `cupy` additionally requires the `cupy` package. Results are saved
            in `self.B_grid`.
        direction_vector
            Vector with the direction of the field to be computed. By default, the
            z-component of the field is used.
        """
        # Ensure normalized dir vector
        dir_vector = np.array(direction_vector)
        dv_norm = np.linalg.norm(dir_vector)
        dir_vector = dir_vector / dv_norm

        # print(self.dip_moments.shape)
        if method == "numba":
            dipole_B(
                self.r,
                self.dip_moments,
                self.Sx,
                self.Sy,
                self.scan_height,
                self.B_grid,
                dir_vector,
            )

        elif method == "cython":
            r = self.r.ravel()
            m = self.dip_moments.ravel()
            mds_clib.dipole_bz_field_C(
                r,
                m,
                self.dip_moments.shape[0],
                self.Sx,
                self.Sy,
                self.Sx.shape[0],
                self.Sy.shape[0],
                self.scan_height,
                self.B_grid,
                dir_vector
            )

        elif method == "numba_cuda":
            if not numba_cuda.is_available():
                raise RuntimeError(
                    "method='numba_cuda' requires a CUDA-capable GPU and a "
                    "working numba.cuda runtime.")

            # The dipole arrays, scan axes and output buffer are static across
            # calls on this instance, so upload them once and cache on
            # self._numba_cuda_cache. Call reset_gpu_cache() if you re-run a
            # reader (which replaces self.r / self.dip_moments) on the same
            # instance. Only dir_vector (3 floats) is re-uploaded per call,
            # since direction_vector may change between invocations.
            cache = self._numba_cuda_cache
            if cache is None:
                # ascontiguousarray(dtype=float64) is a no-op when the array is
                # already float64+contiguous, and defends against float32
                # readers otherwise.
                r_f64 = np.ascontiguousarray(self.r, dtype=np.float64)
                m_f64 = np.ascontiguousarray(self.dip_moments, dtype=np.float64)
                cache = {
                    "d_r": numba_cuda.to_device(r_f64),
                    "d_m": numba_cuda.to_device(m_f64),
                    "d_Sx": numba_cuda.to_device(
                        np.ascontiguousarray(self.Sx, dtype=np.float64)),
                    "d_Sy": numba_cuda.to_device(
                        np.ascontiguousarray(self.Sy, dtype=np.float64)),
                    "d_B": numba_cuda.device_array_like(self.B_grid),
                }
                self._numba_cuda_cache = cache

            d_dir = numba_cuda.to_device(
                np.ascontiguousarray(dir_vector, dtype=np.float64))

            tpb = NUMBA_CUDA_BLOCK
            bpg = ((self.Sx.shape[0] + tpb[0] - 1) // tpb[0],
                   (self.Sy.shape[0] + tpb[1] - 1) // tpb[1])

            dipole_B_numba_cuda[bpg, tpb](
                cache["d_r"], cache["d_m"], cache["d_Sx"], cache["d_Sy"],
                float(self.scan_height), cache["d_B"], d_dir,
            )
            cache["d_B"].copy_to_host(self.B_grid)

        elif method == "cupy":
            if not _HAS_CUPY:
                raise ImportError(
                    "method='cupy' requires the `cupy` package to be installed."
                )
            if cp.cuda.runtime.getDeviceCount() == 0:
                raise RuntimeError(
                    "method='cupy' requires a CUDA-capable GPU; cupy is "
                    "installed but no CUDA device is visible."
                )

            # Static device arrays cached on self._cupy_cache (see the
            # numba_cuda branch above and reset_gpu_cache()). Sx/Sy are stored
            # pre-broadcast to (1, Nx) / (Ny, 1) so cupy expands them to
            # (Ny, Nx) against the output, giving each element its own
            # (sx, sy). dtype=float64 also promotes any float32 reader input.
            cache = self._cupy_cache
            if cache is None:
                cache = {
                    "dip_r": cp.asarray(self.r.ravel(), dtype=np.float64),
                    "dip_m": cp.asarray(self.dip_moments.ravel(),
                                        dtype=np.float64),
                    "Sx": cp.asarray(self.Sx, dtype=np.float64)[None, :],
                    "Sy": cp.asarray(self.Sy, dtype=np.float64)[:, None],
                    "B": cp.empty((self.Sy.shape[0], self.Sx.shape[0]),
                                  dtype=np.float64),
                }
                self._cupy_cache = cache

            dir_d = cp.asarray(dir_vector, dtype=np.float64)

            dipole_B_cupy_kernel(
                cache["dip_r"], cache["dip_m"], dir_d,
                cache["Sx"], cache["Sy"],
                np.float64(self.scan_height),
                np.int32(self.r.shape[0]),
                cache["B"],
            )

            self.B_grid[:] = cp.asnumpy(cache["B"])
