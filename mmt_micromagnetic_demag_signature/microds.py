import numpy as np
from pathlib import Path
import numba
from .clib import mds_clib


nm = 1e-9
µm = 1e-6


# @numba.jit(nopython=True)
@numba.njit(parallel=True)
def dipole_Bz(dip_r, dip_m, Sx_range, Sy_range, Sheight, Bz_grid):
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
    Bz_grid
        M x P array to be populated with the dipole field values
    """
    if Bz_grid.shape[0] * Bz_grid.shape[1] != Sx_range.shape[0] * Sy_range.shape[0]:
        raise Exception('Bz grid array with wrong dimensions')

    pos_r = np.zeros(3)
    pos_r[2] = Sheight
    for j, sy in enumerate(Sy_range):
        for i, sx in enumerate(Sx_range):
            pos_r[0], pos_r[1] = sx, sy
            # Subtract every row of dip_r (Nx3 array) from pos_r (1x3 array):
            r = pos_r - dip_r

            x, y, z = r[:, 0], r[:, 1], r[:, 2]

            rho2 = np.sum(r ** 2, axis=1)
            rho = np.sqrt(rho2)

            mx, my, mz = dip_m[:, 0], dip_m[:, 1], dip_m[:, 2]
            m_dot_r = mx * x + my * y + mz * z
            f = 3e-7 * z * m_dot_r / (rho2 * rho2 * rho)
            g = -1e-7 * mz / (rho2 * rho)

            # Only return Bz
            res = f + g

            Bz_grid[j, i] = np.sum(res)

    return None


class MicroDemagSignature(object):

    def __init__(self,
                 scan_limits,
                 scan_spacing,
                 scan_height,
                 merrill_mag_vbox_file,
                 merrill_energy_log_file=None
                 ):
        """
        This class allows to calculate the dipolar signal from the
        magnetization vectors of a finite element micromagnetic system
        simulated with the MERRILL code. The dipolar signal is recorded in a
        rectangular scan grid that is defined when instatiating this class.
        Information of the micromagnetic model is defined in a vbox file
        produced by MERRILL.

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
        merrill_mag_vbox_file
            Path to a MERRILL's vbox output file with 7 columns:
                x y z mx my mz volume
            where the spatial data is scaled in µm
        merrill_energy_log_file (optional)
            Path to the MERRILL energy file from which magnetic paramters are
            read in the header

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
        self.mag_vbox_file = merrill_mag_vbox_file
        self.energy_log_file = merrill_energy_log_file

        self.Nx = round((scan_limits[1][0] - scan_limits[0][0]) / scan_spacing[0]) + 1
        self.Ny = round((scan_limits[1][1] - scan_limits[0][1]) / scan_spacing[1]) + 1

        self.Sx = scan_limits[0][0] + np.arange(self.Nx) * scan_spacing[0]
        self.Sy = scan_limits[0][1] + np.arange(self.Ny) * scan_spacing[1]
        # self.Sx, self.Sy = np.meshgrid(self.Sx, self.Sy)
        # Sgrid = np.stack(
        #     (self.Sx, self.Sy, np.ones_like(self.Sx) * scan_height), axis=2)

        self.Ms = None

    def _read_magnetic_params(self, log_file):
        """
        Reads the saturation magnetization value from the log file of a MERRILL
        simulation.
        TODO: Might require to check MERRILL keeps the log file format intact

        """
        line = ''
        with open(log_file, 'r') as f:
            while not line.startswith('Material'):
                line = f.readline()
            self.headers = f.readline().strip().split()
            mat_params = np.array(f.readline().strip().split(),
                                  dtype=np.float64)
            self.material_parameters = {}
        # for i, mp in enumerate(headers):
        #     self.material_parameters[mp] = mat_params[i]
        # print('Mat params:', self.material_parameters)
        self.Ms = mat_params[1]

    def read_input_files(self,
                         Ms=None,
                         origin_to_geom_center=False,
                         vbox_file_delimiter=None):
        """
        Parameters
        ----------
        Ms
            If None, the magnetization is computed using the energy log file
            specified in the main class. If Ms is passed or self.Ms is
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

        if self.energy_log_file:
            self._read_magnetic_params(self.energy_log_file)
        elif Ms or self.Ms:
            self.Ms = Ms

        if self.Ms is None:
            raise ValueError('Specify a value for the saturation by either '
                             'setting the Ms argument or via a log file in '
                             'the constructor')

        # Read the vbox file: TODO: we can use a faster method for large files
        self.mag_data = np.loadtxt(self.mag_vbox_file, skiprows=1, ndmin=2,
                                   delimiter=vbox_file_delimiter)
        # Scale spatial data:
        self.mag_data[:, :3] *= µm
        self.mag_data[:, 6] *= (µm ** 3)

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
        self.dip_moments = self.Ms * self.dip_volumes[:, np.newaxis] * self.mag_data[:, 3:6]

        self.Bz_grid = np.zeros((self.Ny, self.Nx), dtype=np.float64)

    def compute_scan_signal(self, method='numba'):
        """
        Computes the dipolar signal at the scan surface

        Parameters
        ----------
        method
            Specify `numba` or `cython` for the calculations. The C method is
            parallelized with OpenMP. Results are saved in the `self.Bz_grid`
            array.

        """
        # print(self.dip_moments.shape)
        if method == 'numba':
            dipole_Bz(self.r, self.dip_moments, self.Sx, self.Sy,
                      self.scan_height, self.Bz_grid)
        elif method == 'cython':
            r = self.r.ravel()
            m = self.dip_moments.ravel()
            mds_clib.dipole_bz_field_C(r, m, self.dip_moments.shape[0],
                                       self.Sx, self.Sy,
                                       self.Sx.shape[0], self.Sy.shape[0],
                                       self.scan_height, self.Bz_grid)
