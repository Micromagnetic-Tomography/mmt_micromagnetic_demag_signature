import numpy as np
from pathlib import Path
import numba


nm = 1e-9
µm = 1e-6


@numba.jit(nopython=True)
def dipole_Bz(dip_r, dip_m, Sx_range, Sy_range, Sheight, Bz_grid):
    """
    Compute the z-component of the dipole field at the pos_r position(s), from
    a group of particles located at the dip_r dipole_positions, and which have
    magnetic dipole moments given in the dip_m array.

    For these arrays, N > 1

    Parameters
    ----------

    dip_r
        N x 3 array with dipole positions (m)
    dip_m
        N x 3 array with dipole moments (Am^2)
    pos_r
        M X P x 3 array (grid) with coordinates of measurement point (m)

    """
    if Bz_grid.shape[0] * Bz_grid.shape[1] != Sx_range.shape[0] * Sy_range.shape[0]:
        raise Exception('Bz grid array with wrong dimensions')

    # For every row of dip_r (Nx3 array), subtract pos_r (1x3 array)
    pos_r = np.zeros(3)
    pos_r[2] = Sheight
    for j, sy in enumerate(Sy_range):
        for i, sx in enumerate(Sx_range):
            pos_r[0], pos_r[1] = sx, sy
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
        """TODO: to be defined.

        :Docstring for MMTDemagsignal.: TODO

        Parameters
        ----------

        merrill_mag_vbox_file
            Path to the MERRILL output file with 7 columns:
                x y z mx my mz volume
            where the spatial data is scaled in µm
        merrill_energy_log_file (optional)
            Path to the MERRILL energy file from which magnetic paramters are
            read from the header

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
        Just assigning Ms for now
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

    def read_input_files(self, Ms=None):
        """
        """

        if self.energy_log_file:
            self._read_magnetic_params(self.energy_log_file)
        elif Ms:
            self.Ms = Ms

        if self.Ms is None:
            raise ValueError('Specify a value for the saturation by either '
                             'setting the Ms argument or via a log file in '
                             'the constructor')

        self.mag_data = np.loadtxt(self.mag_vbox_file, skiprows=1, ndmin=2)
        # Scale spatial data:
        self.mag_data[:, :3] *= nm
        self.mag_data[:, 6] *= (µm ** 3)

        self.dip_volumes = self.mag_data[:, 6]
        # "Unpacking" occurs in the 1st dimension (row) so we transpose
        # The unpacking should generate mem views of the arrays
        self.x, self.y, self.z = self.mag_data[:, :3].T
        self.r = self.mag_data[:, :3]
        self.mx, self.my, self.mz = self.mag_data[:, 3:6].T

        # self.dip_moments = self.Ms * self.dip_volumes
        # # This requires self.dip_moments to have 3 columns :/ so it fails:
        # np.multiply(self.dip_moments[:, np.newaxis], self.mag_data[:, 3:6],
        #             out=self.dip_moments)
        self.dip_moments = self.Ms * self.dip_volumes[:, np.newaxis] * self.mag_data[:, 3:6]

        self.Bz_grid = np.zeros((self.Ny, self.Nx), dtype=np.float64)

    def compute_scan_signal(self):
        """
        """
        # print(self.dip_moments.shape)
        dipole_Bz(self.r, self.dip_moments, self.Sx, self.Sy, self.scan_height,
                  self.Bz_grid)

