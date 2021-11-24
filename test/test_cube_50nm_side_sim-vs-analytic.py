import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
# import json
# import collections
import micromagnetic_demag_signature as mds

µm = 1e-6
nm = 1e-9


def test_cube_50nm_merrill_vs_analytic():

    # M = np.array([1, 1, -1.])
    # print(np.linalg.norm(M))
    # print(M / np.linalg.norm(M))

    # MERRILL Sim -----------------------------------------------------------------

    # Creating scan grid
    scan_spacing = (10 * nm, 10 * nm)
    scan_limits = np.array([[-1.5, -1.5], [1.5, 1.5]]) * µm
    scan_height = 500 * nm

    FILE_energy = Path('cuboid_merrill_test/cube_sim_energy.log')

    # Notice that the cube in the MERRILL simulation is already centered
    # at z = -35 nm
    demag_signal = mds.MicroDemagSignature(
        scan_limits, scan_spacing, scan_height,
        'cuboid_merrill_test/cube_sim_magnetisation_volume.vbox',
        FILE_energy)

    demag_signal.read_input_files()
    demag_signal.compute_scan_signal(method='cython')

    # bz_grid = np.copy(demag_signal.Bz_grid)
    # np.save('./cube_50x50x50_test_demag_signal.npy', bz_grid)

    # -----------------------------------------------------------------------------

    # Analytic solution
    analytic_sol = np.load('analytical_cuboid_code/' +
        'cuboid_50nm_centre-at_-35nm_scan-grid_3microm_m_11-1_Bz.npy')

    # -----------------------------------------------------------------------------
    # Relative error
    err = np.linalg.norm(analytic_sol - demag_signal.Bz_grid * 1e6, ord='fro')
    err = err / np.linalg.norm(analytic_sol, ord='fro')

    # Norm of the simulation with the analytical Bz matrix should be less than:
    assert err < 1e-4

    print('Max MERRILL   :', np.abs(demag_signal.Bz_grid).max() * 1e6)
    print('Max Analytic  :', np.abs(analytic_sol).max())
    print('Max difference:', np.abs(analytic_sol -
                                    demag_signal.Bz_grid * 1e6).max())


if __name__ == "__main__":
    test_cube_50nm_merrill_vs_analytic()
