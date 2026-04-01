import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# import json
# import collections
import mmt_micromagnetic_demag_signature as mds

µm = 1e-6
nm = 1e-9


def test_cube_50nm_merrill_vs_analytic():
    # M = np.array([1, 1, -1.])
    # print(np.linalg.norm(M))
    # print(M / np.linalg.norm(M))
    print('** Testing MERRILL finite element simulation file **')

    # MERRILL Sim -----------------------------------------------------------------

    # Creating scan grid
    scan_spacing = (10 * nm, 10 * nm)
    scan_limits = np.array([[-1.5, -1.5], [1.5, 1.5]]) * µm
    scan_height = 500 * nm

    # Get the absolute path of the script's directory
    script_dir = Path(__file__).resolve().parent
    FILE_energy = script_dir / "cuboid_merrill_test/cube_sim_energy.log"
    FILE_sim = script_dir / "cuboid_merrill_test/cube_sim_magnetisation_volume.vbox"

    # Notice that the cube in the MERRILL simulation is already centered
    # at z = -35 nm
    demag_signal = mds.MicroDemagSignature(scan_limits, scan_spacing, scan_height)
    demag_signal.reader_merrill_vbox(FILE_sim,
                                     FILE_energy,
                                     origin_to_geom_center=False
                                     )
    demag_signal.compute_scan_signal(method="cython")

    # bz_grid = np.copy(demag_signal.Bz_grid)
    # np.save('./cube_50x50x50_test_demag_signal.npy', bz_grid)

    # -----------------------------------------------------------------------------

    # Analytic solution in µT units
    analytic_sol = np.load(
        script_dir / "analytical_cuboid_code/cuboid_50nm_centre-at_-35nm_scan-grid_3microm_m_11-1_Bz.npy"
    )

    # -----------------------------------------------------------------------------
    # Relative error
    err = np.linalg.norm(analytic_sol - demag_signal.B_grid * 1e6, ord="fro")
    relerr = err / np.linalg.norm(analytic_sol, ord="fro")

    # Norm of the simulation with the analytical Bz matrix should be less than:
    assert relerr < 1e-4

    print("Max MERRILL (µT)   :", np.abs(demag_signal.B_grid).max() * 1e6)
    print("Max Analytic (µT)  :", np.abs(analytic_sol).max())
    print(
        "Max difference (µT):", np.abs(analytic_sol - demag_signal.B_grid * 1e6).max()
    )
    print("Error Diff (µT)    :", err)
    print("Rel Error Diff %   :", 100 * relerr)


def test_cube_50nm_fdMicromagnetics_vs_analytic():

    print('** Testing mumax+ finite difference simulation file **')

    # Creating scan grid
    scan_spacing = (10 * nm, 10 * nm)
    scan_limits = np.array([[-1.5, -1.5], [1.5, 1.5]]) * µm
    scan_height = 500 * nm

    script_dir = Path(__file__).resolve().parent
    # FILE_energy = script_dir / "cuboid_merrill_test/cube_sim_energy.log"

    # Notice that the cube in the MERRILL simulation is already centered
    # at z = -35 nm
    demag_signal = mds.MicroDemagSignature(scan_limits, scan_spacing, scan_height)

    sim_file = script_dir / "cuboid_fd_mm_test/cube_mumaxPlus_L_50nm_centerAt_-35nm_dxyz_2nm.npy"
    demag_signal.reader_fd_micromagnetic(sim_file,
                                         Ms=4.8e5,
                                         origin_to_geom_center=False,
                                         dV=[2, 2, 2], n=[25, 25, 25],
                                         units='nanometer',
                                         traslation_vector=[0, 0, 0]
                                         )
    demag_signal.compute_scan_signal(method="cython")

    # -----------------------------------------------------------------------------

    # Analytic solution in µT units
    analytic_sol = np.load(
        script_dir / "analytical_cuboid_code/cuboid_50nm_centre-at_-35nm_scan-grid_3microm_m_11-1_Bz.npy"
    )

    # -----------------------------------------------------------------------------
    # Relative error
    err = np.linalg.norm(analytic_sol - demag_signal.B_grid * 1e6, ord="fro")
    relerr = err / np.linalg.norm(analytic_sol, ord="fro")

    # Norm of the simulation with the analytical Bz matrix should be less than:
    # assert relerr < 1e-4

    print("Max Simulation (µT)  :", np.abs(demag_signal.B_grid).max() * 1e6)
    print("Max Analytic   (µT)  :", np.abs(analytic_sol).max())
    print(
        "Max difference (µT):", np.abs(analytic_sol - demag_signal.B_grid * 1e6).max()
    )
    print("Error Diff (µT)    :", err)
    print("Rel Error Diff %   :", 100 * relerr)


if __name__ == "__main__":
    test_cube_50nm_merrill_vs_analytic()
    print('-' * 80 + '\n')
    test_cube_50nm_fdMicromagnetics_vs_analytic()
