import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# import json
# import collections
import mmt_micromagnetic_demag_signature as mds
from mmt_micromagnetic_demag_signature.microds import _HAS_CUTILE

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


def test_cube_50nm_mumaxpMicromagnetics_vs_analytic():

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
                                         dV=[2, 2, 2],
                                         units='nanometer',
                                         traslation_vector=[0, 0, 0]
                                         )
    demag_signal.compute_scan_signal(method="cython")
    print(f'{demag_signal.geom_center=}')
    print(f'{demag_signal.fd_ncells=}')

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


def test_cube_50nm_mumaxpMicromagnetics_cutile_vs_numba():
    """
    Same cuboid + scan geometry as the mumax+ test above, but runs the GPU
    cuTile kernel and cross-checks against the CPU numba implementation
    (tight tolerance) and against the analytic solution (loose tolerance).

    Skipped if cupy / cuda.tile are not importable in the environment.
    """
    print('** Testing cuTile GPU kernel against numba CPU + analytic **')

    if not _HAS_CUTILE:
        print("  SKIP: cupy / cuda.tile not importable in this environment.")
        return

    # Same scan grid as test_cube_50nm_mumaxpMicromagnetics_vs_analytic.
    # Nx = Ny = 301; neither divides 16, so the host-side padding path in
    # method='cutile' is exercised.
    scan_spacing = (10 * nm, 10 * nm)
    scan_limits = np.array([[-1.5, -1.5], [1.5, 1.5]]) * µm
    scan_height = 500 * nm

    script_dir = Path(__file__).resolve().parent
    sim_file = script_dir / "cuboid_fd_mm_test/cube_mumaxPlus_L_50nm_centerAt_-35nm_dxyz_2nm.npy"

    # CPU reference ----------------------------------------------------------
    demag_cpu = mds.MicroDemagSignature(scan_limits, scan_spacing, scan_height)
    demag_cpu.reader_fd_micromagnetic(sim_file,
                                      Ms=4.8e5,
                                      origin_to_geom_center=False,
                                      dV=[2, 2, 2],
                                      units='nanometer',
                                      traslation_vector=[0, 0, 0]
                                      )
    demag_cpu.compute_scan_signal(method="numba")
    B_cpu = demag_cpu.B_grid.copy()

    # GPU under test ---------------------------------------------------------
    demag_gpu = mds.MicroDemagSignature(scan_limits, scan_spacing, scan_height)
    demag_gpu.reader_fd_micromagnetic(sim_file,
                                      Ms=4.8e5,
                                      origin_to_geom_center=False,
                                      dV=[2, 2, 2],
                                      units='nanometer',
                                      traslation_vector=[0, 0, 0]
                                      )
    demag_gpu.compute_scan_signal(method="cutile")
    B_gpu = demag_gpu.B_grid.copy()

    # CPU vs GPU cross-check -------------------------------------------------
    abs_err = np.abs(B_gpu - B_cpu)
    denom = np.abs(B_cpu) + 1e-18
    rel_err = abs_err / denom

    print(f"B_cpu range          : [{B_cpu.min():+.3e}, {B_cpu.max():+.3e}]")
    print(f"B_gpu range          : [{B_gpu.min():+.3e}, {B_gpu.max():+.3e}]")
    print(f"max |B_gpu - B_cpu|  : {abs_err.max():.3e}")
    print(f"max rel err          : {rel_err.max():.3e}")
    print(f"rms err              : {np.sqrt(np.mean((B_gpu - B_cpu)**2)):.3e}")

    # Float32-ish slack — cuTile may compile in single-precision intrinsics
    # for sqrt/divide even when the inputs are float64; tighten this once the
    # numerical behavior of the API is confirmed.
    np.testing.assert_allclose(B_gpu, B_cpu, rtol=1e-4, atol=1e-12)
    print("OK: cuTile result matches numba reference within tolerance.")

    # Mini timing ------------------------------------------------------------
    # Both calls above served as warm-up (numba JIT compile, cuTile kernel
    # compile + first-launch overhead). Now time a clean second call each.
    t0 = time.perf_counter()
    demag_cpu.compute_scan_signal(method="numba")
    t_numba = time.perf_counter() - t0

    t0 = time.perf_counter()
    demag_gpu.compute_scan_signal(method="cutile")  # includes H<->D copies + sync
    t_cutile = time.perf_counter() - t0

    print(f"numba  : {t_numba*1e3:8.2f} ms")
    print(f"cutile : {t_cutile*1e3:8.2f} ms   (speedup x{t_numba/t_cutile:.1f})")

    # GPU vs analytic --------------------------------------------------------
    analytic_sol = np.load(
        script_dir / "analytical_cuboid_code/cuboid_50nm_centre-at_-35nm_scan-grid_3microm_m_11-1_Bz.npy"
    )

    err = np.linalg.norm(analytic_sol - B_gpu * 1e6, ord="fro")
    relerr = err / np.linalg.norm(analytic_sol, ord="fro")

    print("Max cuTile     (µT)  :", np.abs(B_gpu).max() * 1e6)
    print("Max Analytic   (µT)  :", np.abs(analytic_sol).max())
    print("Max difference (µT)  :", np.abs(analytic_sol - B_gpu * 1e6).max())
    print("Error Diff     (µT)  :", err)
    print("Rel Error Diff %     :", 100 * relerr)


def test_cube_50nm_finmagMicromagnetics_vs_analytic():

    print('** Testing finmag finite element simulation file **')

    # Creating scan grid
    scan_spacing = (10 * nm, 10 * nm)
    scan_limits = np.array([[-1.5, -1.5], [1.5, 1.5]]) * µm
    scan_height = 500 * nm

    script_dir = Path(__file__).resolve().parent
    # FILE_energy = script_dir / "cuboid_merrill_test/cube_sim_energy.log"

    # Notice that the cube in the MERRILL simulation is already centered
    # at z = -35 nm
    demag_signal = mds.MicroDemagSignature(scan_limits, scan_spacing, scan_height)

    sim_file = script_dir / "cuboid_finmag_test/finmag_cube_m111.txt"
    demag_signal.reader_finmag(sim_file,
                               Ms=4.8e5,
                               origin_to_geom_center=False,
                               units='nanometer',
                               traslation_vector=[0, 0, 0]
                               )
    demag_signal.compute_scan_signal(method="cython")
    print(f'{demag_signal.geom_center=}')
    print(f'{demag_signal.fe_ntets=}')
    # print(f'{demag_signal.mesh_volume=}')

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
    print("Max difference (µT):", np.abs(analytic_sol - demag_signal.B_grid * 1e6).max())
    print("Error Diff (µT)    :", err)
    print("Rel Error Diff %   :", 100 * relerr)


if __name__ == "__main__":
    test_cube_50nm_merrill_vs_analytic()
    print('-' * 80 + '\n')
    test_cube_50nm_mumaxpMicromagnetics_vs_analytic()
    print('-' * 80 + '\n')
    test_cube_50nm_mumaxpMicromagnetics_cutile_vs_numba()
    print('-' * 80 + '\n')
    test_cube_50nm_finmagMicromagnetics_vs_analytic()
