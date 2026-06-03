import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# import json
# import collections
import mmt_micromagnetic_demag_signature as mds
from mmt_micromagnetic_demag_signature.microds import _HAS_CUPY, numba_cuda

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


def test_cube_50nm_mumaxpMicromagnetics_gpu_benchmark():
    """
    Same cuboid + scan geometry as the mumax+ test above. For each available
    GPU backend (`numba_cuda`, `cupy`), cross-check against the numba CPU
    reference (tight tolerance) and print a wall-time vs. numba CPU.

    Each backend is skipped individually if its runtime isn't available.
    """
    print('** GPU backend benchmark vs numba CPU reference **')

    scan_spacing = (10 * nm, 10 * nm)
    scan_limits = np.array([[-1.5, -1.5], [1.5, 1.5]]) * µm
    scan_height = 500 * nm

    script_dir = Path(__file__).resolve().parent
    sim_file = script_dir / "cuboid_fd_mm_test/cube_mumaxPlus_L_50nm_centerAt_-35nm_dxyz_2nm.npy"

    def load_signal():
        s = mds.MicroDemagSignature(scan_limits, scan_spacing, scan_height)
        s.reader_fd_micromagnetic(sim_file,
                                  Ms=4.8e5,
                                  origin_to_geom_center=False,
                                  dV=[2, 2, 2],
                                  units='nanometer',
                                  traslation_vector=[0, 0, 0])
        return s

    # CPU reference (also warms up the numba JIT) ----------------------------
    demag_cpu = load_signal()
    demag_cpu.compute_scan_signal(method="numba")
    B_cpu = demag_cpu.B_grid.copy()

    t0 = time.perf_counter()
    demag_cpu.compute_scan_signal(method="numba")
    t_numba = time.perf_counter() - t0
    print(f"numba CPU       : {t_numba*1e3:8.2f} ms   "
          f"(B range [{B_cpu.min():+.3e}, {B_cpu.max():+.3e}])")

    # ---- numba.cuda --------------------------------------------------------
    if not numba_cuda.is_available():
        print("numba_cuda      : SKIP (no CUDA-capable GPU detected)")
    else:
        demag = load_signal()
        demag.compute_scan_signal(method="numba_cuda")   # warm-up
        B_ncu = demag.B_grid.copy()

        np.testing.assert_allclose(B_ncu, B_cpu, rtol=1e-6, atol=1e-15)

        t0 = time.perf_counter()
        demag.compute_scan_signal(method="numba_cuda")
        t_ncu = time.perf_counter() - t0
        print(f"numba_cuda      : {t_ncu*1e3:8.2f} ms   "
              f"(x{t_numba/t_ncu:5.1f} vs numba CPU)   "
              f"max |Δ|={np.abs(B_ncu - B_cpu).max():.2e}")

    # ---- cupy ElementwiseKernel -------------------------------------------
    if not _HAS_CUPY:
        print("cupy            : SKIP (cupy not importable)")
    else:
        demag = load_signal()
        demag.compute_scan_signal(method="cupy")         # warm-up (NVRTC compile)
        B_cup = demag.B_grid.copy()

        np.testing.assert_allclose(B_cup, B_cpu, rtol=1e-6, atol=1e-15)

        t0 = time.perf_counter()
        demag.compute_scan_signal(method="cupy")
        t_cup = time.perf_counter() - t0
        print(f"cupy            : {t_cup*1e3:8.2f} ms   "
              f"(x{t_numba/t_cup:5.1f} vs numba CPU)   "
              f"max |Δ|={np.abs(B_cup - B_cpu).max():.2e}")


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
    test_cube_50nm_mumaxpMicromagnetics_gpu_benchmark()
    print('-' * 80 + '\n')
    test_cube_50nm_finmagMicromagnetics_vs_analytic()
