"""
FP64 verification + timing harness for the optimized GPU kernels.

Drop this next to `test_cube_50nm_side_sim-vs-analytic.py` (it reuses the same
data layout: `cuboid_fd_mm_test/` and `analytical_cuboid_code/`). It loads the
mumax+ 50 nm cube, runs every available backend, and reports for each:

  * relerr vs the analytic Bz (Frobenius)  -> should stay < 1e-4
  * max |Δ| and max relative Δ vs the numba CPU reference -> should be ~1e-12
    or better for FP64 (summation-order differences only)
  * wall time, with device buffers already cached/warmed so the number
    reflects kernel + D2H, not one-off uploads or JIT/NVRTC compile.

Run:  python verify_gpu_fp64.py
"""
import time
from pathlib import Path

import numpy as np

import mmt_micromagnetic_demag_signature as mds
from mmt_micromagnetic_demag_signature.microds import _HAS_CUPY, numba_cuda

nm = 1e-9
um = 1e-6

# Same data layout as the test file. Edit if you keep the harness elsewhere.
DATA_DIR = Path(__file__).resolve().parent
SIM_FILE = DATA_DIR / "cuboid_fd_mm_test/cube_mumaxPlus_L_50nm_centerAt_-35nm_dxyz_2nm.npy"
ANALYTIC = DATA_DIR / "analytical_cuboid_code/cuboid_50nm_centre-at_-35nm_scan-grid_3microm_m_11-1_Bz.npy"

SCAN_SPACING = (10 * nm, 10 * nm)
SCAN_LIMITS = np.array([[-1.5, -1.5], [1.5, 1.5]]) * um
SCAN_HEIGHT = 500 * nm


def load_signal():
    s = mds.MicroDemagSignature(SCAN_LIMITS, SCAN_SPACING, SCAN_HEIGHT)
    s.reader_fd_micromagnetic(SIM_FILE,
                              Ms=4.8e5,
                              origin_to_geom_center=False,
                              dV=[2, 2, 2],
                              units='nanometer',
                              traslation_vector=[0, 0, 0])
    return s


def frob_relerr(B_uT, analytic_uT):
    num = np.linalg.norm(analytic_uT - B_uT, ord="fro")
    return num / np.linalg.norm(analytic_uT, ord="fro")


def vs_reference(B, B_ref):
    absd = np.abs(B - B_ref)
    scale = np.maximum(np.abs(B_ref), 1e-300)
    return absd.max(), (absd / scale).max()


def time_method(demag, method, repeats=3):
    """Warm up (builds the cached device buffers + JIT/NVRTC), then time."""
    demag.compute_scan_signal(method=method)          # warm-up, not timed
    B = demag.B_grid.copy()
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        demag.compute_scan_signal(method=method)
        best = min(best, time.perf_counter() - t0)
    return B, best


def main():
    analytic = np.load(ANALYTIC)                      # already in uT

    # ---- numba CPU reference ------------------------------------------------
    demag_cpu = load_signal()
    B_cpu, t_cpu = time_method(demag_cpu, "numba")
    print(f"Loaded {demag_cpu.r.shape[0]} dipoles, "
          f"scan grid {demag_cpu.B_grid.shape}\n")

    rows = []
    rows.append(("numba (CPU ref)", B_cpu, t_cpu))

    # ---- numba.cuda ---------------------------------------------------------
    if not numba_cuda.is_available():
        print("numba_cuda : SKIP (no CUDA-capable GPU)")
    else:
        demag = load_signal()
        B_ncu, t_ncu = time_method(demag, "numba_cuda")
        rows.append(("numba_cuda", B_ncu, t_ncu))

    # ---- cupy ---------------------------------------------------------------
    if not _HAS_CUPY:
        print("cupy       : SKIP (cupy not importable)")
    else:
        demag = load_signal()
        B_cup, t_cup = time_method(demag, "cupy")
        rows.append(("cupy", B_cup, t_cup))

    # ---- report -------------------------------------------------------------
    hdr = f"{'backend':<16}{'relerr_analytic':>18}{'maxΔ_vs_CPU':>16}" \
          f"{'maxrelΔ_vs_CPU':>18}{'time_ms':>12}{'speedup':>10}"
    print(hdr)
    print("-" * len(hdr))
    for name, B, t in rows:
        relerr = frob_relerr(B * 1e6, analytic)
        maxd, maxrel = vs_reference(B, B_cpu)
        speed = t_cpu / t
        print(f"{name:<16}{relerr:>18.3e}{maxd:>16.3e}"
              f"{maxrel:>18.3e}{t*1e3:>12.2f}{speed:>9.1f}x")

    print("\nExpected: relerr_analytic < 1e-4 for every backend; "
          "maxrelΔ_vs_CPU ~1e-12 or smaller (FP64 summation-order noise).")


if __name__ == "__main__":
    main()
