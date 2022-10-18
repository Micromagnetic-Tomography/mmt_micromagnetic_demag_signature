cimport numpy as cnp

# -----------------------------------------------------------------------------

cdef extern from "clib.h":

    void dipolar_field_C(double * dip_r, double * dip_m, unsigned long n_dip,
                         double * Sx_range, double * Sy_range,
                         unsigned long n_Sx, unsigned long n_Sy,
                         double Sheight, double * Bz_grid)

# -----------------------------------------------------------------------------

def dipole_bz_field_C(double [:] dip_r, double [:] dip_m, unsigned long n_dip,
                      double [:] Sx_range, double [:] Sy_range,
                      unsigned long n_Sx, unsigned long n_Sy,
                      double Sheight, double [:, :] Bz_grid):

    dipolar_field_C(&dip_r[0], &dip_m[0], n_dip,
                    &Sx_range[0], &Sy_range[0],
                    n_Sx, n_Sy,
                    Sheight, &Bz_grid[0, 0])
