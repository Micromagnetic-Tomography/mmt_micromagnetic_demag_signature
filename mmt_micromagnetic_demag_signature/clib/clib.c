#include<stdio.h>
#include<math.h>

/* C function to compute the field from a cluster of dipole sources at a scan
 * signal surface
 *
 * Parameters
 * ----------
 * dip_r
 *     array with the dipole positions in order: x0 y0 z0 x1 y1 z1 ...
 *     The length of the array is `3 * n_dip`
 * dip_m
 *     array with dipole moments in same order than dip_r
 * n_dip
 *     number of dipole sources
 * Sx_range, Sy_range
 *     arrays with scan signal grid x and y positions
 * n_Sx, n_Sy
 *     number of grid positions in x and y directions
 * Sheight
 *     scan height
 * Bz_grid
 *     array with the scan grid positions. The size of this array is
 *     `n_Sx * n_Sy`. x direction is for columns and y direction is for rows,
 *     so it is filled by columns first (see code)
 *     Cython: in cython we pass a 2D array but it is flattened in this code
 *     assuming C order for the array entries
 *
 */

void dipolar_field_C(double * dip_r, double * dip_m, unsigned long n_dip,
                     double * Sx_range, double * Sy_range,
                     unsigned long n_Sx, unsigned long n_Sy,
                     double Sheight, double * Bz_grid) {

    // Loop through the scan positions in y and x
    // Collapsing will parallelize the two for loops
    #pragma omp parallel for collapse(2)
    for (unsigned long j = 0; j < n_Sy; ++j) {
        for (unsigned long i = 0; i < n_Sx; ++i) {

            double pos_r[3];
            pos_r[0] = Sx_range[i];
            pos_r[1] = Sy_range[j];
            pos_r[2] = Sheight;

            double r[3];
            double m_dot_r, rho, rho2;
            double bz_dip = 0;
            // For every dip position and mag moment dip_r, dip_m (Nx3 arrays)
            for (unsigned long dip_idx = 0; dip_idx < n_dip; ++dip_idx) {

                // Compute: distance from scan pos to dipole pos: r
                //          magnitude of distance               : rho
                //          dot product m * r
                m_dot_r = 0; rho = 0; rho2 = 0;
                for (int c = 0; c < 3; ++c) {
                    r[c] = pos_r[c] - dip_r[3 * dip_idx + c];
                    m_dot_r += dip_m[3 * dip_idx + c] * r[c];
                    rho2 += r[c] * r[c];
                }
                rho = sqrt(rho2);

                // Sum to the total dipole field contribution to the scan pos
                // (mu0 / 4 PI) * (3 * z * m.r / rho^5  -  m_z / rho^3)
                bz_dip += 3e-7 * r[2] * m_dot_r / (rho2 * rho2 * rho);
                bz_dip += (-1e-7) * dip_m[3 * dip_idx + 2] / (rho2 * rho);

            } // end dipole loop

            Bz_grid[n_Sx * j + i] = bz_dip;

        } // end Sx loop
    } // end Sy loop
}
