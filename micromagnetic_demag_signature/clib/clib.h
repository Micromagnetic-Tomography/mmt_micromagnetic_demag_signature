#ifndef CLIB_H
#define CLIB_H

void dipolar_field_C(double * dip_r, double * dip_m, unsigned long n_dip,
                     double * Sx_range, double * Sy_range,
                     unsigned long n_Sx, unsigned long n_Sy,
                     double Sheight, double * Bz_grid);

#endif
