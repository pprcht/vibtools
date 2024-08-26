#ifndef vibtools_H
#define vibtools_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

// Declaration of the c_computespec_core function
extern void c_computespec_core(int c_nat, int *c_at, double (*c_xyz)[3],
                               double *c_hess, double (*c_dipd)[3],
                               double *c_ams, double c_fscal, double *c_freq,
                               double *c_ints);

// Declaration of the print routine
void c_print_vib_spectrum_stdout(int nat3, double *freq, double *intens);

// Declaration of the broadening routine
void c_lorentzian_broadening(int c_nmodes, const double *c_freq,
                             const double *c_intens, int c_npoints,
                             double *c_plt, double c_xmin, double c_xmax,
                             double c_dx, double c_fwhm);

#ifdef __cplusplus
}
#endif

#endif // vibtools_H
