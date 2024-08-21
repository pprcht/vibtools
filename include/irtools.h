#ifndef IRTOOLS_H
#define IRTOOLS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

// Declaration of the c_computespec_core function
extern void c_computespec_core(
    int c_nat, int *c_at, double (*c_xyz)[3], double *c_hess,
    double (*c_dipd)[3], double *c_ams, double c_fscal,
    double *c_freq, double *c_ints
);

#ifdef __cplusplus
}
#endif

#endif // IRTOOLS_H

