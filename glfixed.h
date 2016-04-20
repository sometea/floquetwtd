#ifndef __GLFIXED_H
#define __GLFIXED_H

#include <stdlib.h>
#include <gsl/gsl_math.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* Workspace for fixed-order Gauss-Legendre integration */

typedef struct
  {
    size_t n;         /* number of points */
    double *x;        /* Gauss abscissae/points */
    double *w;        /* Gauss weights for each abscissae */
    int precomputed;  /* high precision abscissae/weights precomputed? */
  }
gsl_integration_glfixed_table;


gsl_integration_glfixed_table *
  gsl_integration_glfixed_table_alloc (size_t n);

void
  gsl_integration_glfixed_table_free (gsl_integration_glfixed_table * t);

/* Routine for fixed-order Gauss-Legendre integration */

double
  gsl_integration_glfixed (const gsl_function *f,
                           double a,
                           double b,
                           const gsl_integration_glfixed_table * t);

/* Routine to retrieve the i-th Gauss-Legendre point and weight from t */

int
  gsl_integration_glfixed_point (double a,
                                 double b,
                                 size_t i,
                                 double *xi,
                                 double *wi,
                                 const gsl_integration_glfixed_table * t);

__END_DECLS

#endif
