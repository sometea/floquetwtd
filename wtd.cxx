#include "wtd.h"
#include "glfixed.h"

#include <cmath>
#include <iostream>

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

WTD *objectWhichWillHandle;

double wrapperItp(double t0, void *params) {
  return objectWhichWillHandle->itp2(t0, *((double*)params));
}

WTD::WTD() {
  Nquad = 5;
  Nomega = 6;
  NT = 20;
  eta = 0.1;
  omega = 2.0*M_PI;
  Tp = 2.0*M_PI/omega;

  finalizeConstructor();
}

WTD::~WTD() {
  gsl_matrix_complex_free(preSF);
  gsl_matrix_complex_free(preK);
  delete[] x;
  delete[] w;
}

void WTD::finalizeConstructor() {
  currTau = -1.0;
  N = Nquad * Nomega;
  preSF = gsl_matrix_complex_alloc(N, Nomega);
  preK = gsl_matrix_complex_alloc(N, N);
  x = new double[N];
  w = new double[N];
  precalcQuadrature();
  precalcSF();
}

gsl_complex WTD::K(double tau, double delta) {
  return gsl_complex_mul_real(gsl_complex_exp(gsl_complex_rect(0, -delta*tau/2.0)),
                              tau*gsl_sf_sinc(delta*tau/2.0/M_PI) / 2.0 / M_PI);
}

gsl_complex WTD::SF(double e, int n) {
  double result = 0;
  if (n*omega>-e) {
    if (n==0) result = -exp(-2.0*M_PI*eta);
    else if (n>0) result =  2.0*sinh(2.0*M_PI*eta)*exp(-2.0*M_PI*n*eta);
    return gsl_complex_rect(result, 0.0);
  } else {
    return gsl_complex_rect(0.0, 0.0);
  }
}

void WTD::precalcSF() {
  for (int i=0; i<N; i++)
    for (int j=0; j<Nomega; j++)
      gsl_matrix_complex_set(preSF, i, j, SF(x[i], j));
}

void WTD::precalcQuadrature() {
  // calculate the quadrature points and weights
  gsl_integration_glfixed_table *quadTable = gsl_integration_glfixed_table_alloc(Nquad);
  for (int i=0; i<N; i++) {
    int compartment = i / Nquad;
    gsl_integration_glfixed_point(-compartment*omega, -(compartment+1)*omega,
                                  i % Nquad, &(x[i]), &(w[i]), quadTable);
  }
  gsl_integration_glfixed_table_free(quadTable);
}

void WTD::precalcK() {
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      double e1, e2;
      if (i==0) e1 = -Nomega*omega;
      else e1 = x[N-i];
      if (j==0) e2 = -Nomega*omega;
      else e2 = x[N-j];
      gsl_matrix_complex_set(preK, i, j, K(currTau, e1-e2));
    }
}

gsl_complex WTD::KfromPreK(int i, int j, int n, int m) {
  int i1 = n*Nquad-i, i2 = m*Nquad-j;
  if ((i1<0) || (i2<0)) return gsl_complex_rect(0.0,0.0);
  else return gsl_matrix_complex_get(preK, i1, i2);
}

double WTD::itp2(double t0, double tau) {
  // recalculate the K matrix if tau has changed since the last call
  if ((currTau==-1.0) || (tau!=currTau)) {
    currTau = tau;
    precalcK();
  }

  // build the matrix
  gsl_matrix_complex *M = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex_set_identity(M);
  gsl_matrix_complex *q = gsl_matrix_complex_alloc(N, N);
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      double e1 = x[i], e2 = x[j];
      gsl_complex element = gsl_complex_rect(0.0,0.0);
      for (int n=0; n<Nomega; n++)
        for (int m=0; m<Nomega; m++) {
          gsl_complex inc = gsl_complex_rect(1.0,0.0);
          inc = gsl_complex_mul(inc, gsl_complex_conjugate(gsl_matrix_complex_get(preSF, i, n)));
          inc = gsl_complex_mul(inc, gsl_matrix_complex_get(preSF, j, m));
          if ((GSL_REAL(inc)!=0.0) || (GSL_IMAG(inc)!=0.0)) {
            //inc = gsl_complex_mul(inc, K(tau, e1+n*omega-e2-m*omega));
            inc = gsl_complex_mul(inc, KfromPreK(i,j,n,m));
            inc = gsl_complex_mul(inc, gsl_complex_exp(gsl_complex_rect(0,
                    -t0*(e1+n*omega-e2-m*omega))));
          }
          element = gsl_complex_add(element, inc);
        }
      element = gsl_complex_mul_real(element, sqrt(w[i]*w[j]));
      gsl_matrix_complex_set(q, i, j, element);
    }
  gsl_matrix_complex_sub(M, q);
  gsl_matrix_complex_free(q);

  // compute the determinant
  gsl_permutation *p = gsl_permutation_alloc(N);
  int signum = 0;
  gsl_linalg_complex_LU_decomp(M, p, &signum);
  gsl_complex result = gsl_linalg_complex_LU_det(M, signum);
  gsl_permutation_free(p);

  // clean up and return the real part of the result
  gsl_matrix_complex_free(M);
  return GSL_REAL(result);
}

double WTD::itp(double tau) {
  double result = 0.0;
  for (int i=0; i<NT; i++){
    result += itp2(i*Tp/NT, tau) / NT;
  }
  return result;

  // gsl_function f;
  // double t = tau;
  // objectWhichWillHandle = this;
  // f.function = wrapperItp;
  // f.params = &t;
  // gsl_integration_glfixed_table *quadTable = gsl_integration_glfixed_table_alloc(NT);
  // double result = gsl_integration_glfixed(&f, 0, Tp, quadTable) / double(Tp);
  // gsl_integration_glfixed_table_free(quadTable);
  // return result;
}
