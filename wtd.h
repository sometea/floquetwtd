#ifndef __WTD_H
#define __WTD_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

class WTD {
public:
  WTD();
  virtual ~WTD();
  double itp2(double t0, double tau);
  double itp(double tau);
protected:
  gsl_complex K(double tau, double delta);
  gsl_complex KfromPreK(int i, int j, int n, int m);
  virtual gsl_complex SF(double e, int n);
  void precalcQuadrature();
  void precalcSF();
  void precalcK();

  void finalizeConstructor();  // has to be called at the end of the constructor of derived classes

  int Nquad, Nomega, NT;
  double eta, omega, Tp;
  gsl_matrix_complex *preSF;
  double *x; double *w;

  double currTau;
  gsl_matrix_complex *preK;

  int N;
};

#endif
