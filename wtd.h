#ifndef __WTD_H
#define __WTD_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

class WTD {
public:
  WTD();
  virtual ~WTD();

  // calculate the two-time ITP
  double itp2(double t0, double tau);

  // calculate the one-time ITP (integrate itp2 over one period)
  double itp(double tau);

protected:
  // The sine kernel K
  gsl_complex K(double tau, double delta);

  // get the sine kernel K from a precalulcated table to speed up the calculation
  gsl_complex KfromPreK(int i, int j, int n, int m);

  // The Floquet scattering matrix at energy e and number of exchanged energy
  // quanta n. Overwrite this method in derived classes to calculate the WTD
  // for a particular scatterer
  virtual gsl_complex SF(double e, int n);

  // calculate Gaussian quadrature points and weights used to discretize the
  // Fredholm operator whose determinant is needed
  void precalcQuadrature();

  // precalculate the Floquet scattering matrix
  void precalcSF();

  // precalculate values of the kernel function so they can be used by calling
  // KfromPreK
  void precalcK();

  // This has to be called at the end of the constructor of derived classes
  void finalizeConstructor();

  int Nquad, Nomega, NT;
  double eta, omega, Tp;
  gsl_matrix_complex *preSF;
  double *x; double *w;

  double currTau;
  gsl_matrix_complex *preK;

  int N;
};

#endif
