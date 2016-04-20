#ifndef __SES_H
#define __SES_H

#include "wtd.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>

// Implement the Floquet scattering matrix of the mesoscopic capacitor

class SES: public WTD {
public:
  SES(double _Tc=0.2, double _delta=5.0);
  ~SES();

protected:
  // The Floquet S matrix
  virtual gsl_complex SF(double e, int n);

private:
  // the frozen scattering phase at time t
  gsl_complex phase(double t);

  // the frozen scattering matrix at energy e
  gsl_complex Sfrozen(double e);

  // convert a frequency index that can be smaller than zero to an index
  // betwen 0 and NFT
  int convertIndex(int i);

  double Tc, delta, eU0;
  int NFT;
  gsl_vector_complex *c;
};

#endif
