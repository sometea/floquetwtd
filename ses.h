#ifndef __SES_H
#define __SES_H

#include "wtd.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>

class SES: public WTD {
public:
  SES(double _Tc=0.2, double _delta=5.0);
  ~SES();
protected:
  virtual gsl_complex SF(double e, int n);
private:
  gsl_complex phase(double t);
  gsl_complex Sfrozen(double e);
  int convertIndex(int i);

  double Tc, delta, eU0;
  int NFT;
  gsl_vector_complex *c;
};

#endif
