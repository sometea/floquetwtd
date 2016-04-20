#include "ses.h"

#include <cmath>
#include <iostream>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_fft_complex.h>

SES::SES(double _Tc, double _delta) {
  // constants
  omega = 2.0*M_PI;  // angular frequency of the driving
  Tp = 2.0*M_PI/omega;  // the period of the driving

  Tc = _Tc;  // time duration of the lifiting of the top gate potential
  delta = _delta * omega;  // energy spacing of the quantum dot
  eU0 = -delta / 2.0;  // initial position of the energy level

  Nquad = 5;  // number of quadrature points for the Fredholm operator
  Nomega = int(2.5*delta/omega);  // number of Floquet scattering matrix elements
  NT = 80;  // number of points for the time integration
  N = Nquad * Nomega;
  NFT = 4 * Nomega;

  gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(NFT);
  gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(NFT);

  // calculate the scattering phases
  c = gsl_vector_complex_alloc(NFT);
  for (int i=0; i<NFT; i++)
    gsl_vector_complex_set(c, i, gsl_complex_div_real(phase(Tp/(NFT-1)*i), double(NFT)));

  // do a Fourier transform
  gsl_fft_complex_forward(c->data, c->stride, c->size, wavetable, workspace);

  gsl_fft_complex_wavetable_free(wavetable);
  gsl_fft_complex_workspace_free(workspace);

  // precalculate Floquet scattering matrix elements
  finalizeConstructor();
}

SES::~SES() {
  gsl_vector_complex_free(c);
}

gsl_complex SES::phase(double t) {
  t = t - floor(t);
  if (t<0.5)
    return gsl_complex_exp(gsl_complex_rect(0.0, -eU0*t));
  else
    return gsl_complex_exp(gsl_complex_rect(0.0, (-eU0-delta)*(t-1.0)));
}

gsl_complex SES::Sfrozen(double e) {
  gsl_complex p = gsl_complex_exp(gsl_complex_rect(0.0, e*2.0*M_PI/delta));
  gsl_complex mp = gsl_complex_negative(p);
  double r = sqrt(1.0-Tc);
  gsl_complex mrp = gsl_complex_mul_real(p, -r);
  return gsl_complex_div(gsl_complex_add_real(mp, r), gsl_complex_add_real(mrp, 1.0));
}

int SES::convertIndex(int i) {
  if (i>=0) return i;
  else return NFT+i;
}

gsl_complex SES::SF(double e, int n) {
  if (n*omega>-e) {
    gsl_complex result = gsl_complex_rect(0.0, 0.0);
    for (int m=-Nomega; m<=Nomega; m++) {
      gsl_complex inc = gsl_complex_rect(1.0, 0.0);
      inc = gsl_complex_mul(inc, gsl_vector_complex_get(c, convertIndex(m)));
      inc = gsl_complex_mul(inc, gsl_complex_conjugate(
                                  gsl_vector_complex_get(c, convertIndex(n+m))));
      inc = gsl_complex_mul(inc, Sfrozen(e-m*omega));
      result = gsl_complex_add(result, inc);
    }
    return result;
  } else {
    return gsl_complex_rect(0.0, 0.0);
  }
}
