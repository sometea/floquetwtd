#include <iostream>
#include <cmath>
#include <cstdlib>

#include "ses.h"

using namespace std;

int main(int argc, char **argv) {
  // constants
  int N = 400;  // number of time steps
  double Tmax = 5.5;  // calculate ITP from 0 to Tmax

  // create the scattering object, an instance of a class derived from WTD
  SES w(0.9, .5);

  // calculate the ITP and write to stdout.
  // if one argument n is given, just calculate the n'th time step
  // otherwise, calculate all of them.
  if (argc<2) {
    for (int i=0; i<N; i++) {
      cout << Tmax/(N-1)*i << " " << w.itp(Tmax/(N-1)*i) << "\n";
    }
  } else {
    int i = atoi(argv[1]);
    cout << Tmax/(N-1)*i << " " << w.itp(Tmax/(N-1)*i) << "\n";
  }
  return 0;
}
