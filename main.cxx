#include <iostream>
#include <cmath>
#include <cstdlib>

#include "ses.h"

using namespace std;

int main(int argc, char **argv) {
  SES w(0.9, .5);
  int N = 400;
  double Tmax = 5.5;
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
