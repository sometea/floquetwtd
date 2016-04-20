# Numerical Floquet scattering calculation of electronic waiting time distributions (WTDs)

## Description
This is some code written in C++ to numerically calculate electronic waiting
time distribution (WTD) for periodically driven phase-coherent conductors. The
WTD W(t) is the probability distribution for a time t to pass between two
successive electron transfers in a mesoscopic electronic conductor.

The method to calculate WTDs from the Floquet scattering matrix is described in:
D. Dasenbrook, C. Flindt and M. Büttiker, [Floquet Theory of Electron Waiting Times in Quantum-Coherent Conductors](http://dx.doi.org/10.1103/PhysRevLett.112.146801), Physical Review Letters 112, 146801 (2014)

This code was, among other applications, used to calculate the WTDs in:
P. P. Hofer, D. Dasenbrook and C. Flindt, [Electron waiting times for the mesoscopic capacitor](http://dx.doi.org/10.1016/j.physe.2015.08.034), Physica E (2015), doi:10.1016/j.physe.2015.08.034

A different variant of this code was also used to calculate the WTDs in:
D. Dasenbrook and C. Flindt, "Quantum theory of an electron waiting time clock", [arXiv:1602.07917](arxiv.org/abs/1602.07917) (2016).

The Floquet scattering matrix formalism used by this method is valid for quantum mechanical many-electron systems, as long as interactions between the electrons can be neglected.

## Compilation
This software makes use of the [GNU Scientific Library](www.gnu.org/software/gsl) (GSL). To compile it, you need to have this library installed on your system, as well as [cmake](cmake.org). Then, compile by typing
```
cmake
```
followed by
```
make
```

It is assumed that you are on a Linux or Unix-like system, and that you have
basic build tools such as the GNU compiler collection installed. Windows/Cygwin
might work as well.

## Usage
Currently, you have to implement the scattering matrix by editing the source
code. To implement a particular Floquet scattering matrix, create a class
derived from the WTD class in the wtd.cxx and wtd.h files. As an example, have a
look at ses.cxx and ses.h, which implement the Floquet scattering matrix of the
single electron source used in the experiment by Gwendal Feve et al. in Paris.
Read the above paper in Physica E to find out more about this scattering matrix.
Then, change the main programme in main.cxx to create an instance of your new
class instead of SES. The main programme, if called with no parameters,
calculates idle time probabilities for N different time steps between 0 and Tmax
(currently hard-coded) and outputs them to stdout. If called with one parameter
n, it will just calculate the n'th time step and output it. This is useful for
launching many instances of the programme in parallel for different data points,
e.g. on a cluster.

Future versions of this software might provide a more convenient user interface.

## LICENSE
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
