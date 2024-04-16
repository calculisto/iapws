# calculisto::iapws
A C++ library implementing (some) releases and guidelines by the
International Association for the Properties of Water and Steam,
[IAPWS](http://www.iapws.org).

Implemented formulations are:

1. R6-95(2018) (a.k.a. IAPWS-95): [Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use](http://www.iapws.org/relguide/IAPWS-95.html).
   The main, de-facto state-of-the-art, equation of state for ordinary 
   water (i.e. liquid and vapor). It is expresses its properties in terms of 
   density and temperature. Implemented in the header `r6.hpp`.

2. R6-inverse. This formulation inverts the pressure (density, temperature)
   relation of R6 into density (pressure, temperature). It allows to express the
   properties as functions of the  pressure and temperature.
   Implemented in the header `r6_inverse.hpp`, depends on the 
   [calculisto::root_finding](https://github.com/calculisto/root_finding) library.

3. R7-97(2012) (a.k.a. IF-97): [Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam](http://www.iapws.org/relguide/IF97-Rev.html).
   An approximation of R6 that provides direct calculations of thermodynamic 
   properties as functions of selected independent properties (e.g. 
   pressure-temperature, enthalpy-entropy, etc.)
   This is a faster, albeit less precise, alternative to R6-inverse, designed 
   for industrial use. It is also serves as initial guess for the inversion 
   process of R6-inverse.
   Implemented in the header `r7.hpp`.

4. R10-06(2009): [Revised Release on the Equation of State 2006 for H2O Ice Ih](http://www.iapws.org/relguide/Ice-2009.html).
   Thermodynamic properties of ordinary water in its solid hexagonal phase I as
   functions of the pressure and temperature.
   Implemented in the header `r10.hpp`.

5. R10-inverse. This is R10-06 in terms of density-temperature, via a pressure /
   density inversion (similar to r6-inverse).
   Implemented in the header `r10_inverse.hpp`.

6. R12-08: [Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance](http://www.iapws.org/relguide/viscosity.html).
   The viscosity of water (liquid and vapor).
   Implemented in the header `r12.hpp`.

7. R14-08(2011): [Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance](http://www.iapws.org/relguide/MeltSub.html).
   Sublimation pressure of ice Ih and melting pressure of ices Ih, III, V, VI
   and VII as functions of temperature.
   Implemented in the header `r14.hpp`.

8. G12-15: [Guideline on Thermodynamic Properties of Supercooled Water ](http://iapws.org/relguide/Supercooled.html).
   The properties of supercooled water in terms of density-temperature.
   Implemented in the header `g12.hpp`.

9. G12-inverse. This is G12-15 in terms of pressure-temperature, via a pressure
   / density inversion.
   Implemented in the header `g12_inverse.hpp`.

Some phase diagrams built using this library are showcased here:
https://calcul-isto.cnrs-orleans.fr/misc/diagrams

## Dependencies
- [calculisto::array](https://github.com/le-migou/array),
- [calculisto::root_finding](https://github.com/le-migou/root_finding) for the
  functions that need an inversion.

To run the tests:
- A C++20 capable compiler,
- [GNU Make](https://www.gnu.org/software/make),
- The [doctest/doctest](https://github.com/onqtam/doctest) library.

## Installation
This is a headers-only library. Just put it where your compiler can find it.

## Tests
To run the tests, execute `make check` in the root directory of the project.

## License
SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

## Affiliation
This material is developed by the Numerical Modelling platform at the 
Institut des Sciences de la Terre d'Orléans (https://www.isto-orleans.fr), 
a joint laboratory of the University of Orléans (https://www.univ-orleans.fr), 
the french National Center For Scientific Research (https://www.cnrs.fr) and 
the french Geological Survey (https://www.brgm.fr).

![logo ISTO](https://calcul-isto.cnrs-orleans.fr/logos/isto-156.png) &emsp;
![logo CNRS](https://calcul-isto.cnrs-orleans.fr/logos/cnrs-128.png) &emsp;
![logo UO](https://calcul-isto.cnrs-orleans.fr/logos/uo-180.png) &emsp;
![logo BRGM](https://calcul-isto.cnrs-orleans.fr/logos/brgm-256.png)
