isto::iapws
===========

A C++20 headers-only library implementing (some) guidelines released by the
`International Association for the Properties of Water and Steam <http://www.iapws.org/>`_ 
(IAPWS).

Implemented formulations are:

#. R6-95(2018) (a.k.a. IAPWS-95): `Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use <http://www.iapws.org/relguide/IAPWS-95.html>`_.
   The main, de-facto state-of-the-art, equation of state for ordinary 
   water (i.e. liquid and vapor). It is expresses its properties in terms of 
   density and temperature. Implemented in the header ``r6.hpp``.

#. R6-inverse. This formulation inverts the pressure (density, temperature)
   relation of R6 into density (pressure, temperature). It allows to express the
   properties as functions of the  pressure and temperature.
   Implemented in the header ``r6_inverse.hpp``, depends on the 
   `isto::units <https://github.com/le-migou/root_finding>`_ library.

#. R7-97(2012) (a.k.a. IF-97): `Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam <http://www.iapws.org/relguide/IF97-Rev.html>`_.
   An approximation of R6 that provides direct calculations of thermodynamic 
   properties as functions of selected independent properties (e.g. 
   pressure-temperature, enthalpy-entropy, etc.)
   This is a faster, albeit less precise, alternative to R6-inverse, designed 
   for industrial use. It is also serves as initial guess for the inversion 
   process of R6-inverse.
   Implemented in the header ``r7.hpp``.

#. R10-06(2009): `Revised Release on the Equation of State 2006 for H2O Ice Ih <http://www.iapws.org/relguide/Ice-2009.html>`_.
   Thermodynamic properties of ordinary water in its solid hexagonal phase I.
   Implemented in the header ``r10.hpp``.

#. R12-08: `Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance <http://www.iapws.org/relguide/viscosity.html>`_.
   The viscosity of water (liquid and vapor).
   Implemented in the header ``r12.hpp``.

#. R14-08(2011): `Revised Release on the Pressure along the Melting and Sublimation Curves of Ordinary Water Substance <http://www.iapws.org/relguide/MeltSub.html>`_.
   Sublimation pressure of ice Ih and melting pressure of ices Ih, III, V, VI
   and VII as functions of temperature.
   Implemented in the header ``r14.hpp``.


Dependencies
------------
 
- `isto::array <https://github.com/le-migou/array>`_.
- `isto::root_finding <https://github.com/le-migou/root_finding>`_ for R6-inverse.

Optionally, `isto::units <https://github.com/le-migou/units>`_. See
`Flavors`_ below.

To run the tests:

- `GNU Make <https://www.gnu.org/software/make/>`_.
- The `onqtam/doctest <https://github.com/onqtam/doctest>`_ library.


Installation
------------

This is a headers-only library. Just put it where your compiler can find it.


Flavors
-------

This library comes in two "flavors". The *constrained* flavor enforces physical
dimension correctness, the *relaxed* flavor does not.

If the header ``isto/units/units.hpp`` is found and the macro
``ISTO_IAPWS_FORCE_RELAXED`` is not defined when including any header of this
library, then the flavor is *constrained*.


If the header ``isto/units/units.hpp`` is not found or the macro
``ISTO_IAPWS_FORCE_RELAXED`` is defined when including any header of this
library, then the flavor is *relaxed*.

In the *constrained* flavor, the function templates arguments and return values 
have a physical dimension, as provided by the ``isto::units`` library. A typical
signature for a function template is::

        constexpr auto
    pressure (Density_t auto const& density, Temperature auto const& temperature)

where ``Density`` and ``Temperature`` are concepts defined in the namespace 
``isto::units`` that constrain the parameters into having the correct physical 
dimension. The return value of this function will satisfy the constaint
``Pressure``.

In the *relaxed* flavor, the template functions arguments and return types can
be anything for which the computation makes sense (e.g. a `double`).
A typical signature looks like::

        constexpr auto
    pressure_dt (auto const& density, auto const& temperature)

Note the ``_dt`` suffix which indicates that this function expects a density and
a temperature (in that order) as parameters. (But the compiler will not be able
to enforce it.)


Tests
-----

The tests require the `onqtam/doctest library`_.
Edit the ``config.mk`` file to make the ``DOCTEST_HEADERS`` variable point to 
the directory containing ``doctest/doctest.h``. 

To execute the tests run

    $ make check

in the root directory of the project.

The tests require the `isto::units` library to test the constrained flavor of
the library.


License
-------

SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception


Affiliation
-----------

This material is developed by the Scientific Computations and Modelling
platform at the Institut des Sciences de la Terre d'Orléans
(https://www.isto-orleans.fr/), a joint laboratory of the University of Orléans
(https://www.univ-orleans.fr/), the french National Center For Scientific
Research (https://www.cnrs.fr/) and the french Geological Survey
(https://www.brgm.eu/).
