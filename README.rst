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
   `isto::root_finding <https://github.com/le-migou/root_finding>`_ library.

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

The full list of functions is summarized in `Functions summary`_.

Some phase diagrams built with this library are showcased in this page 
https://calcul-isto.cnrs-orleans.fr/misc/diagrams/

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


Functions summary
-----------------

The table below summarizes the functions exposed by the library.

========================================== ================= ===============
Function                                   Argument          Argument   
========================================== ================= ===============
r6::pressure                               density           temperature
r6::massic_internal_energy                 density           temperature
r6::massic_entropy                         density           temperature
r6::massic_enthalpy                        density           temperature
r6::massic_isochoric_heat_capacity         density           temperature
r6::massic_isobaric_heat_capacity          density           temperature
r6::massic_gibbs_free_energy               density           temperature
r6::speed_of_sound                         density           temperature
r6::isothermal_stress_coefficient          density           temperature
r6::relative_pressure_coefficient          density           temperature
r6_inverse::density                        pressure          temperature
r7::massic_volume                          pressure          temperature
r7::density                                pressure          temperature
r7::massic_enthalpy                        pressure          temperature
r7::massic_internal_energy                 pressure          temperature
r7::massic_entropy                         pressure          temperature
r7::massic_isobaric_heat_capacity          pressure          temperature
r7::massic_isochoric_heat_capacity         pressure          temperature
r7::speed_of_sound                         pressure          temperature
r7::isobaric_cubic_expansion_coefficient   pressure          temperature
r7::isothermal_compressibility             pressure          temperature
r7::relative_pressure_coefficient          pressure          temperature
r7::isothermal_stress_coefficient          pressure          temperature
r7::pressure                               massic_enthalpy   massic_entropy
r7::temperature                            massic_enthalpy   massic_entropy
r7::massic_volume                          massic_enthalpy   massic_entropy
r7::density                                massic_enthalpy   massic_entropy
r7::massic_internal_energy                 massic_enthalpy   massic_entropy
r7::massic_isobaric_heat_capacity          massic_enthalpy   massic_entropy
r7::speed_of_sound                         massic_enthalpy   massic_entropy
r7::temperature                            pressure          massic_enthalpy
r7::density                                pressure          massic_enthalpy
r7::massic_internal_energy                 pressure          massic_enthalpy
r7::massic_entropy                         pressure          massic_enthalpy
r7::massic_isobaric_heat_capacity          pressure          massic_enthalpy
r7::speed_of_sound                         pressure          massic_enthalpy
r7::massic_volume                          pressure          massic_enthalpy
r7::temperature                            pressure          massic_entropy
r7::density                                pressure          massic_entropy
r7::massic_internal_energy                 pressure          massic_entropy
r7::massic_enthalpy                        pressure          massic_entropy
r7::massic_isobaric_heat_capacity          pressure          massic_entropy
r7::speed_of_sound                         pressure          massic_entropy
r7::massic_volume                          pressure          massic_entropy
r10::massic_volume                         pressure          temperature
r10::density                               pressure          temperature
r10::massic_entropy                        pressure          temperature
r10::massic_isobaric_heat_capacity         pressure          temperature
r10::massic_enthalpy                       pressure          temperature
r10::massic_internal_energy                pressure          temperature
r10::massic_helmholtz_energy               pressure          temperature
r10::cubic_expansion_coefficient           pressure          temperature
r10::pressure_coefficient                  pressure          temperature
r10::isothermal_compressibility            pressure          temperature
r10::isentropic_compressibility            pressure          temperature
r12::viscosity                             temperature       density
r14::ih::melting_pressure                  temperature
r14::ih::sublimation_pressure              temperature
r14::iii::melting_pressure                 temperature
r14::v::melting_pressure                   temperature
r14::vi::melting_pressure                  temperature
r14::vii::melting_pressure                 temperature
========================================== ================= ===============


Tests
-----

The tests require the `onqtam/doctest <https://github.com/onqtam/doctest>`_ 
testing framework.

Edit the ``config.mk`` file to make the ``DOCTEST_HEADERS`` variable point to 
the directory containing ``doctest/doctest.h``. 

To execute the tests run

    $ make check

in the root directory of the project.


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

+-------------+-----------+-------------+-------------+
| |logo_isto| | |logo_uo| | |logo_cnrs| | |logo_brgm| |
+-------------+-----------+-------------+-------------+

.. |logo_isto| image:: https://calcul-isto.cnrs-orleans.fr/logos/isto-156.png
   :width: 156px
   :target: https://www.isto-orleans.fr/
   :align: middle
.. |logo_uo| image:: https://calcul-isto.cnrs-orleans.fr/logos/uo-180.png
   :width: 180px
   :target: https://www.univ-orleans.fr/
   :align: middle
.. |logo_cnrs| image:: https://calcul-isto.cnrs-orleans.fr/logos/cnrs-128.png
   :width: 128px
   :target: https://www.cnrs.fr/
   :align: middle
.. |logo_brgm| image:: https://calcul-isto.cnrs-orleans.fr/logos/brgm-256.png
   :width: 256px
   :target: https://www.brgm.fr/
   :align: middle
