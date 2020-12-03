isto::iapws
===========

A C++20 headers-only library implementing (some) formulations released by the
`International Association for the Properties of Water and Steam <http://www.iapws.org/>`_.


Requirements
------------

Otionally, the `isto::units library <https://github.com/le-migou/units>`_. See
`Flavors`_ below.

To run the tests:

- `GNU Make <https://www.gnu.org/software/make/>`_.
- The `onqtam/doctest library <https://github.com/onqtam/doctest>`_.


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
signature for a function template is

.. code-block:: c++
        template <class T> 
        pressure_t
    pressure (density_t <T> const& density, temperature_t <T> const& temperature)

where ``pressure_t``, ``density_t`` and ``temperature_t`` exist in the namespace
``isto::units`` and wrap values with the given physical dimension.

In the *relaxed* flavor, the template functions arguments and return types can
be anything for which the computation make sense (e.g. a `double`).
A typical signature looks like

.. code-block:: c++
        template <class T> 
        T
    pressure (T const& density, T const& temperature)


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
