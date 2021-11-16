FMTLIB_HEADERS=../../external/fmtlib/fmt/include
FMTLIB_LIBRARY=../../external/fmtlib/fmt/build
DOCTEST_HEADERS=../../external/onqtam/doctest
DEPENDENCIES_HEADERS=                         \
	 ../../external/isto/array/include        \
	 ../../external/isto/units/include        \
	 ../../external/isto/root_finding/include \

PROJECT=iapws
LINK.o=${LINK.cc}
CXXFLAGS+=-std=c++2a -Wall -Wextra -I../${FMTLIB_HEADERS} $(foreach dir, ${DEPENDENCIES_HEADERS}, -I../${dir})
LDFLAGS+= -L../${FMTLIB_LIBRARY}
LDLIBS+= -lfmt

