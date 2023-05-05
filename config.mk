DOCTEST_HEADERS=../../external/onqtam/doctest
DEPENDENCIES_HEADERS=                       \
	 ../../external/isto/array/include        \
	 ../../external/isto/root_finding/include \

PROJECT=iapws
LINK.o=${LINK.cc}
CXXFLAGS+=-std=c++2a -Wall -Wextra $(foreach dir, ${DEPENDENCIES_HEADERS}, -I../${dir})
LDFLAGS+=
LDLIBS+= -lfmt

