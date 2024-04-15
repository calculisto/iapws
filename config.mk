DEPENDENCIES_HEADERS=                       \
	 ../../external/calculisto/array/include        \
	 ../../external/calculisto/root_finding/include \

PROJECT=iapws
LINK.o=${LINK.cc}
CXXFLAGS+=-std=c++2a -Wall -Wextra $(foreach dir, ${DEPENDENCIES_HEADERS}, -I../${dir})
LDFLAGS+=
LDLIBS+= -lfmt

