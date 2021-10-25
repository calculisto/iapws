include ../config.mk

CXXFLAGS+=-I../${DOCTEST_HEADERS} -g

.PHONY: all clean check

all: check

check: tests
	./tests

SRC = $(wildcard test_*.cpp)
OBJ = $(subst .cpp,.o, $(SRC))

tests: tests.o $(OBJ)             \
	test-r6-relaxed.o             \
	test-r6-constrained.o         \
	test-r7-relaxed.o             \
	test-r7-constrained.o         \
	test-r6_inverse-relaxed.o     \
	test-r6_inverse-constrained.o \
	test-r10-relaxed.o            \
	test-r10-constrained.o        \
	test-r14-relaxed.o            \
	test-r14-constrained.o        \
	test-r12-relaxed.o            \
	test-r12-constrained.o        \

SEPARATOR=+
NAMES = $(subst test_,,$(basename $(SRC)))
TEMPLATE = $(join test_, $(1)).o: $(join test_, $(1)).cpp $(foreach h, $(subst $(SEPARATOR), ,$(1)), $(join ../include/isto/${PROJECT}/, $(h)).hpp)
$(foreach name, $(NAMES), $(eval $(call TEMPLATE,$(name))))

test-r6-relaxed.o:         test-r6-relaxed.cpp         ../include/isto/iapws/r6.hpp         ../include/isto/iapws/detail/common.hpp
test-r6-constrained.o:     test-r6-constrained.cpp     ../include/isto/iapws/r6.hpp         ../include/isto/iapws/detail/common.hpp
test-r7-relaxed.o:         test-r7-relaxed.cpp         ../include/isto/iapws/r7.hpp         ../include/isto/iapws/detail/common.hpp ../include/isto/iapws/detail/data_for_the_tests.hpp
test-r7-constrained.o:     test-r7-constrained.cpp     ../include/isto/iapws/r7.hpp         ../include/isto/iapws/detail/common.hpp ../include/isto/iapws/detail/data_for_the_tests.hpp
test-r6_inverse-relaxed.o: test-r6_inverse-relaxed.cpp ../include/isto/iapws/r6_inverse.hpp ../include/isto/iapws/detail/common.hpp
test-r10-relaxed.o:        test-r10-relaxed.cpp        ../include/isto/iapws/r10.hpp        ../include/isto/iapws/detail/common.hpp
test-r10-constrained.o:    test-r10-constrained.cpp    ../include/isto/iapws/r10.hpp        ../include/isto/iapws/detail/common.hpp
test-r14-relaxed.o:        test-r14-relaxed.cpp        ../include/isto/iapws/r14.hpp        ../include/isto/iapws/detail/common.hpp
test-r14-constrained.o:    test-r14-constrained.cpp    ../include/isto/iapws/r14.hpp        ../include/isto/iapws/detail/common.hpp
test-r12-relaxed.o:        test-r12-relaxed.cpp        ../include/isto/iapws/r12.hpp        ../include/isto/iapws/detail/common.hpp
test-r12-constrained.o:    test-r12-constrained.cpp    ../include/isto/iapws/r12.hpp        ../include/isto/iapws/detail/common.hpp

clean: 
	rm -f tests *.o 
