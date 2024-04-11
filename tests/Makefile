include ../config.mk

CXXFLAGS+=-g

.PHONY: all clean check

all: check

check: tests
	./tests

SRC = $(wildcard test_*.cpp)
OBJ = $(subst .cpp,.o, $(SRC))

tests: tests.o $(OBJ) \
	test-r6.o           \
	test-r7.o           \
	test-r6_inverse.o   \
	test-r10.o          \
	test-r14.o          \
	test-r12.o          \
	test-g12.o          \
	test-r10_inverse.o  \
	test-g12_inverse.o  \

SEPARATOR=+
NAMES = $(subst test_,,$(basename $(SRC)))
TEMPLATE = $(join test_, $(1)).o: $(join test_, $(1)).cpp $(foreach h, $(subst $(SEPARATOR), ,$(1)), $(join ../include/isto/${PROJECT}/, $(h)).hpp)
$(foreach name, $(NAMES), $(eval $(call TEMPLATE,$(name))))

test-r6.o:          test-r6.cpp          ../include/isto/iapws/r6.hpp          ../include/isto/iapws/detail/common.hpp
test-r7.o:          test-r7.cpp          ../include/isto/iapws/r7.hpp          ../include/isto/iapws/detail/common.hpp ../include/isto/iapws/detail/data_for_the_tests.hpp
test-r6_inverse.o:  test-r6_inverse.cpp  ../include/isto/iapws/r6_inverse.hpp  ../include/isto/iapws/detail/common.hpp
test-r10.o:         test-r10.cpp         ../include/isto/iapws/r10.hpp         ../include/isto/iapws/detail/common.hpp
test-r14.o:         test-r14.cpp         ../include/isto/iapws/r14.hpp         ../include/isto/iapws/detail/common.hpp
test-r12.o:         test-r12.cpp         ../include/isto/iapws/r12.hpp         ../include/isto/iapws/detail/common.hpp
test-g12.o:         test-g12.cpp         ../include/isto/iapws/g12.hpp         ../include/isto/iapws/detail/common.hpp
test-r10_inverse.o: test-r10_inverse.cpp ../include/isto/iapws/r10_inverse.hpp ../include/isto/iapws/detail/common.hpp
test-g12_inverse.o: test-g12_inverse.cpp ../include/isto/iapws/g12_inverse.hpp ../include/isto/iapws/detail/common.hpp

clean: 
	rm -f tests *.o 
