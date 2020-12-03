include ../config.mk

CXXFLAGS+=-I../${DOCTEST_HEADERS}

.PHONY: all clean check

all: check

check: tests
	./tests

SRC = $(wildcard test_*.cpp)
OBJ = $(subst .cpp,.o, $(SRC))

tests: tests.o $(OBJ) test-r6-relaxed.o test-r6-constrained.o

SEPARATOR=+
NAMES = $(subst test_,,$(basename $(SRC)))
TEMPLATE = $(join test_, $(1)).o: $(join test_, $(1)).cpp $(foreach h, $(subst $(SEPARATOR), ,$(1)), $(join ../include/isto/${PROJECT}/, $(h)).hpp)
$(foreach name, $(NAMES), $(eval $(call TEMPLATE,$(name))))

test-r6-relaxed.o: test-r6-relaxed.cpp ../include/isto/iapws/r6.hpp
test-r6-constrained.o: test-r6-constrained.cpp ../include/isto/iapws/r6.hpp

clean: 
	rm -f tests *.o 
