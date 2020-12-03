.PHONY: all check clean

all: check

check:
	${MAKE} -C tests check
clean:
	${MAKE} -C tests clean
