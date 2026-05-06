PROJECT=iapws
LINK.o=${LINK.cc}
CXXFLAGS+=-std=c++20 -Wall -Wextra \
					-I ../git_submodules/calculisto/array/include \
	 			  -I ../git_submodules/calculisto/root_finding/include \

LDLIBS+= -lfmt

