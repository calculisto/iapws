PROJECT=iapws
LINK.o=${LINK.cc}
CXXFLAGS+=-std=c++23 -Wall -Wextra -Wconversion -Wshadow \
  -I ../git_submodules/calculisto/array/include \
  -I ../git_submodules/calculisto/root_finding/include \
  -I ../git_submodules/calculisto/auto_diff/include \
  -I ../git_submodules/calculisto/finite_difference/include \

LDLIBS+= -lfmt

