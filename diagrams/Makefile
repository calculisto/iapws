include ../config.mk
CXXFLAGS+=-I../../../external/isto/root_finding/include    \
		  -I../../../external/isto/template_pow/include    \
		  -O3 -g

include make_vars.mk

gpls_lines = lines_tp.gpl lines_dt.gpl
gpls += $(gpls_lines)

pngs   = $(subst .gpl,.png,$(gpls))
thumbs = $(subst .gpl,-thumb.png,$(gpls))

all: $(thumbs) make_html
	./make_html

include make_rules.mk
$(gpls_lines) &: diagram_lines
	./diagram_lines

$(thumbs): %-thumb.png: %.png
	convert $< -thumbnail 320x240 $@

$(pngs): %.png: %.gpl
	gnuplot -c $<

make_vars.mk make_rules.mk &: make_make
	./make_make

make_make.o: diagram_*.hpp
make_html.o: diagram_*.hpp


clean:
	rm -fr *.dat *.gpl *.png *.txt *.mk *.o diagram_r7 diagram_r6 diagram_r6_inverse diagram_r6_inverse_extended diagram_r7_vs_r6_inverse make_html diagram_lines diagram_r10 diagram_r12 diagram_r12_pt

diagram_%.o: diagram_%.hpp diagram.hpp

