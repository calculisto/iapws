include ../config.mk
CXXFLAGS+=-I../../../external/isto/root_finding/include \
		  -I../../../external/isto/uncertain_value/include \
		  -I../../../external/isto/template_pow/include  \
		  -I../../../external/onqtam/doctest

all: density_r7.png          \
	density_r7_plog.png      \
	density_r7_dlog.png      \
	density_r7_plog_dlog.png \
	lines.png


density_r7.png: density_r7.dat density_r7.gpl
	gnuplot -c density_r7.gpl
density_r7_plog.png: density_r7_plog.dat density_r7_plog.gpl
	gnuplot -c density_r7_plog.gpl
density_r7_dlog.png: density_r7.dat density_r7_dlog.gpl
	gnuplot -c density_r7_dlog.gpl
density_r7_plog_dlog.png: density_r7_plog.dat density_r7_plog_dlog.gpl
	gnuplot -c density_r7_plog_dlog.gpl
lines.png: lines.gpl        \
		melting_ih.dat  \
		melting_iii.dat \
		melting_v.dat   \
		melting_vi.dat  \
		melting_vii.dat \
		saturation.dat  \
		sublimation.dat 
	gnuplot -c lines.gpl

%.dat: diagram
	./diagram

diagram: diagram.cpp

clean:
	rm -f *.png *.dat
