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
	rm -fr *.dat *.gpl *.png *.txt *.mk *.o diagram_r7 diagram_r6 diagram_r6_inverse diagram_r6_inverse_extended diagram_r7_vs_r6_inverse make_html diagram_lines

diagram_r7.o: diagram_r7.hpp diagram.hpp
diagram_r6.o: diagram_r6.hpp diagram.hpp
diagram_r6_inverse.o: diagram_r6_inverse.hpp diagram.hpp 
diagram_r6_inverse_extended.o: diagram_r6_inverse.hpp diagram.hpp
diagram_r7_vs_r6_inverse.o: diagram_r7_vs_r6_inverse.hpp diagram.hpp

#  
#  
#  
#  all: $(thumbs)
#  
#  
#  
#  $(filter %-thumb.png, $(thbs)): %.png: %.gpl
#  	gnuplot -c $<

#gpls=                                                               \
#	density_r6_inv_ylin_zlin.gpl                                    \
#	density_r6_inv_ylin_zlog.gpl                                    \
#	density_r6_inv_ylog_zlin.gpl                                    \
#	density_r6_inv_ylog_zlog.gpl                                    \
#	density_r7_vs_r6_inv_ylin_zlin.gpl                              \
#	density_r7_vs_r6_inv_ylin_zlog.gpl                              \
#	density_r7_vs_r6_inv_ylog_zlin.gpl                              \
#	density_r7_vs_r6_inv_ylog_zlog.gpl                              \
#	density_r7_ylin_zlin.gpl                                        \
#	density_r7_ylin_zlog.gpl                                        \
#	density_r7_ylog_zlin.gpl                                        \
#	density_r7_ylog_zlog.gpl                                        \
#	isobaric_cubic_expansion_coefficient_r6_inv_ylin_zlin.gpl       \
#	isobaric_cubic_expansion_coefficient_r6_inv_ylin_zlog.gpl       \
#	isobaric_cubic_expansion_coefficient_r6_inv_ylog_zlin.gpl       \
#	isobaric_cubic_expansion_coefficient_r6_inv_ylog_zlog.gpl       \
#	isobaric_cubic_expansion_coefficient_r6_ylin_zlin.gpl           \
#	isobaric_cubic_expansion_coefficient_r6_ylin_zlog.gpl           \
#	isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_ylin_zlin.gpl \
#	isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_ylin_zlog.gpl \
#	isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_ylog_zlin.gpl \
#	isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_ylog_zlog.gpl \
#	isobaric_cubic_expansion_coefficient_r7_ylin_zlin.gpl           \
#	isobaric_cubic_expansion_coefficient_r7_ylin_zlog.gpl           \
#	isobaric_cubic_expansion_coefficient_r7_ylog_zlin.gpl           \
#	isobaric_cubic_expansion_coefficient_r7_ylog_zlog.gpl           \
#	isothermal_compressibility_r6_inv_ylin_zlin.gpl                 \
#	isothermal_compressibility_r6_inv_ylin_zlog.gpl                 \
#	isothermal_compressibility_r6_inv_ylog_zlin.gpl                 \
#	isothermal_compressibility_r6_inv_ylog_zlog.gpl                 \
#	isothermal_compressibility_r6_ylin_zlin.gpl                     \
#	isothermal_compressibility_r6_ylin_zlog.gpl                     \
#	isothermal_compressibility_r7_vs_r6_inv_ylin_zlin.gpl           \
#	isothermal_compressibility_r7_vs_r6_inv_ylin_zlog.gpl           \
#	isothermal_compressibility_r7_vs_r6_inv_ylog_zlin.gpl           \
#	isothermal_compressibility_r7_vs_r6_inv_ylog_zlog.gpl           \
#	isothermal_compressibility_r7_ylin_zlin.gpl                     \
#	isothermal_compressibility_r7_ylin_zlog.gpl                     \
#	isothermal_compressibility_r7_ylog_zlin.gpl                     \
#	isothermal_compressibility_r7_ylog_zlog.gpl                     \
#	isothermal_stress_coefficient_r6_inv_ylin_zlin.gpl              \
#	isothermal_stress_coefficient_r6_inv_ylin_zlog.gpl              \
#	isothermal_stress_coefficient_r6_inv_ylog_zlin.gpl              \
#	isothermal_stress_coefficient_r6_inv_ylog_zlog.gpl              \
#	isothermal_stress_coefficient_r6_ylin_zlin.gpl                  \
#	isothermal_stress_coefficient_r6_ylin_zlog.gpl                  \
#	isothermal_stress_coefficient_r7_vs_r6_inv_ylin_zlin.gpl        \
#	isothermal_stress_coefficient_r7_vs_r6_inv_ylin_zlog.gpl        \
#	isothermal_stress_coefficient_r7_vs_r6_inv_ylog_zlin.gpl        \
#	isothermal_stress_coefficient_r7_vs_r6_inv_ylog_zlog.gpl        \
#	isothermal_stress_coefficient_r7_ylin_zlin.gpl                  \
#	isothermal_stress_coefficient_r7_ylin_zlog.gpl                  \
#	isothermal_stress_coefficient_r7_ylog_zlin.gpl                  \
#	isothermal_stress_coefficient_r7_ylog_zlog.gpl                  \
#	massic_enthalpy_r6_inv_ylin_zlin.gpl                            \
#	massic_enthalpy_r6_inv_ylin_zlog.gpl                            \
#	massic_enthalpy_r6_inv_ylog_zlin.gpl                            \
#	massic_enthalpy_r6_inv_ylog_zlog.gpl                            \
#	massic_enthalpy_r6_ylin_zlin.gpl                                \
#	massic_enthalpy_r6_ylin_zlog.gpl                                \
#	massic_enthalpy_r7_vs_r6_inv_ylin_zlin.gpl                      \
#	massic_enthalpy_r7_vs_r6_inv_ylin_zlog.gpl                      \
#	massic_enthalpy_r7_vs_r6_inv_ylog_zlin.gpl                      \
#	massic_enthalpy_r7_vs_r6_inv_ylog_zlog.gpl                      \
#	massic_enthalpy_r7_ylin_zlin.gpl                                \
#	massic_enthalpy_r7_ylin_zlog.gpl                                \
#	massic_enthalpy_r7_ylog_zlin.gpl                                \
#	massic_enthalpy_r7_ylog_zlog.gpl                                \
#	massic_entropy_r6_inv_ylin_zlin.gpl                             \
#	massic_entropy_r6_inv_ylin_zlog.gpl                             \
#	massic_entropy_r6_inv_ylog_zlin.gpl                             \
#	massic_entropy_r6_inv_ylog_zlog.gpl                             \
#	massic_entropy_r6_ylin_zlin.gpl                                 \
#	massic_entropy_r6_ylin_zlog.gpl                                 \
#	massic_entropy_r7_vs_r6_inv_ylin_zlin.gpl                       \
#	massic_entropy_r7_vs_r6_inv_ylin_zlog.gpl                       \
#	massic_entropy_r7_vs_r6_inv_ylog_zlin.gpl                       \
#	massic_entropy_r7_vs_r6_inv_ylog_zlog.gpl                       \
#	massic_entropy_r7_ylin_zlin.gpl                                 \
#	massic_entropy_r7_ylin_zlog.gpl                                 \
#	massic_entropy_r7_ylog_zlin.gpl                                 \
#	massic_entropy_r7_ylog_zlog.gpl                                 \
#	massic_gibbs_free_energy_r6_inv_ylin_zlin.gpl                   \
#	massic_gibbs_free_energy_r6_inv_ylin_zlog.gpl                   \
#	massic_gibbs_free_energy_r6_inv_ylog_zlin.gpl                   \
#	massic_gibbs_free_energy_r6_inv_ylog_zlog.gpl                   \
#	massic_gibbs_free_energy_r6_ylin_zlin.gpl                       \
#	massic_gibbs_free_energy_r6_ylin_zlog.gpl                       \
#	massic_internal_energy_r6_inv_ylin_zlin.gpl                     \
#	massic_internal_energy_r6_inv_ylin_zlog.gpl                     \
#	massic_internal_energy_r6_inv_ylog_zlin.gpl                     \
#	massic_internal_energy_r6_inv_ylog_zlog.gpl                     \
#	massic_internal_energy_r6_ylin_zlin.gpl                         \
#	massic_internal_energy_r6_ylin_zlog.gpl                         \
#	massic_internal_energy_r7_vs_r6_inv_ylin_zlin.gpl               \
#	massic_internal_energy_r7_vs_r6_inv_ylin_zlog.gpl               \
#	massic_internal_energy_r7_vs_r6_inv_ylog_zlin.gpl               \
#	massic_internal_energy_r7_vs_r6_inv_ylog_zlog.gpl               \
#	massic_internal_energy_r7_ylin_zlin.gpl                         \
#	massic_internal_energy_r7_ylin_zlog.gpl                         \
#	massic_internal_energy_r7_ylog_zlin.gpl                         \
#	massic_internal_energy_r7_ylog_zlog.gpl                         \
#	massic_isobaric_heat_capacity_r6_inv_ylin_zlin.gpl              \
#	massic_isobaric_heat_capacity_r6_inv_ylin_zlog.gpl              \
#	massic_isobaric_heat_capacity_r6_inv_ylog_zlin.gpl              \
#	massic_isobaric_heat_capacity_r6_inv_ylog_zlog.gpl              \
#	massic_isobaric_heat_capacity_r6_ylin_zlin.gpl                  \
#	massic_isobaric_heat_capacity_r6_ylin_zlog.gpl                  \
#	massic_isobaric_heat_capacity_r7_vs_r6_inv_ylin_zlin.gpl        \
#	massic_isobaric_heat_capacity_r7_vs_r6_inv_ylin_zlog.gpl        \
#	massic_isobaric_heat_capacity_r7_vs_r6_inv_ylog_zlin.gpl        \
#	massic_isobaric_heat_capacity_r7_vs_r6_inv_ylog_zlog.gpl        \
#	massic_isobaric_heat_capacity_r7_ylin_zlin.gpl                  \
#	massic_isobaric_heat_capacity_r7_ylin_zlog.gpl                  \
#	massic_isobaric_heat_capacity_r7_ylog_zlin.gpl                  \
#	massic_isobaric_heat_capacity_r7_ylog_zlog.gpl                  \
#	massic_isochoric_heat_capacity_r6_inv_ylin_zlin.gpl             \
#	massic_isochoric_heat_capacity_r6_inv_ylin_zlog.gpl             \
#	massic_isochoric_heat_capacity_r6_inv_ylog_zlin.gpl             \
#	massic_isochoric_heat_capacity_r6_inv_ylog_zlog.gpl             \
#	massic_isochoric_heat_capacity_r6_ylin_zlin.gpl                 \
#	massic_isochoric_heat_capacity_r6_ylin_zlog.gpl                 \
#	massic_isochoric_heat_capacity_r7_vs_r6_inv_ylin_zlin.gpl       \
#	massic_isochoric_heat_capacity_r7_vs_r6_inv_ylin_zlog.gpl       \
#	massic_isochoric_heat_capacity_r7_vs_r6_inv_ylog_zlin.gpl       \
#	massic_isochoric_heat_capacity_r7_vs_r6_inv_ylog_zlog.gpl       \
#	massic_isochoric_heat_capacity_r7_ylin_zlin.gpl                 \
#	massic_isochoric_heat_capacity_r7_ylin_zlog.gpl                 \
#	massic_isochoric_heat_capacity_r7_ylog_zlin.gpl                 \
#	massic_isochoric_heat_capacity_r7_ylog_zlog.gpl                 \
#	pressure_r6_ylin_zlin.gpl                                       \
#	pressure_r6_ylin_zlog.gpl                                       \
#	relative_pressure_coefficient_r6_inv_ylin_zlin.gpl              \
#	relative_pressure_coefficient_r6_inv_ylin_zlog.gpl              \
#	relative_pressure_coefficient_r6_inv_ylog_zlin.gpl              \
#	relative_pressure_coefficient_r6_inv_ylog_zlog.gpl              \
#	relative_pressure_coefficient_r6_ylin_zlin.gpl                  \
#	relative_pressure_coefficient_r6_ylin_zlog.gpl                  \
#	relative_pressure_coefficient_r7_vs_r6_inv_ylin_zlin.gpl        \
#	relative_pressure_coefficient_r7_vs_r6_inv_ylin_zlog.gpl        \
#	relative_pressure_coefficient_r7_vs_r6_inv_ylog_zlin.gpl        \
#	relative_pressure_coefficient_r7_vs_r6_inv_ylog_zlog.gpl        \
#	relative_pressure_coefficient_r7_ylin_zlin.gpl                  \
#	relative_pressure_coefficient_r7_ylin_zlog.gpl                  \
#	relative_pressure_coefficient_r7_ylog_zlin.gpl                  \
#	relative_pressure_coefficient_r7_ylog_zlog.gpl                  \
#	speed_of_sound_r6_inv_ylin_zlin.gpl                             \
#	speed_of_sound_r6_inv_ylin_zlog.gpl                             \
#	speed_of_sound_r6_inv_ylog_zlin.gpl                             \
#	speed_of_sound_r6_inv_ylog_zlog.gpl                             \
#	speed_of_sound_r6_ylin_zlin.gpl                                 \
#	speed_of_sound_r6_ylin_zlog.gpl                                 \
#	speed_of_sound_r7_vs_r6_inv_ylin_zlin.gpl                       \
#	speed_of_sound_r7_vs_r6_inv_ylin_zlog.gpl                       \
#	speed_of_sound_r7_vs_r6_inv_ylog_zlin.gpl                       \
#	speed_of_sound_r7_vs_r6_inv_ylog_zlog.gpl                       \
#	speed_of_sound_r7_ylin_zlin.gpl                                 \
#	speed_of_sound_r7_ylin_zlog.gpl                                 \
#	speed_of_sound_r7_ylog_zlin.gpl                                 \
#	speed_of_sound_r7_ylog_zlog.gpl                                 \




### dats=                                                                  \
### 	density_r6_inv_ext_ylin.dat                                        \
### 	density_r6_inv_ext_ylin_failed.dat                                 \
### 	massic_enthalpy_r6_inv_ext_ylin.dat                                \
### 	massic_gibbs_free_energy_r6_inv_ext_ylin.dat                       \
### 	massic_internal_energy_r6_inv_ext_ylin.dat                         \
### 	massic_entropy_r6_inv_ext_ylin.dat                                 \
### 	massic_isobaric_heat_capacity_r6_inv_ext_ylin.dat                  \
### 	massic_isochoric_heat_capacity_r6_inv_ext_ylin.dat                 \
### 	speed_of_sound_r6_inv_ext_ylin.dat                                 \
### 	isobaric_cubic_expansion_coefficient_r6_inv_ext_ylin.dat           \
### 	isothermal_compressibility_r6_inv_ext_ylin.dat                     \
### 	relative_pressure_coefficient_r6_inv_ext_ylin.dat                  \
### 	isothermal_stress_coefficient_r6_inv_ext_ylin.dat                  \
### 	density_r6_inv_ext_ylog.dat                                        \
### 	density_r6_inv_ext_ylog_failed.dat                                 \
### 	massic_enthalpy_r6_inv_ext_ylog.dat                                \
### 	massic_gibbs_free_energy_r6_inv_ext_ylog.dat                       \
### 	massic_internal_energy_r6_inv_ext_ylog.dat                         \
### 	massic_entropy_r6_inv_ext_ylog.dat                                 \
### 	massic_isobaric_heat_capacity_r6_inv_ext_ylog.dat                  \
### 	massic_isochoric_heat_capacity_r6_inv_ext_ylog.dat                 \
### 	speed_of_sound_r6_inv_ext_ylog.dat                                 \
### 	isobaric_cubic_expansion_coefficient_r6_inv_ext_ylog.dat           \
### 	isothermal_compressibility_r6_inv_ext_ylog.dat                     \
### 	relative_pressure_coefficient_r6_inv_ext_ylog.dat                  \
### 	isothermal_stress_coefficient_r6_inv_ext_ylog.dat                  \
### 	                                                                   \
### 	density_r6_inv_ylin.dat                                            \
### 	density_r6_inv_ylin_failed.dat                                     \
### 	massic_enthalpy_r6_inv_ylin.dat                                    \
### 	massic_gibbs_free_energy_r6_inv_ylin.dat                           \
### 	massic_internal_energy_r6_inv_ylin.dat                             \
### 	massic_entropy_r6_inv_ylin.dat                                     \
### 	massic_isobaric_heat_capacity_r6_inv_ylin.dat                      \
### 	massic_isochoric_heat_capacity_r6_inv_ylin.dat                     \
### 	speed_of_sound_r6_inv_ylin.dat                                     \
### 	isobaric_cubic_expansion_coefficient_r6_inv_ylin.dat               \
### 	isothermal_compressibility_r6_inv_ylin.dat                         \
### 	relative_pressure_coefficient_r6_inv_ylin.dat                      \
### 	isothermal_stress_coefficient_r6_inv_ylin.dat                      \
### 	density_r6_inv_ylog.dat                                            \
### 	density_r6_inv_ylog_failed.dat                                     \
### 	massic_enthalpy_r6_inv_ylog.dat                                    \
### 	massic_gibbs_free_energy_r6_inv_ylog.dat                           \
### 	massic_internal_energy_r6_inv_ylog.dat                             \
### 	massic_entropy_r6_inv_ylog.dat                                     \
### 	massic_isobaric_heat_capacity_r6_inv_ylog.dat                      \
### 	massic_isochoric_heat_capacity_r6_inv_ylog.dat                     \
### 	speed_of_sound_r6_inv_ylog.dat                                     \
### 	isobaric_cubic_expansion_coefficient_r6_inv_ylog.dat               \
### 	isothermal_compressibility_r6_inv_ylog.dat                         \
### 	relative_pressure_coefficient_r6_inv_ylog.dat                      \
### 	isothermal_stress_coefficient_r6_inv_ylog.dat                      \
### 	                                                                   \
### 	density_r7_ylin.dat                                                \
### 	massic_enthalpy_r7_ylin.dat                                        \
### 	massic_internal_energy_r7_ylin.dat                                 \
### 	massic_entropy_r7_ylin.dat                                         \
### 	massic_isobaric_heat_capacity_r7_ylin.dat                          \
### 	massic_isochoric_heat_capacity_r7_ylin.dat                         \
### 	speed_of_sound_r7_ylin.dat                                         \
### 	isobaric_cubic_expansion_coefficient_r7_ylin.dat                   \
### 	isothermal_compressibility_r7_ylin.dat                             \
### 	relative_pressure_coefficient_r7_ylin.dat                          \
### 	isothermal_stress_coefficient_r7_ylin.dat                          \
### 	density_r7_ylog.dat                                                \
### 	massic_enthalpy_r7_ylog.dat                                        \
### 	massic_internal_energy_r7_ylog.dat                                 \
### 	massic_entropy_r7_ylog.dat                                         \
### 	massic_isobaric_heat_capacity_r7_ylog.dat                          \
### 	massic_isochoric_heat_capacity_r7_ylog.dat                         \
### 	speed_of_sound_r7_ylog.dat                                         \
### 	isobaric_cubic_expansion_coefficient_r7_ylog.dat                   \
### 	isothermal_compressibility_r7_ylog.dat                             \
### 	relative_pressure_coefficient_r7_ylog.dat                          \
### 	isothermal_stress_coefficient_r7_ylog.dat                          \
### 	                                                                   \
### 	density_r7_vs_r6_inv_error_ylin.dat                                \
### 	density_r7_vs_r6_inv_iter_ylin.dat                                 \
### 	massic_enthalpy_r7_vs_r6_inv_error_ylin.dat                        \
### 	massic_internal_energy_r7_vs_r6_inv_error_ylin.dat                 \
### 	massic_entropy_r7_vs_r6_inv_error_ylin.dat                         \
### 	massic_isobaric_heat_capacity_r7_vs_r6_inv_error_ylin.dat          \
### 	massic_isochoric_heat_capacity_r7_vs_r6_inv_error_ylin.dat         \
### 	speed_of_sound_r7_vs_r6_inv_error_ylin.dat                         \
### 	isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error_ylin.dat   \
### 	isothermal_compressibility_r7_vs_r6_inv_error_ylin.dat             \
### 	relative_pressure_coefficient_r7_vs_r6_inv_error_ylin.dat          \
### 	isothermal_stress_coefficient_r7_vs_r6_inv_error_ylin.dat          \
### 	density_r7_vs_r6_inv_error_ylog.dat                                \
### 	density_r7_vs_r6_inv_iter_ylog.dat                                 \
### 	massic_enthalpy_r7_vs_r6_inv_error_ylog.dat                        \
### 	massic_internal_energy_r7_vs_r6_inv_error_ylog.dat                 \
### 	massic_entropy_r7_vs_r6_inv_error_ylog.dat                         \
### 	massic_isobaric_heat_capacity_r7_vs_r6_inv_error_ylog.dat          \
### 	massic_isochoric_heat_capacity_r7_vs_r6_inv_error_ylog.dat         \
### 	speed_of_sound_r7_vs_r6_inv_error_ylog.dat                         \
### 	isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error_ylog.dat   \
### 	isothermal_compressibility_r7_vs_r6_inv_error_ylog.dat             \
### 	relative_pressure_coefficient_r7_vs_r6_inv_error_ylog.dat          \
### 	isothermal_stress_coefficient_r7_vs_r6_inv_error_ylog.dat          \
### 	                                                                   \
### 	saturation.dat                                                     \
### 	sublimation.dat                                                    \
### 	melting_ih.dat                                                     \
### 	melting_iii.dat                                                    \
### 	melting_v.dat                                                      \
### 	melting_vi.dat                                                     \
### 	melting_vii.dat                                                    \
### 	                                                                   \
###     sublimation_dt.dat                                                 \
###     saturation_gaz_dt.dat                                              \
###     saturation_liq_dt.dat                                              \
###     melting_ih_dt.dat                                                  \
###     melting_iii_dt.dat                                                 \
###     melting_v_dt.dat                                                   \
###     melting_vi_dt.dat                                                  \
###     melting_vii_dt.dat                                                 \
### 	                                                                   \
### 	isobaric_cubic_expansion_coefficient_r6_ylin.dat                   \
### 	isothermal_compressibility_r6_ylin.dat                             \
### 	isothermal_stress_coefficient_r6_ylin.dat                          \
### 	massic_enthalpy_r6_ylin.dat                                        \
### 	massic_entropy_r6_ylin.dat                                         \
### 	massic_gibbs_free_energy_r6_ylin.dat                               \
### 	massic_internal_energy_r6_ylin.dat                                 \
### 	massic_isobaric_heat_capacity_r6_ylin.dat                          \
### 	massic_isochoric_heat_capacity_r6_ylin.dat                         \
### 	pressure_r6_ylin.dat                                               \
### 	relative_pressure_coefficient_r6_ylin.dat                          \
### 	speed_of_sound_r6_ylin.dat                                         \
### 
### pngs =                                                                    \
###     density_r6_inv_ext_ylin_zlin.png                                      \
###     density_r6_inv_ext_ylin_zlog.png                                      \
###     density_r6_inv_ext_ylog_zlin.png                                      \
###     density_r6_inv_ext_ylog_zlog.png                                      \
###     massic_enthalpy_r6_inv_ext_ylin_zlin.png                              \
###     massic_enthalpy_r6_inv_ext_ylin_zlog.png                              \
###     massic_enthalpy_r6_inv_ext_ylog_zlin.png                              \
###     massic_enthalpy_r6_inv_ext_ylog_zlog.png                              \
###     massic_gibbs_free_energy_r6_inv_ext_ylin_zlin.png                     \
###     massic_gibbs_free_energy_r6_inv_ext_ylin_zlog.png                     \
###     massic_gibbs_free_energy_r6_inv_ext_ylog_zlin.png                     \
###     massic_gibbs_free_energy_r6_inv_ext_ylog_zlog.png                     \
###     massic_internal_energy_r6_inv_ext_ylin_zlin.png                       \
###     massic_internal_energy_r6_inv_ext_ylin_zlog.png                       \
###     massic_internal_energy_r6_inv_ext_ylog_zlin.png                       \
###     massic_internal_energy_r6_inv_ext_ylog_zlog.png                       \
###     massic_entropy_r6_inv_ext_ylin_zlin.png                               \
###     massic_entropy_r6_inv_ext_ylin_zlog.png                               \
###     massic_entropy_r6_inv_ext_ylog_zlin.png                               \
###     massic_entropy_r6_inv_ext_ylog_zlog.png                               \
###     massic_isobaric_heat_capacity_r6_inv_ext_ylin_zlin.png                \
###     massic_isobaric_heat_capacity_r6_inv_ext_ylin_zlog.png                \
###     massic_isobaric_heat_capacity_r6_inv_ext_ylog_zlin.png                \
###     massic_isobaric_heat_capacity_r6_inv_ext_ylog_zlog.png                \
###     massic_isochoric_heat_capacity_r6_inv_ext_ylin_zlin.png               \
###     massic_isochoric_heat_capacity_r6_inv_ext_ylin_zlog.png               \
###     massic_isochoric_heat_capacity_r6_inv_ext_ylog_zlin.png               \
###     massic_isochoric_heat_capacity_r6_inv_ext_ylog_zlog.png               \
###     speed_of_sound_r6_inv_ext_ylin_zlin.png                               \
###     speed_of_sound_r6_inv_ext_ylin_zlog.png                               \
###     speed_of_sound_r6_inv_ext_ylog_zlin.png                               \
###     speed_of_sound_r6_inv_ext_ylog_zlog.png                               \
###     isobaric_cubic_expansion_coefficient_r6_inv_ext_ylin_zlin.png         \
###     isobaric_cubic_expansion_coefficient_r6_inv_ext_ylin_zlog.png         \
###     isobaric_cubic_expansion_coefficient_r6_inv_ext_ylog_zlin.png         \
###     isobaric_cubic_expansion_coefficient_r6_inv_ext_ylog_zlog.png         \
###     isothermal_compressibility_r6_inv_ext_ylin_zlin.png                   \
###     isothermal_compressibility_r6_inv_ext_ylin_zlog.png                   \
###     isothermal_compressibility_r6_inv_ext_ylog_zlin.png                   \
###     isothermal_compressibility_r6_inv_ext_ylog_zlog.png                   \
###     relative_pressure_coefficient_r6_inv_ext_ylin_zlin.png                \
###     relative_pressure_coefficient_r6_inv_ext_ylin_zlog.png                \
###     relative_pressure_coefficient_r6_inv_ext_ylog_zlin.png                \
###     relative_pressure_coefficient_r6_inv_ext_ylog_zlog.png                \
###     isothermal_stress_coefficient_r6_inv_ext_ylin_zlin.png                \
###     isothermal_stress_coefficient_r6_inv_ext_ylin_zlog.png                \
###     isothermal_stress_coefficient_r6_inv_ext_ylog_zlin.png                \
###     isothermal_stress_coefficient_r6_inv_ext_ylog_zlog.png                \
###     density_r6_inv_ylin_zlin.png                                          \
###     density_r6_inv_ylin_zlog.png                                          \
###     density_r6_inv_ylog_zlin.png                                          \
###     density_r6_inv_ylog_zlog.png                                          \
###     massic_enthalpy_r6_inv_ylin_zlin.png                                  \
###     massic_enthalpy_r6_inv_ylin_zlog.png                                  \
###     massic_enthalpy_r6_inv_ylog_zlin.png                                  \
###     massic_enthalpy_r6_inv_ylog_zlog.png                                  \
###     massic_gibbs_free_energy_r6_inv_ylin_zlin.png                         \
###     massic_gibbs_free_energy_r6_inv_ylin_zlog.png                         \
###     massic_gibbs_free_energy_r6_inv_ylog_zlin.png                         \
###     massic_gibbs_free_energy_r6_inv_ylog_zlog.png                         \
###     massic_internal_energy_r6_inv_ylin_zlin.png                           \
###     massic_internal_energy_r6_inv_ylin_zlog.png                           \
###     massic_internal_energy_r6_inv_ylog_zlin.png                           \
###     massic_internal_energy_r6_inv_ylog_zlog.png                           \
###     massic_entropy_r6_inv_ylin_zlin.png                                   \
###     massic_entropy_r6_inv_ylin_zlog.png                                   \
###     massic_entropy_r6_inv_ylog_zlin.png                                   \
###     massic_entropy_r6_inv_ylog_zlog.png                                   \
###     massic_isobaric_heat_capacity_r6_inv_ylin_zlin.png                    \
###     massic_isobaric_heat_capacity_r6_inv_ylin_zlog.png                    \
###     massic_isobaric_heat_capacity_r6_inv_ylog_zlin.png                    \
###     massic_isobaric_heat_capacity_r6_inv_ylog_zlog.png                    \
###     massic_isochoric_heat_capacity_r6_inv_ylin_zlin.png                   \
###     massic_isochoric_heat_capacity_r6_inv_ylin_zlog.png                   \
###     massic_isochoric_heat_capacity_r6_inv_ylog_zlin.png                   \
###     massic_isochoric_heat_capacity_r6_inv_ylog_zlog.png                   \
###     speed_of_sound_r6_inv_ylin_zlin.png                                   \
###     speed_of_sound_r6_inv_ylin_zlog.png                                   \
###     speed_of_sound_r6_inv_ylog_zlin.png                                   \
###     speed_of_sound_r6_inv_ylog_zlog.png                                   \
###     isobaric_cubic_expansion_coefficient_r6_inv_ylin_zlin.png             \
###     isobaric_cubic_expansion_coefficient_r6_inv_ylin_zlog.png             \
###     isobaric_cubic_expansion_coefficient_r6_inv_ylog_zlin.png             \
###     isobaric_cubic_expansion_coefficient_r6_inv_ylog_zlog.png             \
###     isothermal_compressibility_r6_inv_ylin_zlin.png                       \
###     isothermal_compressibility_r6_inv_ylin_zlog.png                       \
###     isothermal_compressibility_r6_inv_ylog_zlin.png                       \
###     isothermal_compressibility_r6_inv_ylog_zlog.png                       \
###     relative_pressure_coefficient_r6_inv_ylin_zlin.png                    \
###     relative_pressure_coefficient_r6_inv_ylin_zlog.png                    \
###     relative_pressure_coefficient_r6_inv_ylog_zlin.png                    \
###     relative_pressure_coefficient_r6_inv_ylog_zlog.png                    \
###     isothermal_stress_coefficient_r6_inv_ylin_zlin.png                    \
###     isothermal_stress_coefficient_r6_inv_ylin_zlog.png                    \
###     isothermal_stress_coefficient_r6_inv_ylog_zlin.png                    \
###     isothermal_stress_coefficient_r6_inv_ylog_zlog.png                    \
###     density_r7_ylin_zlin.png                                              \
###     density_r7_ylin_zlog.png                                              \
###     density_r7_ylog_zlin.png                                              \
###     density_r7_ylog_zlog.png                                              \
###     massic_enthalpy_r7_ylin_zlin.png                                      \
###     massic_enthalpy_r7_ylin_zlog.png                                      \
###     massic_enthalpy_r7_ylog_zlin.png                                      \
###     massic_enthalpy_r7_ylog_zlog.png                                      \
###     massic_internal_energy_r7_ylin_zlin.png                               \
###     massic_internal_energy_r7_ylin_zlog.png                               \
###     massic_internal_energy_r7_ylog_zlin.png                               \
###     massic_internal_energy_r7_ylog_zlog.png                               \
###     massic_entropy_r7_ylin_zlin.png                                       \
###     massic_entropy_r7_ylin_zlog.png                                       \
###     massic_entropy_r7_ylog_zlin.png                                       \
###     massic_entropy_r7_ylog_zlog.png                                       \
###     massic_isobaric_heat_capacity_r7_ylin_zlin.png                        \
###     massic_isobaric_heat_capacity_r7_ylin_zlog.png                        \
###     massic_isobaric_heat_capacity_r7_ylog_zlin.png                        \
###     massic_isobaric_heat_capacity_r7_ylog_zlog.png                        \
###     massic_isochoric_heat_capacity_r7_ylin_zlin.png                       \
###     massic_isochoric_heat_capacity_r7_ylin_zlog.png                       \
###     massic_isochoric_heat_capacity_r7_ylog_zlin.png                       \
###     massic_isochoric_heat_capacity_r7_ylog_zlog.png                       \
###     speed_of_sound_r7_ylin_zlin.png                                       \
###     speed_of_sound_r7_ylin_zlog.png                                       \
###     speed_of_sound_r7_ylog_zlin.png                                       \
###     speed_of_sound_r7_ylog_zlog.png                                       \
###     isobaric_cubic_expansion_coefficient_r7_ylin_zlin.png                 \
###     isobaric_cubic_expansion_coefficient_r7_ylin_zlog.png                 \
###     isobaric_cubic_expansion_coefficient_r7_ylog_zlin.png                 \
###     isobaric_cubic_expansion_coefficient_r7_ylog_zlog.png                 \
###     isothermal_compressibility_r7_ylin_zlin.png                           \
###     isothermal_compressibility_r7_ylin_zlog.png                           \
###     isothermal_compressibility_r7_ylog_zlin.png                           \
###     isothermal_compressibility_r7_ylog_zlog.png                           \
###     relative_pressure_coefficient_r7_ylin_zlin.png                        \
###     relative_pressure_coefficient_r7_ylin_zlog.png                        \
###     relative_pressure_coefficient_r7_ylog_zlin.png                        \
###     relative_pressure_coefficient_r7_ylog_zlog.png                        \
###     isothermal_stress_coefficient_r7_ylin_zlin.png                        \
###     isothermal_stress_coefficient_r7_ylin_zlog.png                        \
###     isothermal_stress_coefficient_r7_ylog_zlin.png                        \
###     isothermal_stress_coefficient_r7_ylog_zlog.png                        \
###     density_r7_vs_r6_inv_error_ylin_zlin.png                              \
###     density_r7_vs_r6_inv_error_ylin_zlog.png                              \
###     density_r7_vs_r6_inv_error_ylog_zlin.png                              \
###     density_r7_vs_r6_inv_error_ylog_zlog.png                              \
###     density_r7_vs_r6_inv_iter_ylin_zlin.png                               \
###     density_r7_vs_r6_inv_iter_ylin_zlog.png                               \
###     density_r7_vs_r6_inv_iter_ylog_zlin.png                               \
###     density_r7_vs_r6_inv_iter_ylog_zlog.png                               \
###     massic_enthalpy_r7_vs_r6_inv_error_ylin_zlin.png                      \
###     massic_enthalpy_r7_vs_r6_inv_error_ylin_zlog.png                      \
###     massic_enthalpy_r7_vs_r6_inv_error_ylog_zlin.png                      \
###     massic_enthalpy_r7_vs_r6_inv_error_ylog_zlog.png                      \
###     massic_internal_energy_r7_vs_r6_inv_error_ylin_zlin.png               \
###     massic_internal_energy_r7_vs_r6_inv_error_ylin_zlog.png               \
###     massic_internal_energy_r7_vs_r6_inv_error_ylog_zlin.png               \
###     massic_internal_energy_r7_vs_r6_inv_error_ylog_zlog.png               \
###     massic_entropy_r7_vs_r6_inv_error_ylin_zlin.png                       \
###     massic_entropy_r7_vs_r6_inv_error_ylin_zlog.png                       \
###     massic_entropy_r7_vs_r6_inv_error_ylog_zlin.png                       \
###     massic_entropy_r7_vs_r6_inv_error_ylog_zlog.png                       \
###     massic_isobaric_heat_capacity_r7_vs_r6_inv_error_ylin_zlin.png        \
###     massic_isobaric_heat_capacity_r7_vs_r6_inv_error_ylin_zlog.png        \
###     massic_isobaric_heat_capacity_r7_vs_r6_inv_error_ylog_zlin.png        \
###     massic_isobaric_heat_capacity_r7_vs_r6_inv_error_ylog_zlog.png        \
###     massic_isochoric_heat_capacity_r7_vs_r6_inv_error_ylin_zlin.png       \
###     massic_isochoric_heat_capacity_r7_vs_r6_inv_error_ylin_zlog.png       \
###     massic_isochoric_heat_capacity_r7_vs_r6_inv_error_ylog_zlin.png       \
###     massic_isochoric_heat_capacity_r7_vs_r6_inv_error_ylog_zlog.png       \
###     speed_of_sound_r7_vs_r6_inv_error_ylin_zlin.png                       \
###     speed_of_sound_r7_vs_r6_inv_error_ylin_zlog.png                       \
###     speed_of_sound_r7_vs_r6_inv_error_ylog_zlin.png                       \
###     speed_of_sound_r7_vs_r6_inv_error_ylog_zlog.png                       \
###     isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error_ylin_zlin.png \
###     isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error_ylin_zlog.png \
###     isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error_ylog_zlin.png \
###     isobaric_cubic_expansion_coefficient_r7_vs_r6_inv_error_ylog_zlog.png \
###     isothermal_compressibility_r7_vs_r6_inv_error_ylin_zlin.png           \
###     isothermal_compressibility_r7_vs_r6_inv_error_ylin_zlog.png           \
###     isothermal_compressibility_r7_vs_r6_inv_error_ylog_zlin.png           \
###     isothermal_compressibility_r7_vs_r6_inv_error_ylog_zlog.png           \
###     relative_pressure_coefficient_r7_vs_r6_inv_error_ylin_zlin.png        \
###     relative_pressure_coefficient_r7_vs_r6_inv_error_ylin_zlog.png        \
###     relative_pressure_coefficient_r7_vs_r6_inv_error_ylog_zlin.png        \
###     relative_pressure_coefficient_r7_vs_r6_inv_error_ylog_zlog.png        \
###     isothermal_stress_coefficient_r7_vs_r6_inv_error_ylin_zlin.png        \
###     isothermal_stress_coefficient_r7_vs_r6_inv_error_ylin_zlog.png        \
###     isothermal_stress_coefficient_r7_vs_r6_inv_error_ylog_zlin.png        \
###     isothermal_stress_coefficient_r7_vs_r6_inv_error_ylog_zlog.png        \
###                                                                           \
### 	isobaric_cubic_expansion_coefficient_r6_ylin_zlin.png                 \
### 	isothermal_compressibility_r6_ylin_zlin.png                           \
### 	isothermal_stress_coefficient_r6_ylin_zlin.png                        \
### 	massic_enthalpy_r6_ylin_zlin.png                                      \
### 	massic_entropy_r6_ylin_zlin.png                                       \
### 	massic_gibbs_free_energy_r6_ylin_zlin.png                             \
### 	massic_internal_energy_r6_ylin_zlin.png                               \
### 	massic_isobaric_heat_capacity_r6_ylin_zlin.png                        \
### 	massic_isochoric_heat_capacity_r6_ylin_zlin.png                       \
### 	pressure_r6_ylin_zlin.png                                             \
### 	relative_pressure_coefficient_r6_ylin_zlin.png                        \
### 	speed_of_sound_r6_ylin_zlin.png                                       \
### 	isobaric_cubic_expansion_coefficient_r6_ylin_zlog.png                 \
### 	isothermal_compressibility_r6_ylin_zlog.png                           \
### 	isothermal_stress_coefficient_r6_ylin_zlog.png                        \
### 	massic_enthalpy_r6_ylin_zlog.png                                      \
### 	massic_entropy_r6_ylin_zlog.png                                       \
### 	massic_gibbs_free_energy_r6_ylin_zlog.png                             \
### 	massic_internal_energy_r6_ylin_zlog.png                               \
### 	massic_isobaric_heat_capacity_r6_ylin_zlog.png                        \
### 	massic_isochoric_heat_capacity_r6_ylin_zlog.png                       \
### 	pressure_r6_ylin_zlog.png                                             \
### 	relative_pressure_coefficient_r6_ylin_zlog.png                        \
### 	speed_of_sound_r6_ylin_zlog.png                                       \
### 
### gpls = $(subst .png,.gpl,$(pngs))
### thbs = $(subst .png,-thumb.png,$(pngs)) lines-thumb.png lines_dt-thumb.png
### 
### all: $(thbs)
### 
### lines.png: lines.gpl    \
### 		melting_ih.dat  \
### 		melting_iii.dat \
### 		melting_v.dat   \
### 		melting_vi.dat  \
### 		melting_vii.dat \
### 		saturation.dat  \
### 		sublimation.dat 
### 	gnuplot -c lines.gpl
### 
### lines_dt.png: lines_dt.gpl    \
### 		melting_ih_dt.dat     \
### 		melting_iii_dt.dat    \
### 		melting_v_dt.dat      \
### 		melting_vi_dt.dat     \
### 		melting_vii_dt.dat    \
### 		saturation_liq_dt.dat \
### 		saturation_gaz_dt.dat \
### 		sublimation_dt.dat 
### 	gnuplot -c lines_dt.gpl
### 
### $(filter %_ylin_zlin.png, $(pngs)): %_ylin_zlin.png: %_ylin.dat %_ylin_zlin.gpl
### 	gnuplot -c $*_ylin_zlin.gpl
### $(filter %_ylin_zlog.png, $(pngs)): %_ylin_zlog.png: %_ylin.dat %_ylin_zlog.gpl
### 	gnuplot -c $*_ylin_zlog.gpl
### $(filter %_ylog_zlin.png, $(pngs)): %_ylog_zlin.png: %_ylog.dat %_ylog_zlin.gpl
### 	gnuplot -c $*_ylog_zlin.gpl
### $(filter %_ylog_zlog.png, $(pngs)): %_ylog_zlog.png: %_ylog.dat %_ylog_zlog.gpl
### 	gnuplot -c $*_ylog_zlog.gpl
### 
### $(gpls) lines.gpl lines_dt.gpl &: make_gpl
### 	./make_gpl
### 
### $(dats) &: make_data
### 	./make_data
### 
### $(filter %-thumb.png, $(thbs)): %-thumb.png: %.png
### 	convert $< -thumbnail 320x240 $@
### 
### #thumbs: $(pngs) lines.png lines_dt.png
### #	parallel-moreutils -j 8 -i sh -c 'convert {} -thumbnail 320x240 $$(echo {} | sed s/.png/-thumb.png/)' -- $(pngs) lines.png lines_dt.png
### 
### clean:
### 	rm -f *.png *.gpl
### 
### cleaner: clean
### 	rm -f *.dat *.txt
