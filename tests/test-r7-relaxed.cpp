#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r7.hpp"
    using namespace isto::iapws::r7;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"

    constexpr auto
eps = 1e-8;

TEST_CASE("r7.hpp (relaxed)")
{
    SUBCASE("saturation line")
    {
        CHECK(saturation_pressure_t (300.) == Approx { 0.353658941e-2 * 1e6 }.scale (1e4).epsilon (eps));
        CHECK(saturation_pressure_t (500.) == Approx { 0.263889776e1 * 1e6  }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure_t (600.) == Approx { 0.123443146e2 * 1e6  }.scale (1e8).epsilon (eps));
        CHECK(saturation_temperature_p (0.1 * 1e6) == Approx { 0.372755919e3 }.scale (1e3).epsilon (eps));
        CHECK(saturation_temperature_p (1.  * 1e6) == Approx { 0.453035632e3 }.scale (1e3).epsilon (eps));
        CHECK(saturation_temperature_p (10. * 1e6) == Approx { 0.584149488e3 }.scale (1e3).epsilon (eps));
        CHECK(saturation_pressure_h (1700e3) == Approx { 1.724175718e7 }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure_h (2000e3) == Approx { 2.193442957e7 }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure_h (2400e3) == Approx { 2.018090839e7 }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure_s (3.8e3) == Approx { 1.687755057e7 }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure_s (4.2e3) == Approx { 2.164451789e7 }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure_s (5.2e3) == Approx { 1.668968482e7 }.scale (1e7).epsilon (eps));
        CHECK(saturation_temperature_hs (1800e3, 5.3e3) == Approx { 3.468475498e2 }.scale (1e2).epsilon (eps));
        CHECK(saturation_temperature_hs (2400e3, 6.0e3) == Approx { 4.251373305e2 }.scale (1e2).epsilon (eps));
        CHECK(saturation_temperature_hs (2500e3, 5.5e3) == Approx { 5.225579013e2 }.scale (1e2).epsilon (eps));
    }
    SUBCASE("iapws-r7-region-1")
    {
        CHECK(r1::massic_volume_pt                 (3. * 1e6, 300.) == Approx { 0.100215168e-2 }.scale (1e-2).epsilon (eps));
        CHECK(r1::massic_volume_pt                 (80.* 1e6, 300.) == Approx { 0.971180894e-3 }.scale (1e-2).epsilon (eps));
        CHECK(r1::massic_volume_pt                 (3. * 1e6, 500.) == Approx { 0.120241800e-2 }.scale (1e-2).epsilon (eps));
        CHECK(r1::massic_enthalpy_pt               (3. * 1e6, 300.) == Approx { 0.115331273e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(r1::massic_enthalpy_pt               (80.* 1e6, 300.) == Approx { 0.184142828e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(r1::massic_enthalpy_pt               (3. * 1e6, 500.) == Approx { 0.975542239e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(r1::massic_internal_energy_pt        (3. * 1e6, 300.) == Approx { 0.112324818e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(r1::massic_internal_energy_pt        (80.* 1e6, 300.) == Approx { 0.106448356e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(r1::massic_internal_energy_pt        (3. * 1e6, 500.) == Approx { 0.971934985e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(r1::massic_entropy_pt                (3. * 1e6, 300.) == Approx { 0.392294792   * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(r1::massic_entropy_pt                (80.* 1e6, 300.) == Approx { 0.368563852   * 1e3 }.scale (1e3).epsilon (eps));
        CHECK(r1::massic_entropy_pt                (3. * 1e6, 500.) == Approx { 0.258041912e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r1::massic_isobaric_heat_capacity_pt (3. * 1e6, 300.) == Approx { 0.417301218e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r1::massic_isobaric_heat_capacity_pt (80.* 1e6, 300.) == Approx { 0.401008987e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r1::massic_isobaric_heat_capacity_pt (3. * 1e6, 500.) == Approx { 0.465580682e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r1::speed_of_sound_pt                (3. * 1e6, 300.) == Approx { 0.150773921e4  }.scale (1e4).epsilon (eps));
        CHECK(r1::speed_of_sound_pt                (80.* 1e6, 300.) == Approx { 0.163469054e4  }.scale (1e4).epsilon (eps));
        CHECK(r1::speed_of_sound_pt                (3. * 1e6, 500.) == Approx { 0.124071337e4  }.scale (1e4).epsilon (eps));
        CHECK(r1::temperature_ph ( 3e6,  500e3) == Approx { 0.391798509e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1::temperature_ph (80e6,  500e3) == Approx { 0.378108626e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1::temperature_ph (80e6, 1500e3) == Approx { 0.611041229e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1::temperature_ps ( 3e6,  0.5e3) == Approx { 0.307842258e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1::temperature_ps (80e6,  0.5e3) == Approx { 0.309979785e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1::temperature_ps (80e6,  3.0e3) == Approx { 0.565899909e3 }.scale (1e3).epsilon (eps));
        CHECK(r1::pressure_hs (0.001e3,    0.) == Approx { 9.800980612e-4 }.scale (1e4).epsilon (eps));
        CHECK(r1::pressure_hs (   90e3,    0.) == Approx { 9.192954727e1  }.scale (1e1).epsilon (eps));
        CHECK(r1::pressure_hs ( 1500e3, 3.4e3) == Approx { 5.868294423e1  }.scale (1e1).epsilon (eps));
        // temperature_hs ?
    };
    SUBCASE("iapws-r7-region-2")
    {
        CHECK(r2::massic_volume_pt                 (0.0035e6, 300.) == Approx { 0.394913866e2        }.scale (1e2).epsilon (eps));
        CHECK(r2::massic_volume_pt                 (0.0035e6, 700.) == Approx { 0.923015898e2        }.scale (1e2).epsilon (eps));
        CHECK(r2::massic_volume_pt                 (    30e6, 700.) == Approx { 0.542946619e-2       }.scale (1e-2).epsilon (eps));
        CHECK(r2::massic_enthalpy_pt               (0.0035e6, 300.) == Approx { 0.254991145e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2::massic_enthalpy_pt               (0.0035e6, 700.) == Approx { 0.333568375e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2::massic_enthalpy_pt               (    30e6, 700.) == Approx { 0.263149474e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2::massic_internal_energy_pt        (0.0035e6, 300.) == Approx { 0.241169160e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2::massic_internal_energy_pt        (0.0035e6, 700.) == Approx { 0.301262819e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2::massic_internal_energy_pt        (    30e6, 700.) == Approx { 0.246861076e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2::massic_entropy_pt                (0.0035e6, 300.) == Approx { 0.852238967e1  * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(r2::massic_entropy_pt                (0.0035e6, 700.) == Approx { 0.101749996e2  * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(r2::massic_entropy_pt                (    30e6, 700.) == Approx { 0.517540298e1  * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(r2::massic_isobaric_heat_capacity_pt (0.0035e6, 300.) == Approx { 0.191300162e1  * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(r2::massic_isobaric_heat_capacity_pt (0.0035e6, 700.) == Approx { 0.208141274e1  * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(r2::massic_isobaric_heat_capacity_pt (    30e6, 700.) == Approx { 0.103505092e2  * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(r2::speed_of_sound_pt                (0.0035e6, 300.) == Approx { 0.427920172e3        }.scale (1e3).epsilon (eps));
        CHECK(r2::speed_of_sound_pt                (0.0035e6, 700.) == Approx { 0.644289068e3        }.scale (1e3).epsilon (eps));
        CHECK(r2::speed_of_sound_pt                (    30e6, 700.) == Approx { 0.480386523e3        }.scale (1e3).epsilon (eps));

        CHECK(r2::b2bc_p (0.100000000e9)  == Approx { 0.3516004323e7 }.scale (1e7).epsilon (eps));
        CHECK(r2::b2bc_h (0.3516004323e7) == Approx { 0.100000000e9  }.scale (1e9).epsilon (eps));

        CHECK(r2::subregion_ph (0.001e6, 3000e3) == 1);
        CHECK(r2::subregion_ph (    3e6, 3000e3) == 1);
        CHECK(r2::subregion_ph (    3e6, 4000e3) == 1);
        CHECK(r2::subregion_ph (    5e6, 3500e3) == 2);
        CHECK(r2::subregion_ph (    5e6, 4000e3) == 2);
        CHECK(r2::subregion_ph (   25e6, 3500e3) == 2);
        CHECK(r2::subregion_ph (   40e6, 2700e3) == 3);
        CHECK(r2::subregion_ph (   60e6, 2700e3) == 3);
        CHECK(r2::subregion_ph (   60e6, 3200e3) == 3);

        CHECK(r2::subregion_ps (  0.1e6,  7.5e3) == 1);
        CHECK(r2::subregion_ps (  0.1e6,    8e3) == 1);
        CHECK(r2::subregion_ps (  2.5e6,    8e3) == 1);
        CHECK(r2::subregion_ps (    8e6,    6e3) == 2);
        CHECK(r2::subregion_ps (    8e6,  7.5e3) == 2);
        CHECK(r2::subregion_ps (   90e6,    6e3) == 2);
        CHECK(r2::subregion_ps (   20e6, 5.75e3) == 3);
        CHECK(r2::subregion_ps (   80e6, 5.25e3) == 3);
        CHECK(r2::subregion_ps (   80e6, 5.75e3) == 3);

        CHECK(r2::temperature_ph (0.001e6, 3000e3) == Approx { 0.534433241e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ph (    3e6, 3000e3) == Approx { 0.575373370e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ph (    3e6, 4000e3) == Approx { 0.101077577e4 }.scale (1e4).epsilon (eps));
        CHECK(r2::temperature_ph (    5e6, 3500e3) == Approx { 0.801299102e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ph (    5e6, 4000e3) == Approx { 0.101531583e4 }.scale (1e4).epsilon (eps));
        CHECK(r2::temperature_ph (   25e6, 3500e3) == Approx { 0.875279054e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ph (   40e6, 2700e3) == Approx { 0.743056411e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ph (   60e6, 2700e3) == Approx { 0.791137067e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ph (   60e6, 3200e3) == Approx { 0.882756860e3 }.scale (1e3).epsilon (eps));

        CHECK(r2::temperature_ps (  0.1e6,  7.5e3) == Approx { 0.399517097e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ps (  0.1e6,    8e3) == Approx { 0.514127081e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ps (  2.5e6,    8e3) == Approx { 0.103984917e4 }.scale (1e4).epsilon (eps));
        CHECK(r2::temperature_ps (    8e6,    6e3) == Approx { 0.600484040e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ps (    8e6,  7.5e3) == Approx { 0.106495556e4 }.scale (1e4).epsilon (eps));
        CHECK(r2::temperature_ps (   90e6,    6e3) == Approx { 0.103801126e4 }.scale (1e4).epsilon (eps));
        CHECK(r2::temperature_ps (   20e6, 5.75e3) == Approx { 0.697992849e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ps (   80e6, 5.25e3) == Approx { 0.854011484e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::temperature_ps (   80e6, 5.75e3) == Approx { 0.949017998e3 }.scale (1e3).epsilon (eps));

        CHECK(r2::b2ab_s (7e3) == Approx { 3376.437884e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::subregion_hs (2800e3, 6.5e3) == 1);
        CHECK(r2::subregion_hs (2800e3, 9.5e3) == 1);
        CHECK(r2::subregion_hs (4100e3, 9.5e3) == 1);
        CHECK(r2::subregion_hs (2800e3, 6.0e3) == 2);
        CHECK(r2::subregion_hs (3600e3, 6.0e3) == 2);
        CHECK(r2::subregion_hs (3600e3, 7.0e3) == 2);
        CHECK(r2::subregion_hs (2800e3, 5.1e3) == 3);
        CHECK(r2::subregion_hs (2800e3, 5.8e3) == 3);
        CHECK(r2::subregion_hs (3400e3, 5.8e3) == 3);
        CHECK(r2::pressure_hs (2800e3, 6.5e3) == Approx { 1.371012767e6 }.scale (1e6).epsilon (eps));
        CHECK(r2::pressure_hs (2800e3, 9.5e3) == Approx { 1.879743844e3 }.scale (1e3).epsilon (eps));
        CHECK(r2::pressure_hs (4100e3, 9.5e3) == Approx { 1.024788997e5 }.scale (1e5).epsilon (eps));
        CHECK(r2::pressure_hs (2800e3, 6.0e3) == Approx { 4.793911442e6 }.scale (1e6).epsilon (eps));
        CHECK(r2::pressure_hs (3600e3, 6.0e3) == Approx { 8.395519209e7 }.scale (1e7).epsilon (eps));
        CHECK(r2::pressure_hs (3600e3, 7.0e3) == Approx { 7.527161441e6 }.scale (1e6).epsilon (eps));
        CHECK(r2::pressure_hs (2800e3, 5.1e3) == Approx { 9.439202060e7 }.scale (1e7).epsilon (eps));
        CHECK(r2::pressure_hs (2800e3, 5.8e3) == Approx { 8.414574124e6 }.scale (1e6).epsilon (eps));
        CHECK(r2::pressure_hs (3400e3, 5.8e3) == Approx { 8.376903879e7 }.scale (1e7).epsilon (eps));
        // temperature_hs ?
    };
    SUBCASE("iapws-r7-region-2-metastable-vapor")
    {
        CHECK(r2_metastable_vapor::massic_volume_pt                 (1.0e6, 450.) == Approx { 0.192516540         }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_volume_pt                 (1.0e6, 440.) == Approx { 0.186212297         }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_volume_pt                 (1.5e6, 450.) == Approx { 0.121685206         }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_enthalpy_pt               (1.0e6, 450.) == Approx { 0.276881115e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_enthalpy_pt               (1.0e6, 440.) == Approx { 0.274015123e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_enthalpy_pt               (1.5e6, 450.) == Approx { 0.272134539e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_internal_energy_pt        (1.0e6, 450.) == Approx { 0.257629461e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_internal_energy_pt        (1.0e6, 440.) == Approx { 0.255393894e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_internal_energy_pt        (1.5e6, 450.) == Approx { 0.253881758e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_entropy_pt                (1.0e6, 450.) == Approx { 0.656660377e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_entropy_pt                (1.0e6, 440.) == Approx { 0.650218759e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_entropy_pt                (1.5e6, 450.) == Approx { 0.629170440e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_isobaric_heat_capacity_pt (1.0e6, 450.) == Approx { 0.276349265e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_isobaric_heat_capacity_pt (1.0e6, 440.) == Approx { 0.298166443e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::massic_isobaric_heat_capacity_pt (1.5e6, 450.) == Approx { 0.362795578e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::speed_of_sound_pt                (1.0e6, 450.) == Approx { 0.498408101e3       }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::speed_of_sound_pt                (1.0e6, 440.) == Approx { 0.489363295e3       }.scale (1e4).epsilon (eps));
        CHECK(r2_metastable_vapor::speed_of_sound_pt                (1.5e6, 450.) == Approx { 0.481941819e3       }.scale (1e4).epsilon (eps));
    };
    SUBCASE("iapws-r7-region-3")
    {
        CHECK(r3::pressure_dt                      (500., 650.) == Approx { 0.255837018e2 * 1e6 }.scale (1e2).epsilon (eps));
        CHECK(r3::pressure_dt                      (200., 650.) == Approx { 0.222930643e2 * 1e6 }.scale (1e2).epsilon (eps));
        CHECK(r3::pressure_dt                      (500., 750.) == Approx { 0.783095639e2 * 1e6 }.scale (1e2).epsilon (eps));
        CHECK(r3::massic_enthalpy_dt               (500., 650.) == Approx { 0.186343019e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r3::massic_enthalpy_dt               (200., 650.) == Approx { 0.237512401e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r3::massic_enthalpy_dt               (500., 750.) == Approx { 0.225868845e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r3::massic_internal_energy_dt        (500., 650.) == Approx { 0.181226279e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r3::massic_internal_energy_dt        (200., 650.) == Approx { 0.226365868e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r3::massic_internal_energy_dt        (500., 750.) == Approx { 0.210206932e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r3::massic_entropy_dt                (500., 650.) == Approx { 0.405427273e1 * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(r3::massic_entropy_dt                (200., 650.) == Approx { 0.485438792e1 * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(r3::massic_entropy_dt                (500., 750.) == Approx { 0.446971906e1 * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(r3::massic_isobaric_heat_capacity_dt (500., 650.) == Approx { 0.138935717e2 * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(r3::massic_isobaric_heat_capacity_dt (200., 650.) == Approx { 0.446579342e2 * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(r3::massic_isobaric_heat_capacity_dt (500., 750.) == Approx { 0.634165359e1 * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(r3::speed_of_sound_dt                (500., 650.) == Approx { 0.502005554e3       }.scale (1e3).epsilon (eps));
        CHECK(r3::speed_of_sound_dt                (200., 650.) == Approx { 0.383444594e3       }.scale (1e3).epsilon (eps));
        CHECK(r3::speed_of_sound_dt                (500., 750.) == Approx { 0.760696041e3       }.scale (1e3).epsilon (eps));

        CHECK(r3::b3ab_p (25e6) == Approx { 2.095936454e6 }.scale (1e6).epsilon (eps));
        CHECK(r3::temperature_ph ( 20e6, 1700e3) == Approx { 6.293083892e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ph ( 50e6, 2000e3) == Approx { 6.905718338e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ph (100e6, 2100e3) == Approx { 7.336163014e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ph ( 20e6, 2500e3) == Approx { 6.418418053e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ph ( 50e6, 2400e3) == Approx { 7.351848618e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ph (100e6, 2700e3) == Approx { 8.420460876e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::massic_volume_ph ( 20e6, 1700e3) == Approx { 1.749903962e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ph ( 50e6, 2000e3) == Approx { 1.908139035e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ph (100e6, 2100e3) == Approx { 1.676229776e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ph ( 20e6, 2500e3) == Approx { 6.670547043e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ph ( 50e6, 2400e3) == Approx { 2.801244590e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ph (100e6, 2700e3) == Approx { 2.404234998e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::temperature_ps ( 20e6, 3.8e3) == Approx { 6.282959869e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ps ( 50e6, 3.6e3) == Approx { 6.297158726e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ps (100e6, 4.0e3) == Approx { 7.056880237e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ps ( 20e6, 5.0e3) == Approx { 6.401176443e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ps ( 50e6, 4.5e3) == Approx { 7.163687517e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::temperature_ps (100e6, 5.0e3) == Approx { 8.474332825e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::massic_volume_ps ( 20e6,  3.8e3) == Approx { 1.733791463e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ps ( 50e6,  3.6e3) == Approx { 1.469680170e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ps (100e6,  4.0e3) == Approx { 1.555893131e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ps ( 20e6,  5.0e3) == Approx { 6.262101987e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ps ( 50e6,  4.5e3) == Approx { 2.332634294e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_ps (100e6,  5.0e3) == Approx { 2.449610757e-3 }.scale (1e-3).epsilon (eps));

        CHECK(r3::pressure_hs (1700e3, 3.8e3) == Approx { 2.555703246e7 }.scale (1e7).epsilon (eps));
        CHECK(r3::pressure_hs (2000e3, 4.2e3) == Approx { 4.540873468e7 }.scale (1e7).epsilon (eps));
        CHECK(r3::pressure_hs (2100e3, 4.3e3) == Approx { 6.078123340e7 }.scale (1e7).epsilon (eps));
        CHECK(r3::pressure_hs (2600e3, 5.1e3) == Approx { 3.434999263e7 }.scale (1e7).epsilon (eps));
        CHECK(r3::pressure_hs (2400e3, 4.7e3) == Approx { 6.363924887e7 }.scale (1e7).epsilon (eps));
        CHECK(r3::pressure_hs (2700e3, 5.0e3) == Approx { 8.839043281e7 }.scale (1e7).epsilon (eps));

        CHECK(r3::t_3ab_p (40.0e6) == Approx { 6.930341408e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3cd_p (25.0e6) == Approx { 6.493659208e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3ef_p (40.0e6) == Approx { 7.139593992e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3gh_p (23.0e6) == Approx { 6.498873759e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3ij_p (23.0e6) == Approx { 6.515778091e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3jk_p (23.0e6) == Approx { 6.558338344e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3mn_p (22.8e6) == Approx { 6.496054133e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3op_p (22.8e6) == Approx { 6.500106943e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3qu_p (22.0e6) == Approx { 6.456355027e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3rx_p (22.0e6) == Approx { 6.482622754e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3uv_p (22.3e6) == Approx { 6.477996121e2 }.scale (1e2).epsilon (eps));
        CHECK(r3::t_3wx_p (22.3e6) == Approx { 6.482049480e2 }.scale (1e2).epsilon (eps));

        CHECK(r3::detail::subregion_pt (50.0e6,   630.00) == r3::detail::a);
        CHECK(r3::detail::subregion_pt (80.0e6,   670.00) == r3::detail::a);
        CHECK(r3::detail::subregion_pt (50.0e6,   710.00) == r3::detail::b);
        CHECK(r3::detail::subregion_pt (80.0e6,   750.00) == r3::detail::b);
        CHECK(r3::detail::subregion_pt (20.0e6,   630.00) == r3::detail::c);
        CHECK(r3::detail::subregion_pt (30.0e6,   650.00) == r3::detail::c);
        CHECK(r3::detail::subregion_pt (26.0e6,   656.00) == r3::detail::d);
        CHECK(r3::detail::subregion_pt (30.0e6,   670.00) == r3::detail::d);
        CHECK(r3::detail::subregion_pt (26.0e6,   661.00) == r3::detail::e);
        CHECK(r3::detail::subregion_pt (30.0e6,   675.00) == r3::detail::e);
        CHECK(r3::detail::subregion_pt (26.0e6,   671.00) == r3::detail::f);
        CHECK(r3::detail::subregion_pt (30.0e6,   690.00) == r3::detail::f);
        CHECK(r3::detail::subregion_pt (23.6e6,   649.00) == r3::detail::g);
        CHECK(r3::detail::subregion_pt (24.0e6,   650.00) == r3::detail::g);
        CHECK(r3::detail::subregion_pt (23.6e6,   652.00) == r3::detail::h);
        CHECK(r3::detail::subregion_pt (24.0e6,   654.00) == r3::detail::h);
        CHECK(r3::detail::subregion_pt (23.6e6,   653.00) == r3::detail::i);
        CHECK(r3::detail::subregion_pt (24.0e6,   655.00) == r3::detail::i);
        CHECK(r3::detail::subregion_pt (23.5e6,   655.00) == r3::detail::j);
        CHECK(r3::detail::subregion_pt (24.0e6,   660.00) == r3::detail::j);
        CHECK(r3::detail::subregion_pt (23.0e6,   660.00) == r3::detail::k);
        CHECK(r3::detail::subregion_pt (24.0e6,   670.00) == r3::detail::k);
        CHECK(r3::detail::subregion_pt (22.6e6,   646.00) == r3::detail::l);
        CHECK(r3::detail::subregion_pt (23.0e6,   646.00) == r3::detail::l);
        CHECK(r3::detail::subregion_pt (22.6e6,   648.60) == r3::detail::m);
        CHECK(r3::detail::subregion_pt (22.8e6,   649.30) == r3::detail::m);
        CHECK(r3::detail::subregion_pt (22.6e6,   649.00) == r3::detail::n);
        CHECK(r3::detail::subregion_pt (22.8e6,   649.70) == r3::detail::n);
        CHECK(r3::detail::subregion_pt (22.6e6,   649.10) == r3::detail::o);
        CHECK(r3::detail::subregion_pt (22.8e6,   649.90) == r3::detail::o);
        CHECK(r3::detail::subregion_pt (22.6e6,   649.40) == r3::detail::p);
        CHECK(r3::detail::subregion_pt (22.8e6,   650.20) == r3::detail::p);
        CHECK(r3::detail::subregion_pt (21.1e6,   640.00) == r3::detail::q);
        CHECK(r3::detail::subregion_pt (21.8e6,   643.00) == r3::detail::q);
        CHECK(r3::detail::subregion_pt (21.1e6,   644.00) == r3::detail::r);
        CHECK(r3::detail::subregion_pt (21.8e6,   648.00) == r3::detail::r);
        CHECK(r3::detail::subregion_pt (19.1e6,   635.00) == r3::detail::s);
        CHECK(r3::detail::subregion_pt (20.0e6,   638.00) == r3::detail::s);
        CHECK(r3::detail::subregion_pt (17.0e6,   626.00) == r3::detail::t);
        CHECK(r3::detail::subregion_pt (20.0e6,   640.00) == r3::detail::t);
        CHECK(r3::detail::subregion_pt (21.5e6,   644.60) == r3::detail::u);
        CHECK(r3::detail::subregion_pt (22.0e6,   646.10) == r3::detail::u);
        CHECK(r3::detail::subregion_pt (22.5e6,   648.60) == r3::detail::v);
        CHECK(r3::detail::subregion_pt (22.3e6,   647.90) == r3::detail::v);
        CHECK(r3::detail::subregion_pt (22.15e6,  647.50) == r3::detail::w);
        CHECK(r3::detail::subregion_pt (22.3e6,   648.10) == r3::detail::w);
        CHECK(r3::detail::subregion_pt (22.11e6,  648.00) == r3::detail::x);
        CHECK(r3::detail::subregion_pt (22.3e6,   649.00) == r3::detail::x);
        CHECK(r3::detail::subregion_pt (22.0e6,   646.84) == r3::detail::y);
        CHECK(r3::detail::subregion_pt (22.064e6, 647.05) == r3::detail::y);
        CHECK(r3::detail::subregion_pt (22.0e6,   646.89) == r3::detail::z);
        CHECK(r3::detail::subregion_pt (22.064e6, 647.15) == r3::detail::z);

        CHECK(r3::massic_volume_pt (50.0e6,   630.00) == Approx { 1.470853100e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (80.0e6,   670.00) == Approx { 1.503831359e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (50.0e6,   710.00) == Approx { 2.204728587e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (80.0e6,   750.00) == Approx { 1.973692940e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (20.0e6,   630.00) == Approx { 1.761696406e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (30.0e6,   650.00) == Approx { 1.819560617e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (26.0e6,   656.00) == Approx { 2.245587720e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (30.0e6,   670.00) == Approx { 2.506897702e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (26.0e6,   661.00) == Approx { 2.970225962e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (30.0e6,   675.00) == Approx { 3.004627086e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (26.0e6,   671.00) == Approx { 5.019029401e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (30.0e6,   690.00) == Approx { 4.656470142e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (23.6e6,   649.00) == Approx { 2.163198378e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (24.0e6,   650.00) == Approx { 2.166044161e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (23.6e6,   652.00) == Approx { 2.651081407e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (24.0e6,   654.00) == Approx { 2.967802335e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (23.6e6,   653.00) == Approx { 3.273916816e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (24.0e6,   655.00) == Approx { 3.550329864e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (23.5e6,   655.00) == Approx { 4.545001142e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (24.0e6,   660.00) == Approx { 5.100267704e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (23.0e6,   660.00) == Approx { 6.109525997e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (24.0e6,   670.00) == Approx { 6.427325645e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.6e6,   646.00) == Approx { 2.117860851e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (23.0e6,   646.00) == Approx { 2.062374674e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.6e6,   648.60) == Approx { 2.533063780e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.8e6,   649.30) == Approx { 2.572971781e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.6e6,   649.00) == Approx { 2.923432711e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.8e6,   649.70) == Approx { 2.913311494e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.6e6,   649.10) == Approx { 3.131208996e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.8e6,   649.90) == Approx { 3.221160278e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.6e6,   649.40) == Approx { 3.715596186e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.8e6,   650.20) == Approx { 3.664754790e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (21.1e6,   640.00) == Approx { 1.970999272e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (21.8e6,   643.00) == Approx { 2.043919161e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (21.1e6,   644.00) == Approx { 5.251009921e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (21.8e6,   648.00) == Approx { 5.256844741e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (19.1e6,   635.00) == Approx { 1.932829079e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (20.0e6,   638.00) == Approx { 1.985387227e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (17.0e6,   626.00) == Approx { 8.483262001e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (20.0e6,   640.00) == Approx { 6.227528101e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (21.5e6,   644.60) == Approx { 2.268366647e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.0e6,   646.10) == Approx { 2.296350553e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.5e6,   648.60) == Approx { 2.832373260e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.3e6,   647.90) == Approx { 2.811424405e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.15e6,  647.50) == Approx { 3.694032281e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.3e6,   648.10) == Approx { 3.622226305e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.11e6,  648.00) == Approx { 4.528072649e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.3e6,   649.00) == Approx { 4.556905799e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.0e6,   646.84) == Approx { 2.698354719e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.064e6, 647.05) == Approx { 2.717655648e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.0e6,   646.89) == Approx { 3.798732962e-3 }.scale (1e-3).epsilon (eps));
        CHECK(r3::massic_volume_pt (22.064e6, 647.15) == Approx { 3.701940010e-3 }.scale (1e-3).epsilon (eps));
    };
    /*
    SUBCASE("iapws-r7-region-4")
    {
            using namespace r4;
        CHECK(saturation_pressure_t    (300. ) == Approx { 0.353658941e-2 * 1e6 }.scale (1e4).epsilon (eps));
        CHECK(saturation_pressure_t    (500. ) == Approx { 0.263889776e1  * 1e6 }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure_t    (600. ) == Approx { 0.123443146e2  * 1e6 }.scale (1e8).epsilon (eps));
        CHECK(saturation_temperature_p (0.1e6) == Approx { 0.372755919e3        }.scale (1e3).epsilon (eps));
        CHECK(saturation_temperature_p (1.0e6) == Approx { 0.453035632e3        }.scale (1e3).epsilon (eps));
        CHECK(saturation_temperature_p ( 10e6) == Approx { 0.584149488e3        }.scale (1e3).epsilon (eps));
    };
    */
    SUBCASE("iapws-r7-region-5")
    {
        CHECK(r5::massic_volume_pt                 (0.5e6, 1500.) == Approx { 0.138455090e1        }.scale (1e1).epsilon (eps));
        CHECK(r5::massic_volume_pt                 (30.e6, 1500.) == Approx { 0.230761299e-1       }.scale (1e-1).epsilon (eps));
        CHECK(r5::massic_volume_pt                 (30.e6, 2000.) == Approx { 0.311385219e-1       }.scale (1e-1).epsilon (eps));
        CHECK(r5::massic_enthalpy_pt               (0.5e6, 1500.) == Approx { 0.521976855e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(r5::massic_enthalpy_pt               (30.e6, 1500.) == Approx { 0.516723514e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(r5::massic_enthalpy_pt               (30.e6, 2000.) == Approx { 0.657122604e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(r5::massic_internal_energy_pt        (0.5e6, 1500.) == Approx { 0.452749310e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(r5::massic_internal_energy_pt        (30.e6, 1500.) == Approx { 0.447495124e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(r5::massic_internal_energy_pt        (30.e6, 2000.) == Approx { 0.563707038e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(r5::massic_entropy_pt                (0.5e6, 1500.) == Approx { 0.965408875e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r5::massic_entropy_pt                (30.e6, 1500.) == Approx { 0.772970133e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r5::massic_entropy_pt                (30.e6, 2000.) == Approx { 0.853640523e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r5::massic_isobaric_heat_capacity_pt (0.5e6, 1500.) == Approx { 0.261609445e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r5::massic_isobaric_heat_capacity_pt (30.e6, 1500.) == Approx { 0.272724317e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r5::massic_isobaric_heat_capacity_pt (30.e6, 2000.) == Approx { 0.288569882e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(r5::speed_of_sound_pt                (0.5e6, 1500.) == Approx { 0.917068690e3        }.scale (1e3).epsilon (eps));
        CHECK(r5::speed_of_sound_pt                (30.e6, 1500.) == Approx { 0.928548002e3        }.scale (1e3).epsilon (eps));
        CHECK(r5::speed_of_sound_pt                (30.e6, 2000.) == Approx { 0.106736948e4        }.scale (1e4).epsilon (eps));
    };
    SUBCASE("region P-T")
    {
        CHECK(p_b23_t (0.62315e3)     == Approx { 0.165291643e8 }.scale (1e8).epsilon (eps));
        CHECK(t_b23_p (0.165291643e8) == Approx { 0.62315e3     }.scale (1e3).epsilon (eps));
        CHECK(region_pt (50e6, 280.) == 1);
        CHECK(region_pt (50e6, 1070.) == 2);
        CHECK(region_pt (50e6, 630.) == 3);
        //CHECK(region_pt (0.353658941e-2 * 1e6, 300.) == 4);
        CHECK(region_pt (10e6, 1100.) == 5);
    }
    SUBCASE("region H-S")
    {
        CHECK(detail::h_p_1_s (1.0e3) == Approx { 3.085509647e5 }.scale (1e5).epsilon (eps));
        CHECK(detail::h_p_1_s (2.0e3) == Approx { 7.006304472e5 }.scale (1e5).epsilon (eps));
        CHECK(detail::h_p_1_s (3.0e3) == Approx { 1.198359754e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_p_3a_s (3.8e3) == Approx { 1.685025565e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_p_3a_s (4.0e3) == Approx { 1.816891476e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_p_3a_s (4.2e3) == Approx { 1.949352563e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_pp_2ab_s (7.0e3) == Approx { 2.723729985e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_pp_2ab_s (8.0e3) == Approx { 2.599047210e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_pp_2ab_s (9.0e3) == Approx { 2.511861477e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_pp_2c3b_s (5.5e3) == Approx { 2.687693850e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_pp_2c3b_s (5.0e3) == Approx { 2.451623609e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_pp_2c3b_s (4.5e3) == Approx { 2.144360448e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_b13_s (3.7e3) == Approx { 1.632525047e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_b13_s (3.6e3) == Approx { 1.593027214e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::h_b13_s (3.5e3) == Approx { 1.566104611e6 }.scale (1e6).epsilon (eps));
        CHECK(detail::t_b23_hs (2600e3, 5.10e3) == Approx { 7.135259364e2 }.scale (1e2).epsilon (eps));
        CHECK(detail::t_b23_hs (2700e3, 5.15e3) == Approx { 7.685345532e2 }.scale (1e2).epsilon (eps));
        CHECK(detail::t_b23_hs (2800e3, 5.20e3) == Approx { 8.176202120e2 }.scale (1e2).epsilon (eps));
        CHECK(region_hs (0.001e3,    0.) == 1);
        CHECK(region_hs (   90e3,    0.) == 1);
        CHECK(region_hs ( 1500e3, 3.4e3) == 1);
        CHECK(region_hs ( 2800e3, 6.5e3) == 2);
        CHECK(region_hs ( 2800e3, 9.5e3) == 2);
        CHECK(region_hs ( 4100e3, 9.5e3) == 2);
        CHECK(region_hs ( 2800e3, 6.0e3) == 2);
        CHECK(region_hs ( 3600e3, 6.0e3) == 2);
        CHECK(region_hs ( 3600e3, 7.0e3) == 2);
        CHECK(region_hs ( 2800e3, 5.1e3) == 2);
        CHECK(region_hs ( 2800e3, 5.8e3) == 2);
        CHECK(region_hs ( 3400e3, 5.8e3) == 2);
        CHECK(region_hs ( 1700e3, 3.8e3) == 3);
        CHECK(region_hs ( 2000e3, 4.2e3) == 3);
        CHECK(region_hs ( 2100e3, 4.3e3) == 3);
        CHECK(region_hs ( 2600e3, 5.1e3) == 3);
        CHECK(region_hs ( 2400e3, 4.7e3) == 3);
        CHECK(region_hs ( 2700e3, 5.0e3) == 3);
        CHECK(region_hs ( 1800e3, 5.3e3) == 4);
        CHECK(region_hs ( 2400e3, 6.0e3) == 4);
        CHECK(region_hs ( 2500e3, 5.5e3) == 4);
    }
    SUBCASE("region P-H")
    {
        CHECK(region_ph (    3e6,  500e3) == 1);
        CHECK(region_ph (   80e6,  500e3) == 1);
        CHECK(region_ph (   80e6, 1500e3) == 1);
        CHECK(region_ph (0.001e6, 3000e3) == 2);
        CHECK(region_ph (    3e6, 3000e3) == 2);
        CHECK(region_ph (    3e6, 4000e3) == 2);
        CHECK(region_ph (    5e6, 3500e3) == 2);
        CHECK(region_ph (    5e6, 4000e3) == 2);
        CHECK(region_ph (   25e6, 3500e3) == 2);
        CHECK(region_ph (   40e6, 2700e3) == 2);
        CHECK(region_ph (   60e6, 2700e3) == 2);
        CHECK(region_ph (   60e6, 3200e3) == 2);
        CHECK(region_ph (   20e6, 1700e3) == 3);
        CHECK(region_ph (   50e6, 2000e3) == 3);
        CHECK(region_ph (  100e6, 2100e3) == 3);
        CHECK(region_ph (   20e6, 2500e3) == 3);
        CHECK(region_ph (   50e6, 2400e3) == 3);
        CHECK(region_ph (  100e6, 2700e3) == 3);
        // 4?
    }
    SUBCASE("region P-S")
    {
        CHECK(region_ps (    3e6,  0.5e3) == 1);
        CHECK(region_ps (   80e6,  0.5e3) == 1);
        CHECK(region_ps (   80e6,  3.0e3) == 1);
        CHECK(region_ps (  0.1e6,  7.5e3) == 2);
        CHECK(region_ps (  0.1e6,    8e3) == 2);
        CHECK(region_ps (  2.5e6,    8e3) == 2);
        CHECK(region_ps (    8e6,    6e3) == 2);
        CHECK(region_ps (    8e6,  7.5e3) == 2);
        CHECK(region_ps (   90e6,    6e3) == 2);
        CHECK(region_ps (   20e6, 5.75e3) == 2);
        CHECK(region_ps (   80e6, 5.25e3) == 2);
        CHECK(region_ps (   80e6, 5.75e3) == 2);
        CHECK(region_ps (   20e6,  3.8e3) == 3);
        CHECK(region_ps (   50e6,  3.6e3) == 3);
        CHECK(region_ps (  100e6,  4.0e3) == 3);
        CHECK(region_ps (   20e6,  5.0e3) == 3);
        CHECK(region_ps (   50e6,  4.5e3) == 3);
        CHECK(region_ps (  100e6,  5.0e3) == 3);
        // 4?
    }
    SUBCASE("main API")
    {
        auto const t = 300.;
        auto const p = 3e6;
        auto const h = 1e6;
        auto const s = 1e3;
        // dp ?
        // pt
        massic_volume_pt                 (p, t);
        density_pt                       (p, t);
        massic_enthalpy_pt               (p, t);
        massic_internal_energy_pt        (p, t);
        massic_entropy_pt                (p, t);
        massic_isobaric_heat_capacity_pt (p, t);
        speed_of_sound_pt                (p, t);
        // hs
        pressure_hs                      (h, s);
        temperature_hs                   (h, s);
        massic_volume_hs                 (h, s);
        density_hs                       (h, s);
        massic_internal_energy_hs        (h, s);
        massic_isobaric_heat_capacity_hs (h, s);
        speed_of_sound_hs                (h, s);
        // ph
        temperature_ph                   (p, h);
        massic_volume_ph                 (p, h);
        density_ph                       (p, h);
        massic_internal_energy_ph        (p, h);
        massic_entropy_ph                (p, h);
        massic_isobaric_heat_capacity_ph (p, h);
        speed_of_sound_ph                (p, h);
        // ps
        temperature_ps                   (p, s);
        massic_volume_ps                 (p, s);
        density_ps                       (p, s);
        massic_enthalpy_ps               (p, s);
        massic_internal_energy_ps        (p, s);
        massic_isobaric_heat_capacity_ps (p, s);
        speed_of_sound_ps                (p, s);

        CHECK(massic_enthalpy_pt               (0.0035e6, 300.  ) == Approx { 0.254991145e4 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_enthalpy_pt               (0.0035e6, 700.  ) == Approx { 0.333568375e4 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_enthalpy_pt               (0.5e6   , 1500. ) == Approx { 0.521976855e4 * 1e3 } .scale (1e7 ).epsilon (eps));
        CHECK(massic_enthalpy_pt               (3. * 1e6, 300.  ) == Approx { 0.115331273e3 * 1e3 } .scale (1e6 ).epsilon (eps));
        CHECK(massic_enthalpy_pt               (3. * 1e6, 500.  ) == Approx { 0.975542239e3 * 1e3 } .scale (1e6 ).epsilon (eps));
        CHECK(massic_enthalpy_pt               (30.e6   , 1500. ) == Approx { 0.516723514e4 * 1e3 } .scale (1e7 ).epsilon (eps));
        CHECK(massic_enthalpy_pt               (30.e6   , 2000. ) == Approx { 0.657122604e4 * 1e3 } .scale (1e7 ).epsilon (eps));
        CHECK(massic_enthalpy_pt               (30e6    , 700.  ) == Approx { 0.263149474e4 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_enthalpy_pt               (80.* 1e6, 300.  ) == Approx { 0.184142828e3 * 1e3 } .scale (1e6 ).epsilon (eps));
        CHECK(massic_entropy_pt                (0.0035e6, 300.  ) == Approx { 0.852238967e1 * 1e3 } .scale (1e1 ).epsilon (eps));
        CHECK(massic_entropy_pt                (0.0035e6, 700.  ) == Approx { 0.101749996e2 * 1e3 } .scale (1e2 ).epsilon (eps));
        CHECK(massic_entropy_pt                (0.5e6   , 1500. ) == Approx { 0.965408875e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_entropy_pt                (3. * 1e6, 300.  ) == Approx { 0.392294792   * 1e3 } .scale (1e2 ).epsilon (eps));
        CHECK(massic_entropy_pt                (3. * 1e6, 500.  ) == Approx { 0.258041912e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_entropy_pt                (30.e6   , 1500. ) == Approx { 0.772970133e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_entropy_pt                (30.e6   , 2000. ) == Approx { 0.853640523e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_entropy_pt                (30e6    , 700.  ) == Approx { 0.517540298e1 * 1e3 } .scale (1e1 ).epsilon (eps));
        CHECK(massic_entropy_pt                (80.* 1e6, 300.  ) == Approx { 0.368563852   * 1e3 } .scale (1e3 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (0.0035e6, 300.  ) == Approx { 0.241169160e4 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (0.0035e6, 700.  ) == Approx { 0.301262819e4 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (0.5e6   , 1500. ) == Approx { 0.452749310e4 * 1e3 } .scale (1e7 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (3. * 1e6, 300.  ) == Approx { 0.112324818e3 * 1e3 } .scale (1e6 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (3. * 1e6, 500.  ) == Approx { 0.971934985e3 * 1e3 } .scale (1e6 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (30.e6   , 1500. ) == Approx { 0.447495124e4 * 1e3 } .scale (1e7 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (30.e6   , 2000. ) == Approx { 0.563707038e4 * 1e3 } .scale (1e7 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (30e6    , 700.  ) == Approx { 0.246861076e4 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_internal_energy_pt        (80.* 1e6, 300.  ) == Approx { 0.106448356e3 * 1e3 } .scale (1e6 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (0.0035e6, 300.  ) == Approx { 0.191300162e1 * 1e3 } .scale (1e1 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (0.0035e6, 700.  ) == Approx { 0.208141274e1 * 1e3 } .scale (1e1 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (0.5e6   , 1500. ) == Approx { 0.261609445e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (3. * 1e6, 300.  ) == Approx { 0.417301218e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (3. * 1e6, 500.  ) == Approx { 0.465580682e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (30.e6   , 1500. ) == Approx { 0.272724317e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (30.e6   , 2000. ) == Approx { 0.288569882e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (30e6    , 700.  ) == Approx { 0.103505092e2 * 1e3 } .scale (1e2 ).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity_pt (80.* 1e6, 300.  ) == Approx { 0.401008987e1 * 1e3 } .scale (1e4 ).epsilon (eps));
        CHECK(massic_volume_ph                 (100e6   , 2100e3) == Approx { 1.676229776e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ph                 (100e6   , 2700e3) == Approx { 2.404234998e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ph                 (20e6    , 1700e3) == Approx { 1.749903962e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ph                 (20e6    , 2500e3) == Approx { 6.670547043e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ph                 (50e6    , 2000e3) == Approx { 1.908139035e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ph                 (50e6    , 2400e3) == Approx { 2.801244590e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ps                 (100e6   , 4.0e3 ) == Approx { 1.555893131e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ps                 (100e6   , 5.0e3 ) == Approx { 2.449610757e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ps                 (20e6    , 3.8e3 ) == Approx { 1.733791463e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ps                 (20e6    , 5.0e3 ) == Approx { 6.262101987e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ps                 (50e6    , 3.6e3 ) == Approx { 1.469680170e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_ps                 (50e6    , 4.5e3 ) == Approx { 2.332634294e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (0.0035e6, 300.  ) == Approx { 0.394913866e2       } .scale (1e2 ).epsilon (eps));
        CHECK(massic_volume_pt                 (0.0035e6, 700.  ) == Approx { 0.923015898e2       } .scale (1e2 ).epsilon (eps));
        CHECK(massic_volume_pt                 (0.5e6   , 1500. ) == Approx { 0.138455090e1       } .scale (1e1 ).epsilon (eps));
        CHECK(massic_volume_pt                 (17.0e6  , 626.00) == Approx { 8.483262001e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (19.1e6  , 635.00) == Approx { 1.932829079e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (20.0e6  , 630.00) == Approx { 1.761696406e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (20.0e6  , 638.00) == Approx { 1.985387227e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (20.0e6  , 640.00) == Approx { 6.227528101e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (21.1e6  , 640.00) == Approx { 1.970999272e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (21.1e6  , 644.00) == Approx { 5.251009921e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (21.5e6  , 644.60) == Approx { 2.268366647e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (21.8e6  , 643.00) == Approx { 2.043919161e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (21.8e6  , 648.00) == Approx { 5.256844741e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.064e6, 647.05) == Approx { 2.717655648e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.064e6, 647.15) == Approx { 3.701940010e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.0e6  , 646.10) == Approx { 2.296350553e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.0e6  , 646.84) == Approx { 2.698354719e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.0e6  , 646.89) == Approx { 3.798732962e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.11e6 , 648.00) == Approx { 4.528072649e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.15e6 , 647.50) == Approx { 3.694032281e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.3e6  , 647.90) == Approx { 2.811424405e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.3e6  , 648.10) == Approx { 3.622226305e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.3e6  , 649.00) == Approx { 4.556905799e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.5e6  , 648.60) == Approx { 2.832373260e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.6e6  , 646.00) == Approx { 2.117860851e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.6e6  , 648.60) == Approx { 2.533063780e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.6e6  , 649.00) == Approx { 2.923432711e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.6e6  , 649.10) == Approx { 3.131208996e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.6e6  , 649.40) == Approx { 3.715596186e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.8e6  , 649.30) == Approx { 2.572971781e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.8e6  , 649.70) == Approx { 2.913311494e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.8e6  , 649.90) == Approx { 3.221160278e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (22.8e6  , 650.20) == Approx { 3.664754790e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (23.0e6  , 646.00) == Approx { 2.062374674e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (23.0e6  , 660.00) == Approx { 6.109525997e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (23.5e6  , 655.00) == Approx { 4.545001142e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (23.6e6  , 649.00) == Approx { 2.163198378e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (23.6e6  , 652.00) == Approx { 2.651081407e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (23.6e6  , 653.00) == Approx { 3.273916816e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (24.0e6  , 650.00) == Approx { 2.166044161e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (24.0e6  , 654.00) == Approx { 2.967802335e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (24.0e6  , 655.00) == Approx { 3.550329864e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (24.0e6  , 660.00) == Approx { 5.100267704e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (24.0e6  , 670.00) == Approx { 6.427325645e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (26.0e6  , 656.00) == Approx { 2.245587720e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (26.0e6  , 661.00) == Approx { 2.970225962e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (26.0e6  , 671.00) == Approx { 5.019029401e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (3. * 1e6, 300.  ) == Approx { 0.100215168e-2      } .scale (1e-2).epsilon (eps));
        CHECK(massic_volume_pt                 (3. * 1e6, 500.  ) == Approx { 0.120241800e-2      } .scale (1e-2).epsilon (eps));
        CHECK(massic_volume_pt                 (30.0e6  , 650.00) == Approx { 1.819560617e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (30.0e6  , 670.00) == Approx { 2.506897702e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (30.0e6  , 675.00) == Approx { 3.004627086e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (30.0e6  , 690.00) == Approx { 4.656470142e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (30.e6   , 1500. ) == Approx { 0.230761299e-1      } .scale (1e-1).epsilon (eps));
        CHECK(massic_volume_pt                 (30.e6   , 2000. ) == Approx { 0.311385219e-1      } .scale (1e-1).epsilon (eps));
        CHECK(massic_volume_pt                 (30e6    , 700.  ) == Approx { 0.542946619e-2      } .scale (1e-2).epsilon (eps));
        CHECK(massic_volume_pt                 (50.0e6  , 630.00) == Approx { 1.470853100e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (50.0e6  , 710.00) == Approx { 2.204728587e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (80.* 1e6, 300.  ) == Approx { 0.971180894e-3      } .scale (1e-2).epsilon (eps));
        CHECK(massic_volume_pt                 (80.0e6  , 670.00) == Approx { 1.503831359e-3      } .scale (1e-3).epsilon (eps));
        CHECK(massic_volume_pt                 (80.0e6  , 750.00) == Approx { 1.973692940e-3      } .scale (1e-3).epsilon (eps));
        CHECK(pressure_hs                      (0.001e3 , 0.    ) == Approx { 9.800980612e-4      } .scale (1e4 ).epsilon (eps));
        CHECK(pressure_hs                      (1500e3  , 3.4e3 ) == Approx { 5.868294423e1       } .scale (1e1 ).epsilon (eps));
        CHECK(pressure_hs                      (1700e3  , 3.8e3 ) == Approx { 2.555703246e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (2000e3  , 4.2e3 ) == Approx { 4.540873468e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (2100e3  , 4.3e3 ) == Approx { 6.078123340e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (2400e3  , 4.7e3 ) == Approx { 6.363924887e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (2600e3  , 5.1e3 ) == Approx { 3.434999263e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (2700e3  , 5.0e3 ) == Approx { 8.839043281e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (2800e3  , 5.1e3 ) == Approx { 9.439202060e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (2800e3  , 5.8e3 ) == Approx { 8.414574124e6       } .scale (1e6 ).epsilon (eps));
        CHECK(pressure_hs                      (2800e3  , 6.0e3 ) == Approx { 4.793911442e6       } .scale (1e6 ).epsilon (eps));
        CHECK(pressure_hs                      (2800e3  , 6.5e3 ) == Approx { 1.371012767e6       } .scale (1e6 ).epsilon (eps));
        CHECK(pressure_hs                      (2800e3  , 9.5e3 ) == Approx { 1.879743844e3       } .scale (1e3 ).epsilon (eps));
        CHECK(pressure_hs                      (3400e3  , 5.8e3 ) == Approx { 8.376903879e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (3600e3  , 6.0e3 ) == Approx { 8.395519209e7       } .scale (1e7 ).epsilon (eps));
        CHECK(pressure_hs                      (3600e3  , 7.0e3 ) == Approx { 7.527161441e6       } .scale (1e6 ).epsilon (eps));
        CHECK(pressure_hs                      (4100e3  , 9.5e3 ) == Approx { 1.024788997e5       } .scale (1e5 ).epsilon (eps));
        CHECK(pressure_hs                      (90e3    , 0.    ) == Approx { 9.192954727e1       } .scale (1e1 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (0.0035e6, 300.  ) == Approx { 0.427920172e3       } .scale (1e3 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (0.0035e6, 700.  ) == Approx { 0.644289068e3       } .scale (1e3 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (0.5e6   , 1500. ) == Approx { 0.917068690e3       } .scale (1e3 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (3. * 1e6, 300.  ) == Approx { 0.150773921e4       } .scale (1e4 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (3. * 1e6, 500.  ) == Approx { 0.124071337e4       } .scale (1e4 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (30.e6   , 1500. ) == Approx { 0.928548002e3       } .scale (1e3 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (30.e6   , 2000. ) == Approx { 0.106736948e4       } .scale (1e4 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (30e6    , 700.  ) == Approx { 0.480386523e3       } .scale (1e3 ).epsilon (eps));
        CHECK(speed_of_sound_pt                (80.* 1e6, 300.  ) == Approx { 0.163469054e4       } .scale (1e4 ).epsilon (eps));
        CHECK(temperature_ph                   (0.001e6 , 3000e3) == Approx { 0.534433241e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (100e6   , 2100e3) == Approx { 7.336163014e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ph                   (100e6   , 2700e3) == Approx { 8.420460876e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ph                   (20e6    , 1700e3) == Approx { 6.293083892e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ph                   (20e6    , 2500e3) == Approx { 6.418418053e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ph                   (25e6    , 3500e3) == Approx { 0.875279054e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (3e6     , 3000e3) == Approx { 0.575373370e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (3e6     , 4000e3) == Approx { 0.101077577e4       } .scale (1e4 ).epsilon (eps));
        CHECK(temperature_ph                   (3e6     , 500e3 ) == Approx { 0.391798509e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (40e6    , 2700e3) == Approx { 0.743056411e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (50e6    , 2000e3) == Approx { 6.905718338e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ph                   (50e6    , 2400e3) == Approx { 7.351848618e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ph                   (5e6     , 3500e3) == Approx { 0.801299102e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (5e6     , 4000e3) == Approx { 0.101531583e4       } .scale (1e4 ).epsilon (eps));
        CHECK(temperature_ph                   (60e6    , 2700e3) == Approx { 0.791137067e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (60e6    , 3200e3) == Approx { 0.882756860e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (80e6    , 1500e3) == Approx { 0.611041229e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ph                   (80e6    , 500e3 ) == Approx { 0.378108626e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (0.1e6   , 7.5e3 ) == Approx { 0.399517097e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (0.1e6   , 8e3   ) == Approx { 0.514127081e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (100e6   , 4.0e3 ) == Approx { 7.056880237e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ps                   (100e6   , 5.0e3 ) == Approx { 8.474332825e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ps                   (2.5e6   , 8e3   ) == Approx { 0.103984917e4       } .scale (1e4 ).epsilon (eps));
        CHECK(temperature_ps                   (20e6    , 3.8e3 ) == Approx { 6.282959869e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ps                   (20e6    , 5.0e3 ) == Approx { 6.401176443e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ps                   (20e6    , 5.75e3) == Approx { 0.697992849e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (3e6     , 0.5e3 ) == Approx { 0.307842258e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (50e6    , 3.6e3 ) == Approx { 6.297158726e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ps                   (50e6    , 4.5e3 ) == Approx { 7.163687517e2       } .scale (1e2 ).epsilon (eps));
        CHECK(temperature_ps                   (80e6    , 0.5e3 ) == Approx { 0.309979785e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (80e6    , 3.0e3 ) == Approx { 0.565899909e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (80e6    , 5.25e3) == Approx { 0.854011484e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (80e6    , 5.75e3) == Approx { 0.949017998e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (8e6     , 6e3   ) == Approx { 0.600484040e3       } .scale (1e3 ).epsilon (eps));
        CHECK(temperature_ps                   (8e6     , 7.5e3 ) == Approx { 0.106495556e4       } .scale (1e4 ).epsilon (eps));
        CHECK(temperature_ps                   (90e6    , 6e3   ) == Approx { 0.103801126e4       } .scale (1e4 ).epsilon (eps));
    }
    SUBCASE("expansion, compressibility, etc.")
    {
            using namespace detail;
        for (auto it = 0u; it != T.size (); ++it)
        {
                auto const
            t = T.at (it) + 273.15;
            for (auto ip = 0u; ip != P.size (); ++ip)
            {
                    auto const
                p = P.at (ip) * 1e5;
                    auto const
                a = table_9.at (it).at (ip) * 1e-6;
                    auto const
                b = table_10.at (it).at (ip) * 1e-9;
                    auto const
                c = table_19.at (it).at (ip) * 1e-3;
                    auto const
                d = table_20.at (it).at (ip);
                    auto const
                r = region_pt (p, t);
                INFO ("P = ", p * 1e-5, ", T = ", t - 273.15, ", region = ", r);
                CHECK(isobaric_cubic_expansion_coefficient_pt (p, t) == Approx { a }.scale (fabs (a)).epsilon (1e-4));
                CHECK(isothermal_compressibility_pt           (p, t) == Approx { b }.scale (b).epsilon (1e-4));
                CHECK(relative_pressure_coefficient_pt        (p, t) == Approx { c }.scale (fabs (c)).epsilon (1e-4));
                CHECK(isothermal_stress_coefficient_pt        (p, t) == Approx { d }.scale (fabs (d)).epsilon (1e-4));
            }
        }
    }
}

