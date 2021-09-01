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
        CHECK(saturation_pressure (300.) == Approx { 0.353658941e-2 * 1e6 }.scale (1e4).epsilon (eps));
        CHECK(saturation_pressure (500.) == Approx { 0.263889776e1 * 1e6  }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure (600.) == Approx { 0.123443146e2 * 1e6  }.scale (1e8).epsilon (eps));
        CHECK(saturation_temperature (0.1 * 1e6) == Approx { 0.372755919e3 }.scale (1e3).epsilon (eps));
        CHECK(saturation_temperature (1.  * 1e6) == Approx { 0.453035632e3 }.scale (1e3).epsilon (eps));
        CHECK(saturation_temperature (10. * 1e6) == Approx { 0.584149488e3 }.scale (1e3).epsilon (eps));
    }
    SUBCASE("regions")
    {
        CHECK(detail::b23  (0.62315e3)     == Approx { 0.165291643e8 }.scale (1e8).epsilon (eps));
        CHECK(detail::b23i (0.165291643e8) == Approx { 0.62315e3     }.scale (1e3).epsilon (1e-9));
        CHECK(detail::region (50e6, 280.) == 1);
        CHECK(detail::region (50e6, 1070.) == 2);
        CHECK(detail::region (50e6, 630.) == 3);
        //CHECK(detail::region (0.353658941e-2 * 1e6, 300.) == 4);
        CHECK(detail::region (10e6, 1100.) == 5);
    }
    SUBCASE("iapws-r7-region-1")
    {
            using namespace detail::r1;
        CHECK(massic_volume                 (3. * 1e6, 300.) == Approx { 0.100215168e-2 }.scale (1e-2).epsilon (eps));
        CHECK(massic_volume                 (80.* 1e6, 300.) == Approx { 0.971180894e-3 }.scale (1e-2).epsilon (eps));
        CHECK(massic_volume                 (3. * 1e6, 500.) == Approx { 0.120241800e-2 }.scale (1e-2).epsilon (eps));
        CHECK(massic_enthalpy               (3. * 1e6, 300.) == Approx { 0.115331273e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(massic_enthalpy               (80.* 1e6, 300.) == Approx { 0.184142828e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(massic_enthalpy               (3. * 1e6, 500.) == Approx { 0.975542239e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(massic_internal_energy        (3. * 1e6, 300.) == Approx { 0.112324818e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(massic_internal_energy        (80.* 1e6, 300.) == Approx { 0.106448356e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(massic_internal_energy        (3. * 1e6, 500.) == Approx { 0.971934985e3 * 1e3 }.scale (1e6).epsilon (eps));
        CHECK(massic_entropy                (3. * 1e6, 300.) == Approx { 0.392294792   * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(massic_entropy                (80.* 1e6, 300.) == Approx { 0.368563852   * 1e3 }.scale (1e3).epsilon (eps));
        CHECK(massic_entropy                (3. * 1e6, 500.) == Approx { 0.258041912e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (3. * 1e6, 300.) == Approx { 0.417301218e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (80.* 1e6, 300.) == Approx { 0.401008987e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (3. * 1e6, 500.) == Approx { 0.465580682e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(speed_of_sound                (3. * 1e6, 300.) == Approx { 0.150773921e4  }.scale (1e4).epsilon (eps));
        CHECK(speed_of_sound                (80.* 1e6, 300.) == Approx { 0.163469054e4  }.scale (1e4).epsilon (eps));
        CHECK(speed_of_sound                (3. * 1e6, 500.) == Approx { 0.124071337e4  }.scale (1e4).epsilon (eps));

        CHECK(pressure (0.001e3,    0.) == Approx { 9.800980612e-4  }.scale (1e4).epsilon (eps));
        CHECK(pressure (   90e3,    0.) == Approx { 9.192954727e1  }.scale (1e1).epsilon (eps));
        CHECK(pressure ( 1500e3, 3.4e3) == Approx { 5.868294423e1  }.scale (1e1).epsilon (eps));
    };
    SUBCASE("iapws-r7-region-1-backward")
    {
            using namespace detail;
        CHECK(r1b1::temperature ( 3e6,  500e3) == Approx { 0.391798509e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1b1::temperature (80e6,  500e3) == Approx { 0.378108626e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1b1::temperature (80e6, 1500e3) == Approx { 0.611041229e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1b2::temperature ( 3e6,  0.5e3) == Approx { 0.307842258e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1b2::temperature (80e6,  0.5e3) == Approx { 0.309979785e3 }.scale (1e3).epsilon (eps)); 
        CHECK(r1b2::temperature (80e6,  3.0e3) == Approx { 0.565899909e3 }.scale (1e3).epsilon (eps));
    };
    SUBCASE("iapws-r7-region-2")
    {
            using namespace detail::r2;
        CHECK(massic_volume                 (0.0035e6, 300.) == Approx { 0.394913866e2        }.scale (1e2).epsilon (eps));
        CHECK(massic_volume                 (0.0035e6, 700.) == Approx { 0.923015898e2        }.scale (1e2).epsilon (eps));
        CHECK(massic_volume                 (    30e6, 700.) == Approx { 0.542946619e-2       }.scale (1e-2).epsilon (eps));
        CHECK(massic_enthalpy               (0.0035e6, 300.) == Approx { 0.254991145e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_enthalpy               (0.0035e6, 700.) == Approx { 0.333568375e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_enthalpy               (    30e6, 700.) == Approx { 0.263149474e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (0.0035e6, 300.) == Approx { 0.241169160e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (0.0035e6, 700.) == Approx { 0.301262819e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (    30e6, 700.) == Approx { 0.246861076e4  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_entropy                (0.0035e6, 300.) == Approx { 0.852238967e1  * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(massic_entropy                (0.0035e6, 700.) == Approx { 0.101749996e2  * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(massic_entropy                (    30e6, 700.) == Approx { 0.517540298e1  * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (0.0035e6, 300.) == Approx { 0.191300162e1  * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (0.0035e6, 700.) == Approx { 0.208141274e1  * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (    30e6, 700.) == Approx { 0.103505092e2  * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(speed_of_sound                (0.0035e6, 300.) == Approx { 0.427920172e3        }.scale (1e3).epsilon (eps));
        CHECK(speed_of_sound                (0.0035e6, 700.) == Approx { 0.644289068e3        }.scale (1e3).epsilon (eps));
        CHECK(speed_of_sound                (    30e6, 700.) == Approx { 0.480386523e3        }.scale (1e3).epsilon (eps));
    };
    SUBCASE("iapws-r7-region-2-metastable-vapor")
    {
            using namespace detail::r2_metastable_vapor;
        CHECK(massic_volume                 (1.0e6, 450.) == Approx { 0.192516540         }.scale (1e4).epsilon (eps));
        CHECK(massic_volume                 (1.0e6, 440.) == Approx { 0.186212297         }.scale (1e4).epsilon (eps));
        CHECK(massic_volume                 (1.5e6, 450.) == Approx { 0.121685206         }.scale (1e4).epsilon (eps));
        CHECK(massic_enthalpy               (1.0e6, 450.) == Approx { 0.276881115e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_enthalpy               (1.0e6, 440.) == Approx { 0.274015123e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_enthalpy               (1.5e6, 450.) == Approx { 0.272134539e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (1.0e6, 450.) == Approx { 0.257629461e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (1.0e6, 440.) == Approx { 0.255393894e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (1.5e6, 450.) == Approx { 0.253881758e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_entropy                (1.0e6, 450.) == Approx { 0.656660377e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_entropy                (1.0e6, 440.) == Approx { 0.650218759e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_entropy                (1.5e6, 450.) == Approx { 0.629170440e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (1.0e6, 450.) == Approx { 0.276349265e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (1.0e6, 440.) == Approx { 0.298166443e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (1.5e6, 450.) == Approx { 0.362795578e1 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(speed_of_sound                (1.0e6, 450.) == Approx { 0.498408101e3       }.scale (1e4).epsilon (eps));
        CHECK(speed_of_sound                (1.0e6, 440.) == Approx { 0.489363295e3       }.scale (1e4).epsilon (eps));
        CHECK(speed_of_sound                (1.5e6, 450.) == Approx { 0.481941819e3       }.scale (1e4).epsilon (eps));

    };
    SUBCASE("iapws-r7-region-2-backward")
    {
            using namespace detail::r2b;
        CHECK(B2bc_p (0.100000000e9)  == Approx { 0.3516004323e7 }.scale (1e7).epsilon (eps));
        CHECK(B2bc_h (0.3516004323e7) == Approx { 0.100000000e9  }.scale (1e9).epsilon (eps));

        CHECK(region_h (0.001e6, 3000e3) == 1);
        CHECK(region_h (    3e6, 3000e3) == 1);
        CHECK(region_h (    3e6, 4000e3) == 1);
        CHECK(region_h (    5e6, 3500e3) == 2);
        CHECK(region_h (    5e6, 4000e3) == 2);
        CHECK(region_h (   25e6, 3500e3) == 2);
        CHECK(region_h (   40e6, 2700e3) == 3);
        CHECK(region_h (   60e6, 2700e3) == 3);
        CHECK(region_h (   60e6, 3200e3) == 3);

        CHECK(region_s (  0.1e6,  7.5e3) == 1);
        CHECK(region_s (  0.1e6,    8e3) == 1);
        CHECK(region_s (  2.5e6,    8e3) == 1);
        CHECK(region_s (    8e6,    6e3) == 2);
        CHECK(region_s (    8e6,  7.5e3) == 2);
        CHECK(region_s (   90e6,    6e3) == 2);
        CHECK(region_s (   20e6, 5.75e3) == 3);
        CHECK(region_s (   80e6, 5.25e3) == 3);
        CHECK(region_s (   80e6, 5.75e3) == 3);

        CHECK(temperature_h (0.001e6, 3000e3) == Approx { 0.534433241e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_h (    3e6, 3000e3) == Approx { 0.575373370e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_h (    3e6, 4000e3) == Approx { 0.101077577e4 }.scale (1e4).epsilon (eps));
        CHECK(temperature_h (    5e6, 3500e3) == Approx { 0.801299102e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_h (    5e6, 4000e3) == Approx { 0.101531583e4 }.scale (1e4).epsilon (eps));
        CHECK(temperature_h (   25e6, 3500e3) == Approx { 0.875279054e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_h (   40e6, 2700e3) == Approx { 0.743056411e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_h (   60e6, 2700e3) == Approx { 0.791137067e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_h (   60e6, 3200e3) == Approx { 0.882756860e3 }.scale (1e3).epsilon (eps));

        CHECK(temperature_s (  0.1e6,  7.5e3) == Approx { 0.399517097e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_s (  0.1e6,    8e3) == Approx { 0.514127081e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_s (  2.5e6,    8e3) == Approx { 0.103984917e4 }.scale (1e4).epsilon (eps));
        CHECK(temperature_s (    8e6,    6e3) == Approx { 0.600484040e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_s (    8e6,  7.5e3) == Approx { 0.106495556e4 }.scale (1e4).epsilon (eps));
        CHECK(temperature_s (   90e6,    6e3) == Approx { 0.103801126e4 }.scale (1e4).epsilon (eps));
        CHECK(temperature_s (   20e6, 5.75e3) == Approx { 0.697992849e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_s (   80e6, 5.25e3) == Approx { 0.854011484e3 }.scale (1e3).epsilon (eps));
        CHECK(temperature_s (   80e6, 5.75e3) == Approx { 0.949017998e3 }.scale (1e3).epsilon (eps));
    };
    SUBCASE("iapws-r7-region-3")
    {
            using namespace detail::r3;
        CHECK(pressure                      (500., 650.) == Approx { 0.255837018e2 * 1e6 }.scale (1e2).epsilon (eps));
        CHECK(pressure                      (200., 650.) == Approx { 0.222930643e2 * 1e6 }.scale (1e2).epsilon (eps));
        CHECK(pressure                      (500., 750.) == Approx { 0.783095639e2 * 1e6 }.scale (1e2).epsilon (eps));
        CHECK(massic_enthalpy               (500., 650.) == Approx { 0.186343019e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_enthalpy               (200., 650.) == Approx { 0.237512401e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_enthalpy               (500., 750.) == Approx { 0.225868845e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (500., 650.) == Approx { 0.181226279e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (200., 650.) == Approx { 0.226365868e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_internal_energy        (500., 750.) == Approx { 0.210206932e4 * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_entropy                (500., 650.) == Approx { 0.405427273e1 * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(massic_entropy                (200., 650.) == Approx { 0.485438792e1 * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(massic_entropy                (500., 750.) == Approx { 0.446971906e1 * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (500., 650.) == Approx { 0.138935717e2 * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (200., 650.) == Approx { 0.446579342e2 * 1e3 }.scale (1e2).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (500., 750.) == Approx { 0.634165359e1 * 1e3 }.scale (1e1).epsilon (eps));
        CHECK(speed_of_sound                (500., 650.) == Approx { 0.502005554e3       }.scale (1e3).epsilon (eps));
        CHECK(speed_of_sound                (200., 650.) == Approx { 0.383444594e3       }.scale (1e3).epsilon (eps));
        CHECK(speed_of_sound                (500., 750.) == Approx { 0.760696041e3       }.scale (1e3).epsilon (eps));
    };
    SUBCASE("iapws-r7-region-4")
    {
            using namespace detail::r4;
        CHECK(saturation_pressure    (300. ) == Approx { 0.353658941e-2 * 1e6 }.scale (1e4).epsilon (eps));
        CHECK(saturation_pressure    (500. ) == Approx { 0.263889776e1  * 1e6 }.scale (1e7).epsilon (eps));
        CHECK(saturation_pressure    (600. ) == Approx { 0.123443146e2  * 1e6 }.scale (1e8).epsilon (eps));
        CHECK(saturation_temperature (0.1e6) == Approx { 0.372755919e3        }.scale (1e3).epsilon (eps));
        CHECK(saturation_temperature (1.0e6) == Approx { 0.453035632e3        }.scale (1e3).epsilon (eps));
        CHECK(saturation_temperature ( 10e6) == Approx { 0.584149488e3        }.scale (1e3).epsilon (eps));
    };
    SUBCASE("iapws-r7-region-5")
    {
            using namespace detail::r5;
        CHECK(massic_volume                 (0.5e6, 1500.) == Approx { 0.138455090e1        }.scale (1e1).epsilon (eps));
        CHECK(massic_volume                 (30.e6, 1500.) == Approx { 0.230761299e-1       }.scale (1e-1).epsilon (eps));
        CHECK(massic_volume                 (30.e6, 2000.) == Approx { 0.311385219e-1       }.scale (1e-1).epsilon (eps));
        CHECK(massic_enthalpy               (0.5e6, 1500.) == Approx { 0.521976855e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(massic_enthalpy               (30.e6, 1500.) == Approx { 0.516723514e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(massic_enthalpy               (30.e6, 2000.) == Approx { 0.657122604e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(massic_internal_energy        (0.5e6, 1500.) == Approx { 0.452749310e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(massic_internal_energy        (30.e6, 1500.) == Approx { 0.447495124e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(massic_internal_energy        (30.e6, 2000.) == Approx { 0.563707038e4  * 1e3 }.scale (1e7).epsilon (eps));
        CHECK(massic_entropy                (0.5e6, 1500.) == Approx { 0.965408875e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_entropy                (30.e6, 1500.) == Approx { 0.772970133e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_entropy                (30.e6, 2000.) == Approx { 0.853640523e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (0.5e6, 1500.) == Approx { 0.261609445e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (30.e6, 1500.) == Approx { 0.272724317e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(massic_isobaric_heat_capacity (30.e6, 2000.) == Approx { 0.288569882e1  * 1e3 }.scale (1e4).epsilon (eps));
        CHECK(speed_of_sound                (0.5e6, 1500.) == Approx { 0.917068690e3        }.scale (1e3).epsilon (eps));
        CHECK(speed_of_sound                (30.e6, 1500.) == Approx { 0.928548002e3        }.scale (1e3).epsilon (eps));
        CHECK(speed_of_sound                (30.e6, 2000.) == Approx { 0.106736948e4        }.scale (1e4).epsilon (eps));

    };
}

