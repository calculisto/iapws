#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#define ISTO_IAPWS_FORCE_RELAXED 1
#include "../include/isto/iapws/r7.hpp"
    using namespace isto::iapws::r7;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"

TEST_CASE("r7.hpp (relaxed)")
{
    SUBCASE("saturation line")
    {
        CHECK(saturation_pressure (300.) == Approx { 0.353658941e-2 * 1e6 }.epsilon (1e-8));
        CHECK(saturation_pressure (500.) == Approx { 0.263889776e1 * 1e6 }.epsilon (1e-8));
        CHECK(saturation_pressure (600.) == Approx { 0.123443146e2 * 1e6 }.epsilon (1e-8));
        CHECK(saturation_temperature (0.1 * 1e6) == Approx { 0.372755919e3 }.epsilon (1e-8));
        CHECK(saturation_temperature (1.  * 1e6) == Approx { 0.453035632e3 }.epsilon (1e-8));
        CHECK(saturation_temperature (10. * 1e6) == Approx { 0.584149488e3 }.epsilon (1e-8));
    }
    SUBCASE("regions")
    {
        CHECK(detail::b23 (0.62315e3) == Approx { 0.165291643e8 }.epsilon (1e-8));
        CHECK(detail::b23i (0.165291643e8) == Approx { 0.62315e3 }.epsilon (1e-9));
        CHECK(region (50e6, 280.) == 1);
        CHECK(region (50e6, 1070.) == 2);
        CHECK(region (50e6, 630.) == 3);
        //CHECK(region (0.353658941e-2 * 1e6, 300.) == 4);
        CHECK(region (10e6, 1100.) == 5);
    }
#if 0
    SUBCASE("iapws-r7-region-1")
    {
        eq_up_t <6> (massic_volume                 (300., 3. ), 0.100215168e-2, 1e-6);
        eq_up_t <6> (massic_volume                 (300., 80.), 0.971180894e-3, 1e-6);
        eq_up_t <6> (massic_volume                 (500., 3. ), 0.120241800e-2, 1e-6);
        eq_up_t <6> (massic_enthalpy               (300., 3. ), 0.115331273e3, 1e-5);
        eq_up_t <6> (massic_enthalpy               (300., 80.), 0.184142828e3, 1e-5);
        eq_up_t <6> (massic_enthalpy               (500., 3. ), 0.975542239e3, 1e-6);
        eq_up_t <6> (massic_internal_energy        (300., 3. ), 0.112324818e3, 1e-5);
        eq_up_t <6> (massic_internal_energy        (300., 80.), 0.106448356e3, 1e-5);
        eq_up_t <6> (massic_internal_energy        (500., 3. ), 0.971934985e3, 1e-6);
        eq_up_t <6> (massic_entropy                (300., 3. ), 0.392294792, 1e-5);
        eq_up_t <6> (massic_entropy                (300., 80.), 0.368563852, 1e-5);
        eq_up_t <6> (massic_entropy                (500., 3. ), 0.258041912e1, 1e-6);
        eq_up_t <6> (massic_isobaric_heat_capacity (300., 3. ), 0.417301218e1, 1e-6);
        eq_up_t <6> (massic_isobaric_heat_capacity (300., 80.), 0.401008987e1, 1e-6);
        eq_up_t <6> (massic_isobaric_heat_capacity (500., 3. ), 0.465580682e1, 1e-6);
        eq_up_t <6> (speed_of_sound                (300., 3. ), 0.150773921e4, 1e-6);
        eq_up_t <6> (speed_of_sound                (300., 80.), 0.163469054e4, 1e-6);
        eq_up_t <6> (speed_of_sound                (500., 3. ), 0.124071337e4, 1e-6);
    };
    SUBCASE("iapws-r7-region-1-backward")
    {
        // How do we deal with regions here?
        // A posteriori?
        /*
        eq_up_to <6> (temperature ( 3.,  500.), 0.391798509e3, 1e-6); 
        eq_up_to <6> (temperature (80.,  500.), 0.378108626e3, 1e-6); 
        eq_up_to <6> (temperature (80., 1500.), 0.611041229e3, 1e-6); 
        eq_up_to <6> (temperature ( 3., 0.5), 0.307842258e3, 1e-6); 
        eq_up_to <6> (temperature (80., 0.5), 0.309979785e3, 1e-6); 
        eq_up_to <6> (temperature (80., 3.0), 0.565899909e3, 1e-6);
        */
    };
    SUBCASE("iapws-r7-region-2")
    {
        eq_up_to <6> (massic_volume                 (300., 0.0035), 0.394913866e2, 1e-6);
        eq_up_to <6> (massic_volume                 (700., 0.0035), 0.923015898e2, 1e-6);
        eq_up_to <6> (massic_volume                 (700., 30.   ), 0.542946619e-2, 1e-6);
        eq_up_to <6> (massic_enthalpy               (300., 0.0035), 0.254991145e4, 1e-5);
        eq_up_to <6> (massic_enthalpy               (700., 0.0035), 0.333568375e4, 1e-5);
        eq_up_to <6> (massic_enthalpy               (700., 30.   ), 0.263149474e4, 1e-6);
        eq_up_to <6> (massic_internal_energy        (300., 0.0035), 0.241169160e4, 1e-5);
        eq_up_to <6> (massic_internal_energy        (700., 0.0035), 0.301262819e4, 1e-5);
        eq_up_to <6> (massic_internal_energy        (700., 30.   ), 0.246861076e4, 1e-6);
        eq_up_to <6> (massic_entropy                (300., 0.0035), 0.852238967e1, 1e-5);
        eq_up_to <6> (massic_entropy                (700., 0.0035), 0.101749996e2, 1e-5);
        eq_up_to <6> (massic_entropy                (700., 30.   ), 0.517540298e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (300., 0.0035), 0.191300162e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (700., 0.0035), 0.208141274e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (700., 30.   ), 0.103505092e2, 1e-6);
        eq_up_to <6> (speed_of_sound                (300., 0.0035), 0.427920172e3, 1e-6);
        eq_up_to <6> (speed_of_sound                (700., 0.0035), 0.644289068e3, 1e-6);
        eq_up_to <6> (speed_of_sound                (700., 30.   ), 0.480386523e3, 1e-6);
    };
    SUBCASE("iapws-r7-region-2-metastable-vapor")
    {
        eq_up_to <6> (massic_volume                 (450., 1.  ), 0.192516540, 1e-6);
        eq_up_to <6> (massic_volume                 (440., 1.  ), 0.186212297, 1e-6);
        eq_up_to <6> (massic_volume                 (450., 1.5),  0.121685206, 1e-6);
        eq_up_to <6> (massic_enthalpy               (450., 1.  ), 0.276881115e4, 1e-6);
        eq_up_to <6> (massic_enthalpy               (440., 1.  ), 0.274015123e4, 1e-6);
        eq_up_to <6> (massic_enthalpy               (450., 1.5),  0.272134539e4, 1e-6);
        eq_up_to <6> (massic_internal_energy        (450., 1.  ), 0.257629461e4, 1e-6);
        eq_up_to <6> (massic_internal_energy        (440., 1.  ), 0.255393894e4, 1e-6);
        eq_up_to <6> (massic_internal_energy        (450., 1.5),  0.253881758e4, 1e-6);
        eq_up_to <6> (massic_entropy                (450., 1.  ), 0.656660377e1, 1e-6);
        eq_up_to <6> (massic_entropy                (440., 1.  ), 0.650218759e1, 1e-6);
        eq_up_to <6> (massic_entropy                (450., 1.5),  0.629170440e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (450., 1.  ), 0.276349265e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (440., 1.  ), 0.298166443e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (450., 1.5),  0.362795578e1, 1e-6);
        eq_up_to <6> (speed_of_sound                (450., 1.  ), 0.498408101e3, 1e-6);
        eq_up_to <6> (speed_of_sound                (440., 1.  ), 0.489363295e3, 1e-6);
        eq_up_to <6> (speed_of_sound                (450., 1.5),  0.481941819e3, 1e-6);

    };
    SUBCASE("iapws-r7-region-2-backward")
    {
        // Region?
        /*
        eq_up_to <6> (temperature (0.001, 3000.), 0.534433241e3, 1e-6);
        eq_up_to <6> (temperature (3.,    3000.), 0.575373370e3, 1e-6);
        eq_up_to <6> (temperature (3.,    4000.), 0.101077577e4, 1e-6);
        eq_up_to <6> (temperature (5.,    3500.), 0.801299102e3, 1e-6);
        eq_up_to <6> (temperature (5.,    4000.), 0.101531583e4, 1e-6);
        eq_up_to <6> (temperature (25.,   3500.), 0.875279054e3, 1e-6);
        eq_up_to <6> (temperature (40.,   2700.), 0.743056411e3, 1e-6);
        eq_up_to <6> (temperature (60.,   2700.), 0.791137067e3, 1e-6);
        eq_up_to <6> (temperature (60.,   3200.), 0.882756860e3, 1e-6);
        eq_up_to <6> (temperature (0.1, 7.5 ), 0.399517097e3, 1e-6);
        eq_up_to <6> (temperature (0.1, 8.  ), 0.514127081e3, 1e-6);
        eq_up_to <6> (temperature (2.5, 8.  ), 0.103984917e4, 1e-6);
        eq_up_to <6> (temperature (8.,  6.  ), 0.600484040e3, 1e-6);
        eq_up_to <6> (temperature (8.,  7.5 ), 0.106495556e4, 1e-6);
        eq_up_to <6> (temperature (90., 6.  ), 0.103801126e4, 1e-6);
        eq_up_to <6> (temperature (20., 5.75), 0.697992849e3, 1e-6);
        eq_up_to <6> (temperature (80., 5.25), 0.854011484e3, 1e-6);
        eq_up_to <6> (temperature (80., 5.75), 0.949017998e3, 1e-6);
        */
    };
    SUBCASE("iapws-r7-region-3")
    {
        // Region?
        /*
        eq_up_to <6> (pressure                      (650., 500.), 0.255837018e2, 1e-6);
        eq_up_to <6> (pressure                      (650., 200.), 0.222930643e2, 1e-6);
        eq_up_to <6> (pressure                      (750., 500.), 0.783095639e2, 1e-6);
        eq_up_to <6> (massic_enthalpy               (650., 500.), 0.186343019e4, 1e-5);
        eq_up_to <6> (massic_enthalpy               (650., 200.), 0.237512401e4, 1e-5);
        eq_up_to <6> (massic_enthalpy               (750., 500.), 0.225868845e4, 1e-6);
        eq_up_to <6> (massic_internal_energy        (650., 500.), 0.181226279e4, 1e-5);
        eq_up_to <6> (massic_internal_energy        (650., 200.), 0.226365868e4, 1e-5);
        eq_up_to <6> (massic_internal_energy        (750., 500.), 0.210206932e4, 1e-6);
        eq_up_to <6> (massic_entropy                (650., 500.), 0.405427273e1, 1e-5);
        eq_up_to <6> (massic_entropy                (650., 200.), 0.485438792e1, 1e-5);
        eq_up_to <6> (massic_entropy                (750., 500.), 0.446971906e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (650., 500.), 0.138935717e2, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (650., 200.), 0.446579342e2, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (750., 500.), 0.634165359e1, 1e-6);
        eq_up_to <6> (speed_of_sound                (650., 500.), 0.502005554e3, 1e-6);
        eq_up_to <6> (speed_of_sound                (650., 200.), 0.383444594e3, 1e-6);
        eq_up_to <6> (speed_of_sound                (750., 500.), 0.760696041e3, 1e-6);
        */
    };
    SUBCASE("iapws-r7-region-4")
    {
        eq_up_to <6> (saturation_pressure (300.), 0.353658941e-2, 1e-6);
        eq_up_to <6> (saturation_pressure (500.), 0.263889776e1, 1e-6);
        eq_up_to <6> (saturation_pressure (600.), 0.123443146e2, 1e-6);
        eq_up_to <6> (saturation_temperature (0.1), 0.372755919e3, 1e-6);
        eq_up_to <6> (saturation_temperature (1. ), 0.453035632e3, 1e-6);
        eq_up_to <6> (saturation_temperature (10.), 0.584149488e3, 1e-6);
    };
    SUBCASE("iapws-r7-region-5")
    {
        eq_up_to <6> (massic_volume                 (1500., 0.5), 0.138455090e1, 1e-6);
        eq_up_to <6> (massic_volume                 (1500., 30.), 0.230761299e-1, 1e-6);
        eq_up_to <6> (massic_volume                 (2000., 30.), 0.311385219e-1, 1e-6);
        eq_up_to <6> (massic_enthalpy               (1500., 0.5), 0.521976855e4, 1e-6);
        eq_up_to <6> (massic_enthalpy               (1500., 30.), 0.516723514e4, 1e-6);
        eq_up_to <6> (massic_enthalpy               (2000., 30.), 0.657122604e4, 1e-6);
        eq_up_to <6> (massic_internal_energy        (1500., 0.5), 0.452749310e4, 1e-6);
        eq_up_to <6> (massic_internal_energy        (1500., 30.), 0.447495124e4, 1e-6);
        eq_up_to <6> (massic_internal_energy        (2000., 30.), 0.563707038e4, 1e-6);
        eq_up_to <6> (massic_entropy                (1500., 0.5), 0.965408875e1, 1e-6);
        eq_up_to <6> (massic_entropy                (1500., 30.), 0.772970133e1, 1e-6);
        eq_up_to <6> (massic_entropy                (2000., 30.), 0.853640523e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (1500., 0.5), 0.261609445e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (1500., 30.), 0.272724317e1, 1e-6);
        eq_up_to <6> (massic_isobaric_heat_capacity (2000., 30.), 0.288569882e1, 1e-6);
        eq_up_to <6> (speed_of_sound                (1500., 0.5), 0.917068690e3, 1e-6);
        eq_up_to <6> (speed_of_sound                (1500., 30.), 0.928548002e3, 1e-6);
        eq_up_to <6> (speed_of_sound                (2000., 30.), 0.106736948e4, 1e-6);

    };
#endif
}

