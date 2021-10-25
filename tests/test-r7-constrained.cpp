#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/r7.hpp"
    using namespace isto::iapws::r7;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"

    using namespace isto::units::unit::literals;

    constexpr auto
operator ""_kg_per_m3 (long double v)
{
    return v * unit::kilogram <> * pow <-3> (unit::metre <>);
}

    constexpr auto
eps = 1e-8;

TEST_CASE("r7.hpp (constrained)")
{
    SUBCASE("saturation line")
    {
        CHECK_U(saturation_pressure (300._K), 0.353658941e-2_MPa, eps);
        CHECK_U(saturation_pressure (500._K), 0.263889776e1_MPa, eps);
        CHECK_U(saturation_pressure (600._K), 0.123443146e2_MPa, eps);
        CHECK_U(saturation_temperature (0.1_MPa), 0.372755919e3_K, eps);
        CHECK_U(saturation_temperature (1._MPa ), 0.453035632e3_K, eps);
        CHECK_U(saturation_temperature (10._MPa), 0.584149488e3_K, eps);
        CHECK_US(saturation_pressure_h (massic_enthalpy_t { 1700e3 }), pressure_t { 1.724175718e7 }, 1e7,eps);
        CHECK_US(saturation_pressure_h (massic_enthalpy_t { 2000e3 }), pressure_t { 2.193442957e7 }, 1e7,eps);
        CHECK_US(saturation_pressure_h (massic_enthalpy_t { 2400e3 }), pressure_t { 2.018090839e7 }, 1e7,eps);
        CHECK_US(saturation_pressure_s (massic_entropy_t { 3.8e3 }), pressure_t { 1.687755057e7 }, 1e7,eps);
        CHECK_US(saturation_pressure_s (massic_entropy_t { 4.2e3 }), pressure_t { 2.164451789e7 }, 1e7,eps);
        CHECK_US(saturation_pressure_s (massic_entropy_t { 5.2e3 }), pressure_t { 1.668968482e7 }, 1e7,eps);
        CHECK_US(saturation_temperature_hs (massic_enthalpy_t { 1800e3 }, massic_entropy_t { 5.3e3 }), temperature_t { 3.468475498e2 }, 1e2,eps);
        CHECK_US(saturation_temperature_hs (massic_enthalpy_t { 2400e3 }, massic_entropy_t { 6.0e3 }), temperature_t { 4.251373305e2 }, 1e2,eps);
        CHECK_US(saturation_temperature_hs (massic_enthalpy_t { 2500e3 }, massic_entropy_t { 5.5e3 }), temperature_t { 5.225579013e2 }, 1e2,eps);
    }
    SUBCASE("iapws-r7-region-1")
    {
        CHECK_US(r1::massic_volume                 (3._MPa , 300._K), massic_volume_t        { 0.100215168e-2 }, 1e-2,  eps);
        CHECK_US(r1::massic_volume                 (80._MPa, 300._K), massic_volume_t        { 0.971180894e-3 }, 1e-2,  eps);
        CHECK_US(r1::massic_volume                 (3._MPa , 500._K), massic_volume_t        { 0.120241800e-2 }, 1e-2,  eps);
        CHECK_US(r1::massic_enthalpy               (3._MPa , 300._K), massic_enthalpy_t      { 0.115331273e3 * 1e3  }, 1e-6, eps);
        CHECK_US(r1::massic_enthalpy               (80._MPa, 300._K), massic_enthalpy_t      { 0.184142828e3 * 1e3  }, 1e-6, eps);
        CHECK_US(r1::massic_enthalpy               (3._MPa , 500._K), massic_enthalpy_t      { 0.975542239e3 * 1e3  }, 1e-6, eps);
        CHECK_US(r1::massic_internal_energy        (3._MPa , 300._K), massic_energy_t        { 0.112324818e3 * 1e3  }, 1e-6, eps);
        CHECK_US(r1::massic_internal_energy        (80._MPa, 300._K), massic_energy_t        { 0.106448356e3 * 1e3  }, 1e-6, eps);
        CHECK_US(r1::massic_internal_energy        (3._MPa , 500._K), massic_energy_t        { 0.971934985e3 * 1e3  }, 1e-6, eps);
        CHECK_US(r1::massic_entropy                (3._MPa , 300._K), massic_entropy_t       { 0.392294792   * 1e3  }, 1e-2, eps);
        CHECK_US(r1::massic_entropy                (80._MPa, 300._K), massic_entropy_t       { 0.368563852   * 1e3  }, 1e-3, eps);
        CHECK_US(r1::massic_entropy                (3._MPa , 500._K), massic_entropy_t       { 0.258041912e1 * 1e3  }, 1e-4, eps);
        CHECK_US(r1::massic_isobaric_heat_capacity (3._MPa , 300._K), massic_heat_capacity_t { 0.417301218e1 * 1e3  }, 1e-4, eps);
        CHECK_US(r1::massic_isobaric_heat_capacity (80._MPa, 300._K), massic_heat_capacity_t { 0.401008987e1 * 1e3  }, 1e-4, eps);
        CHECK_US(r1::massic_isobaric_heat_capacity (3._MPa , 500._K), massic_heat_capacity_t { 0.465580682e1 * 1e3  }, 1e-4, eps);
        CHECK_US(r1::speed_of_sound                (3._MPa , 300._K), velocity_t             { 0.150773921e4  }, 1e4, eps);
        CHECK_US(r1::speed_of_sound                (80._MPa, 300._K), velocity_t             { 0.163469054e4  }, 1e4, eps);
        CHECK_US(r1::speed_of_sound                (3._MPa , 500._K), velocity_t             { 0.124071337e4  }, 1e4, eps);
        CHECK_US(r1::temperature (pressure_t {  3e6 }, massic_enthalpy_t { 500e3 }), 0.391798509e3_K, 1e3, eps); 
        CHECK_US(r1::temperature (pressure_t { 80e6 }, massic_enthalpy_t { 500e3 }), 0.378108626e3_K, 1e3, eps); 
        CHECK_US(r1::temperature (pressure_t { 80e6 }, massic_enthalpy_t {1500e3 }), 0.611041229e3_K, 1e3, eps); 
        CHECK_US(r1::temperature (pressure_t {  3e6 }, massic_entropy_t  { 0.5e3 }), 0.307842258e3_K, 1e3, eps); 
        CHECK_US(r1::temperature (pressure_t { 80e6 }, massic_entropy_t  { 0.5e3 }), 0.309979785e3_K, 1e3, eps); 
        CHECK_US(r1::temperature (pressure_t { 80e6 }, massic_entropy_t  { 3.0e3 }), 0.565899909e3_K, 1e3, eps);
        CHECK_US(r1::pressure (massic_enthalpy_t {0.001e3 }, massic_entropy_t {   0. }), pressure_t { 9.800980612e-4  }, 1e4, eps);
        CHECK_US(r1::pressure (massic_enthalpy_t {   90e3 }, massic_entropy_t {   0. }), pressure_t { 9.192954727e1   }, 1e1, eps);
        CHECK_US(r1::pressure (massic_enthalpy_t { 1500e3 }, massic_entropy_t {3.4e3 }), pressure_t { 5.868294423e1   }, 1e1, eps);
    };
    SUBCASE("iapws-r7-region-2")
    {
        CHECK_US(r2::massic_volume                 (0.0035_MPa, 300._K), massic_volume_t        { 0.394913866e2         },  1e2, eps);
        CHECK_US(r2::massic_volume                 (0.0035_MPa, 700._K), massic_volume_t        { 0.923015898e2         },  1e2, eps);
        CHECK_US(r2::massic_volume                 (   30._MPa, 700._K), massic_volume_t        { 0.542946619e-2        }, 1e-2, eps);
        CHECK_US(r2::massic_enthalpy               (0.0035_MPa, 300._K), massic_enthalpy_t      { 0.254991145e4  * 1e3  },  1e7, eps);
        CHECK_US(r2::massic_enthalpy               (0.0035_MPa, 700._K), massic_enthalpy_t      { 0.333568375e4  * 1e3  },  1e7, eps);
        CHECK_US(r2::massic_enthalpy               (   30._MPa, 700._K), massic_enthalpy_t      { 0.263149474e4  * 1e3  },  1e7, eps);
        CHECK_US(r2::massic_internal_energy        (0.0035_MPa, 300._K), massic_energy_t        { 0.241169160e4  * 1e3  },  1e7, eps);
        CHECK_US(r2::massic_internal_energy        (0.0035_MPa, 700._K), massic_energy_t        { 0.301262819e4  * 1e3  },  1e7, eps);
        CHECK_US(r2::massic_internal_energy        (   30._MPa, 700._K), massic_energy_t        { 0.246861076e4  * 1e3  },  1e7, eps);
        CHECK_US(r2::massic_entropy                (0.0035_MPa, 300._K), massic_entropy_t       { 0.852238967e1  * 1e3  },  1e4, eps);
        CHECK_US(r2::massic_entropy                (0.0035_MPa, 700._K), massic_entropy_t       { 0.101749996e2  * 1e3  },  1e5, eps);
        CHECK_US(r2::massic_entropy                (   30._MPa, 700._K), massic_entropy_t       { 0.517540298e1  * 1e3  },  1e4, eps);
        CHECK_US(r2::massic_isobaric_heat_capacity (0.0035_MPa, 300._K), massic_heat_capacity_t { 0.191300162e1  * 1e3  },  1e4, eps);
        CHECK_US(r2::massic_isobaric_heat_capacity (0.0035_MPa, 700._K), massic_heat_capacity_t { 0.208141274e1  * 1e3  },  1e4, eps);
        CHECK_US(r2::massic_isobaric_heat_capacity (   30._MPa, 700._K), massic_heat_capacity_t { 0.103505092e2  * 1e3  },  1e5, eps);
        CHECK_US(r2::speed_of_sound                (0.0035_MPa, 300._K), velocity_t             { 0.427920172e3         },  1e3, eps);
        CHECK_US(r2::speed_of_sound                (0.0035_MPa, 700._K), velocity_t             { 0.644289068e3         },  1e3, eps);
        CHECK_US(r2::speed_of_sound                (   30._MPa, 700._K), velocity_t             { 0.480386523e3         },  1e3, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 2800e3 }, massic_entropy_t { 6.5e3 }), pressure_t { 1.371012767e6 }, 1e6, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 2800e3 }, massic_entropy_t { 9.5e3 }), pressure_t { 1.879743844e3 }, 1e9, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 4100e3 }, massic_entropy_t { 9.5e3 }), pressure_t { 1.024788997e5 }, 1e7, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 2800e3 }, massic_entropy_t { 6.0e3 }), pressure_t { 4.793911442e6 }, 1e6, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 3600e3 }, massic_entropy_t { 6.0e3 }), pressure_t { 8.395519209e7 }, 1e7, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 3600e3 }, massic_entropy_t { 7.0e3 }), pressure_t { 7.527161441e6 }, 1e6, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 2800e3 }, massic_entropy_t { 5.1e3 }), pressure_t { 9.439202060e7 }, 1e7, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 2800e3 }, massic_entropy_t { 5.8e3 }), pressure_t { 8.414574124e6 }, 1e6, eps);
        CHECK_US(r2::pressure (massic_enthalpy_t { 3400e3 }, massic_entropy_t { 5.8e3 }), pressure_t { 8.376903879e7 }, 1e7, eps);
        CHECK_US(r2::b2bc_p (pressure_t        { 0.100000000e9  }), massic_enthalpy_t { 0.3516004323e7 }, 1e7, eps);
        CHECK_US(r2::b2bc_h (massic_enthalpy_t { 0.3516004323e7 }), pressure_t        { 0.100000000e9 }, 1e9, eps);
        CHECK_US(r2::temperature (pressure_t { 0.001e6 }, massic_enthalpy_t { 3000e3 }), 0.534433241e3_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {     3e6 }, massic_enthalpy_t { 3000e3 }), 0.575373370e3_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {     3e6 }, massic_enthalpy_t { 4000e3 }), 0.101077577e4_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {     5e6 }, massic_enthalpy_t { 3500e3 }), 0.801299102e3_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {     5e6 }, massic_enthalpy_t { 4000e3 }), 0.101531583e4_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {    25e6 }, massic_enthalpy_t { 3500e3 }), 0.875279054e3_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {    40e6 }, massic_enthalpy_t { 2700e3 }), 0.743056411e3_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {    60e6 }, massic_enthalpy_t { 2700e3 }), 0.791137067e3_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {    60e6 }, massic_enthalpy_t { 3200e3 }), 0.882756860e3_K, 1e3, eps);
        CHECK_US(r2::temperature (pressure_t {   0.1e6 }, massic_entropy_t  {  7.5e3 }), 0.399517097e3_K, 1e3, 1e-6);
        CHECK_US(r2::temperature (pressure_t {   0.1e6 }, massic_entropy_t  {    8e3 }), 0.514127081e3_K, 1e3, 1e-6);
        CHECK_US(r2::temperature (pressure_t {   2.5e6 }, massic_entropy_t  {    8e3 }), 0.103984917e4_K, 1e4, 1e-6);
        CHECK_US(r2::temperature (pressure_t {     8e6 }, massic_entropy_t  {    6e3 }), 0.600484040e3_K, 1e3, 1e-6);
        CHECK_US(r2::temperature (pressure_t {     8e6 }, massic_entropy_t  {  7.5e3 }), 0.106495556e4_K, 1e4, 1e-6);
        CHECK_US(r2::temperature (pressure_t {    90e6 }, massic_entropy_t  {    6e3 }), 0.103801126e4_K, 1e4, 1e-6);
        CHECK_US(r2::temperature (pressure_t {    20e6 }, massic_entropy_t  { 5.75e3 }), 0.697992849e3_K, 1e3, 1e-6);
        CHECK_US(r2::temperature (pressure_t {    80e6 }, massic_entropy_t  { 5.25e3 }), 0.854011484e3_K, 1e3, 1e-6);
        CHECK_US(r2::temperature (pressure_t {    80e6 }, massic_entropy_t  { 5.75e3 }), 0.949017998e3_K, 1e3, 1e-6);
    };
    SUBCASE("iapws-r7-region-2-metastable-vapor")
    {
        CHECK_US(r2_metastable_vapor::massic_volume                 (1.0_MPa, 450._K), massic_volume_t        { 0.192516540          }, 1e0, eps);
        CHECK_US(r2_metastable_vapor::massic_volume                 (1.0_MPa, 440._K), massic_volume_t        { 0.186212297          }, 1e0, eps);
        CHECK_US(r2_metastable_vapor::massic_volume                 (1.5_MPa, 450._K), massic_volume_t        { 0.121685206          }, 1e0, eps);
        CHECK_US(r2_metastable_vapor::massic_enthalpy               (1.0_MPa, 450._K), massic_enthalpy_t      { 0.276881115e4 * 1e3  }, 1e7, eps);
        CHECK_US(r2_metastable_vapor::massic_enthalpy               (1.0_MPa, 440._K), massic_enthalpy_t      { 0.274015123e4 * 1e3  }, 1e7, eps);
        CHECK_US(r2_metastable_vapor::massic_enthalpy               (1.5_MPa, 450._K), massic_enthalpy_t      { 0.272134539e4 * 1e3  }, 1e7, eps);
        CHECK_US(r2_metastable_vapor::massic_internal_energy        (1.0_MPa, 450._K), massic_energy_t        { 0.257629461e4 * 1e3  }, 1e7, eps);
        CHECK_US(r2_metastable_vapor::massic_internal_energy        (1.0_MPa, 440._K), massic_energy_t        { 0.255393894e4 * 1e3  }, 1e7, eps);
        CHECK_US(r2_metastable_vapor::massic_internal_energy        (1.5_MPa, 450._K), massic_energy_t        { 0.253881758e4 * 1e3  }, 1e7, eps);
        CHECK_US(r2_metastable_vapor::massic_entropy                (1.0_MPa, 450._K), massic_entropy_t       { 0.656660377e1 * 1e3  }, 1e4, eps);
        CHECK_US(r2_metastable_vapor::massic_entropy                (1.0_MPa, 440._K), massic_entropy_t       { 0.650218759e1 * 1e3  }, 1e4, eps);
        CHECK_US(r2_metastable_vapor::massic_entropy                (1.5_MPa, 450._K), massic_entropy_t       { 0.629170440e1 * 1e3  }, 1e4, eps);
        CHECK_US(r2_metastable_vapor::massic_isobaric_heat_capacity (1.0_MPa, 450._K), massic_heat_capacity_t { 0.276349265e1 * 1e3  }, 1e4, eps);
        CHECK_US(r2_metastable_vapor::massic_isobaric_heat_capacity (1.0_MPa, 440._K), massic_heat_capacity_t { 0.298166443e1 * 1e3  }, 1e4, eps);
        CHECK_US(r2_metastable_vapor::massic_isobaric_heat_capacity (1.5_MPa, 450._K), massic_heat_capacity_t { 0.362795578e1 * 1e3  }, 1e4, eps);
        CHECK_US(r2_metastable_vapor::speed_of_sound                (1.0_MPa, 450._K), velocity_t             { 0.498408101e3        }, 1e3, eps);
        CHECK_US(r2_metastable_vapor::speed_of_sound                (1.0_MPa, 440._K), velocity_t             { 0.489363295e3        }, 1e3, eps);
        CHECK_US(r2_metastable_vapor::speed_of_sound                (1.5_MPa, 450._K), velocity_t             { 0.481941819e3        }, 1e3, eps);

    };
    SUBCASE("iapws-r7-region-3")
    {
        CHECK_US(r3::pressure                      (500._kg_per_m3, 650._K), pressure_t             { 0.255837018e2 * 1e6 }, 1e8, eps);
        CHECK_US(r3::pressure                      (200._kg_per_m3, 650._K), pressure_t             { 0.222930643e2 * 1e6 }, 1e8, eps);
        CHECK_US(r3::pressure                      (500._kg_per_m3, 750._K), pressure_t             { 0.783095639e2 * 1e6 }, 1e8, eps);
        CHECK_US(r3::massic_enthalpy               (500._kg_per_m3, 650._K), massic_enthalpy_t      { 0.186343019e4 * 1e3 }, 1e7, eps);
        CHECK_US(r3::massic_enthalpy               (200._kg_per_m3, 650._K), massic_enthalpy_t      { 0.237512401e4 * 1e3 }, 1e7, eps);
        CHECK_US(r3::massic_enthalpy               (500._kg_per_m3, 750._K), massic_enthalpy_t      { 0.225868845e4 * 1e3 }, 1e7, eps);
        CHECK_US(r3::massic_internal_energy        (500._kg_per_m3, 650._K), massic_energy_t        { 0.181226279e4 * 1e3 }, 1e7, eps);
        CHECK_US(r3::massic_internal_energy        (200._kg_per_m3, 650._K), massic_energy_t        { 0.226365868e4 * 1e3 }, 1e7, eps);
        CHECK_US(r3::massic_internal_energy        (500._kg_per_m3, 750._K), massic_energy_t        { 0.210206932e4 * 1e3 }, 1e7, eps);
        CHECK_US(r3::massic_entropy                (500._kg_per_m3, 650._K), massic_entropy_t       { 0.405427273e1 * 1e3 }, 1e4, eps);
        CHECK_US(r3::massic_entropy                (200._kg_per_m3, 650._K), massic_entropy_t       { 0.485438792e1 * 1e3 }, 1e4, eps);
        CHECK_US(r3::massic_entropy                (500._kg_per_m3, 750._K), massic_entropy_t       { 0.446971906e1 * 1e3 }, 1e4, eps);
        CHECK_US(r3::massic_isobaric_heat_capacity (500._kg_per_m3, 650._K), massic_heat_capacity_t { 0.138935717e2 * 1e3 }, 1e6, eps);
        CHECK_US(r3::massic_isobaric_heat_capacity (200._kg_per_m3, 650._K), massic_heat_capacity_t { 0.446579342e2 * 1e3 }, 1e5, eps);
        CHECK_US(r3::massic_isobaric_heat_capacity (500._kg_per_m3, 750._K), massic_heat_capacity_t { 0.634165359e1 * 1e3 }, 1e4, eps);
        CHECK_US(r3::speed_of_sound                (500._kg_per_m3, 650._K), velocity_t             { 0.502005554e3       }, 1e3, eps);
        CHECK_US(r3::speed_of_sound                (200._kg_per_m3, 650._K), velocity_t             { 0.383444594e3       }, 1e3, eps);
        CHECK_US(r3::speed_of_sound                (500._kg_per_m3, 750._K), velocity_t             { 0.760696041e3       }, 1e3, eps);
        CHECK_US(r3::temperature (pressure_t {  20e6 }, massic_enthalpy_t {  1700e3 }), temperature_t { 6.293083892e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t {  50e6 }, massic_enthalpy_t {  2000e3 }), temperature_t { 6.905718338e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t { 100e6 }, massic_enthalpy_t {  2100e3 }), temperature_t { 7.336163014e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t {  20e6 }, massic_enthalpy_t {  2500e3 }), temperature_t { 6.418418053e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t {  50e6 }, massic_enthalpy_t {  2400e3 }), temperature_t { 7.351848618e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t { 100e6 }, massic_enthalpy_t {  2700e3 }), temperature_t { 8.420460876e2 }, 1e2, eps);
        CHECK_US(r3::massic_volume (pressure_t {  20e6 }, massic_enthalpy_t {  1700e3 }), massic_volume_t { 1.749903962e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t {  50e6 }, massic_enthalpy_t {  2000e3 }), massic_volume_t { 1.908139035e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t { 100e6 }, massic_enthalpy_t {  2100e3 }), massic_volume_t { 1.676229776e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t {  20e6 }, massic_enthalpy_t {  2500e3 }), massic_volume_t { 6.670547043e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t {  50e6 }, massic_enthalpy_t {  2400e3 }), massic_volume_t { 2.801244590e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t { 100e6 }, massic_enthalpy_t {  2700e3 }), massic_volume_t { 2.404234998e-3 }, 1e-3, eps);
        CHECK_US(r3::temperature (pressure_t {  20e6 }, massic_entropy_t {  3.8e3 }), temperature_t { 6.282959869e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t {  50e6 }, massic_entropy_t {  3.6e3 }), temperature_t { 6.297158726e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t { 100e6 }, massic_entropy_t {  4.0e3 }), temperature_t { 7.056880237e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t {  20e6 }, massic_entropy_t {  5.0e3 }), temperature_t { 6.401176443e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t {  50e6 }, massic_entropy_t {  4.5e3 }), temperature_t { 7.163687517e2 }, 1e2, eps);
        CHECK_US(r3::temperature (pressure_t { 100e6 }, massic_entropy_t {  5.0e3 }), temperature_t { 8.474332825e2 }, 1e2, eps);
        CHECK_US(r3::massic_volume (pressure_t {  20e6 }, massic_entropy_t {   3.8e3 }), massic_volume_t { 1.733791463e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t {  50e6 }, massic_entropy_t {   3.6e3 }), massic_volume_t { 1.469680170e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t { 100e6 }, massic_entropy_t {   4.0e3 }), massic_volume_t { 1.555893131e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t {  20e6 }, massic_entropy_t {   5.0e3 }), massic_volume_t { 6.262101987e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t {  50e6 }, massic_entropy_t {   4.5e3 }), massic_volume_t { 2.332634294e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (pressure_t { 100e6 }, massic_entropy_t {   5.0e3 }), massic_volume_t { 2.449610757e-3 }, 1e-3, eps);
        CHECK_US(r3::saturation_pressure (massic_enthalpy_t { 1700e3 }), pressure_t { 1.724175718e7 }, 1e7, eps);
        CHECK_US(r3::saturation_pressure (massic_enthalpy_t { 2000e3 }), pressure_t { 2.193442957e7 }, 1e7, eps);
        CHECK_US(r3::saturation_pressure (massic_enthalpy_t { 2400e3 }), pressure_t { 2.018090839e7 }, 1e7, eps);
        CHECK_US(r3::saturation_pressure (massic_entropy_t { 3.8e3 }), pressure_t { 1.687755057e7 }, 1e7, eps);
        CHECK_US(r3::saturation_pressure (massic_entropy_t { 4.2e3 }), pressure_t { 2.164451789e7 }, 1e7, eps);
        CHECK_US(r3::saturation_pressure (massic_entropy_t { 5.2e3 }), pressure_t { 1.668968482e7 }, 1e7, eps);
        CHECK_US(r3::pressure (massic_enthalpy_t { 1700e3 }, massic_entropy_t { 3.8e3 }), pressure_t { 2.555703246e7 }, 1e7, eps);
        CHECK_US(r3::pressure (massic_enthalpy_t { 2000e3 }, massic_entropy_t { 4.2e3 }), pressure_t { 4.540873468e7 }, 1e7, eps);
        CHECK_US(r3::pressure (massic_enthalpy_t { 2100e3 }, massic_entropy_t { 4.3e3 }), pressure_t { 6.078123340e7 }, 1e7, eps);
        CHECK_US(r3::pressure (massic_enthalpy_t { 2600e3 }, massic_entropy_t { 5.1e3 }), pressure_t { 3.434999263e7 }, 1e7, eps);
        CHECK_US(r3::pressure (massic_enthalpy_t { 2400e3 }, massic_entropy_t { 4.7e3 }), pressure_t { 6.363924887e7 }, 1e7, eps);
        CHECK_US(r3::pressure (massic_enthalpy_t { 2700e3 }, massic_entropy_t { 5.0e3 }), pressure_t { 8.839043281e7 }, 1e7, eps);
        CHECK_US(r3::saturation_temperature (massic_enthalpy_t { 1800e3 }, massic_entropy_t { 5.3e3 }), temperature_t { 3.468475498e2 }, 1e2, eps);
        CHECK_US(r3::saturation_temperature (massic_enthalpy_t { 2400e3 }, massic_entropy_t { 6.0e3 }), temperature_t { 4.251373305e2 }, 1e2, eps);
        CHECK_US(r3::saturation_temperature (massic_enthalpy_t { 2500e3 }, massic_entropy_t { 5.5e3 }), temperature_t { 5.225579013e2 }, 1e2, eps);
        CHECK_US(r3::massic_volume (50.0_MPa,   630.00_K), massic_volume_t { 1.470853100e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (80.0_MPa,   670.00_K), massic_volume_t { 1.503831359e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (50.0_MPa,   710.00_K), massic_volume_t { 2.204728587e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (80.0_MPa,   750.00_K), massic_volume_t { 1.973692940e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (20.0_MPa,   630.00_K), massic_volume_t { 1.761696406e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (30.0_MPa,   650.00_K), massic_volume_t { 1.819560617e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (26.0_MPa,   656.00_K), massic_volume_t { 2.245587720e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (30.0_MPa,   670.00_K), massic_volume_t { 2.506897702e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (26.0_MPa,   661.00_K), massic_volume_t { 2.970225962e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (30.0_MPa,   675.00_K), massic_volume_t { 3.004627086e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (26.0_MPa,   671.00_K), massic_volume_t { 5.019029401e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (30.0_MPa,   690.00_K), massic_volume_t { 4.656470142e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (23.6_MPa,   649.00_K), massic_volume_t { 2.163198378e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (24.0_MPa,   650.00_K), massic_volume_t { 2.166044161e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (23.6_MPa,   652.00_K), massic_volume_t { 2.651081407e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (24.0_MPa,   654.00_K), massic_volume_t { 2.967802335e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (23.6_MPa,   653.00_K), massic_volume_t { 3.273916816e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (24.0_MPa,   655.00_K), massic_volume_t { 3.550329864e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (23.5_MPa,   655.00_K), massic_volume_t { 4.545001142e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (24.0_MPa,   660.00_K), massic_volume_t { 5.100267704e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (23.0_MPa,   660.00_K), massic_volume_t { 6.109525997e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (24.0_MPa,   670.00_K), massic_volume_t { 6.427325645e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.6_MPa,   646.00_K), massic_volume_t { 2.117860851e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (23.0_MPa,   646.00_K), massic_volume_t { 2.062374674e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.6_MPa,   648.60_K), massic_volume_t { 2.533063780e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.8_MPa,   649.30_K), massic_volume_t { 2.572971781e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.6_MPa,   649.00_K), massic_volume_t { 2.923432711e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.8_MPa,   649.70_K), massic_volume_t { 2.913311494e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.6_MPa,   649.10_K), massic_volume_t { 3.131208996e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.8_MPa,   649.90_K), massic_volume_t { 3.221160278e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.6_MPa,   649.40_K), massic_volume_t { 3.715596186e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.8_MPa,   650.20_K), massic_volume_t { 3.664754790e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (21.1_MPa,   640.00_K), massic_volume_t { 1.970999272e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (21.8_MPa,   643.00_K), massic_volume_t { 2.043919161e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (21.1_MPa,   644.00_K), massic_volume_t { 5.251009921e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (21.8_MPa,   648.00_K), massic_volume_t { 5.256844741e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (19.1_MPa,   635.00_K), massic_volume_t { 1.932829079e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (20.0_MPa,   638.00_K), massic_volume_t { 1.985387227e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (17.0_MPa,   626.00_K), massic_volume_t { 8.483262001e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (20.0_MPa,   640.00_K), massic_volume_t { 6.227528101e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (21.5_MPa,   644.60_K), massic_volume_t { 2.268366647e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.0_MPa,   646.10_K), massic_volume_t { 2.296350553e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.5_MPa,   648.60_K), massic_volume_t { 2.832373260e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.3_MPa,   647.90_K), massic_volume_t { 2.811424405e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.15_MPa,  647.50_K), massic_volume_t { 3.694032281e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.3_MPa,   648.10_K), massic_volume_t { 3.622226305e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.11_MPa,  648.00_K), massic_volume_t { 4.528072649e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.3_MPa,   649.00_K), massic_volume_t { 4.556905799e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.0_MPa,   646.84_K), massic_volume_t { 2.698354719e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.064_MPa, 647.05_K), massic_volume_t { 2.717655648e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.0_MPa,   646.89_K), massic_volume_t { 3.798732962e-3 }, 1e-3, eps);
        CHECK_US(r3::massic_volume (22.064_MPa, 647.15_K), massic_volume_t { 3.701940010e-3 }, 1e-3, eps);
    };
    /*
    SUBCASE("iapws-r7-region-4")
    {
            using namespace isto::iapws::r7::r4;
        CHECK_US(saturation_pressure    (300._K)  , 0.353658941e-2_MPa, 1e-2, eps);
        CHECK_US(saturation_pressure    (500._K)  , 0.263889776e1_MPa,  1e1, eps);
        CHECK_US(saturation_pressure    (600._K)  , 0.123443146e2_MPa,  1e2, eps);
        CHECK_US(saturation_temperature ( 0.1_MPa), 0.372755919e3_K,    1e3, eps);
        CHECK_US(saturation_temperature (  1._MPa), 0.453035632e3_K,    1e3, eps);
        CHECK_US(saturation_temperature ( 10._MPa), 0.584149488e3_K,    1e3, eps);
    };
    */
    SUBCASE("iapws-r7-region-5")
    {
        CHECK_US(r5::massic_volume                 (0.5_MPa, 1500._K), massic_volume_t        { 0.138455090e1         }, 1e1,  eps);
        CHECK_US(r5::massic_volume                 (30._MPa, 1500._K), massic_volume_t        { 0.230761299e-1        }, 1e-1, eps);
        CHECK_US(r5::massic_volume                 (30._MPa, 2000._K), massic_volume_t        { 0.311385219e-1        }, 1e-1, eps);
        CHECK_US(r5::massic_enthalpy               (0.5_MPa, 1500._K), massic_enthalpy_t      { 0.521976855e4  * 1e3  }, 1e7,  eps);
        CHECK_US(r5::massic_enthalpy               (30._MPa, 1500._K), massic_enthalpy_t      { 0.516723514e4  * 1e3  }, 1e7,  eps);
        CHECK_US(r5::massic_enthalpy               (30._MPa, 2000._K), massic_enthalpy_t      { 0.657122604e4  * 1e3  }, 1e7,  eps);
        CHECK_US(r5::massic_internal_energy        (0.5_MPa, 1500._K), massic_energy_t        { 0.452749310e4  * 1e3  }, 1e7,  eps);
        CHECK_US(r5::massic_internal_energy        (30._MPa, 1500._K), massic_energy_t        { 0.447495124e4  * 1e3  }, 1e7,  eps);
        CHECK_US(r5::massic_internal_energy        (30._MPa, 2000._K), massic_energy_t        { 0.563707038e4  * 1e3  }, 1e7,  eps);
        CHECK_US(r5::massic_entropy                (0.5_MPa, 1500._K), massic_entropy_t       { 0.965408875e1  * 1e3  }, 1e4,  eps);
        CHECK_US(r5::massic_entropy                (30._MPa, 1500._K), massic_entropy_t       { 0.772970133e1  * 1e3  }, 1e4,  eps);
        CHECK_US(r5::massic_entropy                (30._MPa, 2000._K), massic_entropy_t       { 0.853640523e1  * 1e3  }, 1e4,  eps);
        CHECK_US(r5::massic_isobaric_heat_capacity (0.5_MPa, 1500._K), massic_heat_capacity_t { 0.261609445e1  * 1e3  }, 1e4,  eps);
        CHECK_US(r5::massic_isobaric_heat_capacity (30._MPa, 1500._K), massic_heat_capacity_t { 0.272724317e1  * 1e3  }, 1e4,  eps);
        CHECK_US(r5::massic_isobaric_heat_capacity (30._MPa, 2000._K), massic_heat_capacity_t { 0.288569882e1  * 1e3  }, 1e4,  eps);
        CHECK_US(r5::speed_of_sound                (0.5_MPa, 1500._K), velocity_t             { 0.917068690e3         }, 1e3,  eps);
        CHECK_US(r5::speed_of_sound                (30._MPa, 1500._K), velocity_t             { 0.928548002e3         }, 1e3,  eps);
        CHECK_US(r5::speed_of_sound                (30._MPa, 2000._K), velocity_t             { 0.106736948e4         }, 1e4,  eps);
    };
    SUBCASE("region P-T")
    {
        CHECK_U(b23 (0.62315e3_K),       0.165291643e2_MPa, eps);
        CHECK_U(b23 (0.165291643e2_MPa), 0.62315e3_K,       eps);
        CHECK(region (50._MPa, 280._K)  == 1);
        CHECK(region (50._MPa, 1070._K) == 2);
        CHECK(region (50._MPa, 630._K)  == 3);
        //CHECK(region (0.353658941e-2_MPa, 300._K) == 4);
        CHECK(region (10._MPa, 1100._K)  == 5);
    }
    SUBCASE("region H-S")
    {
        CHECK(region (massic_enthalpy_t { 0.001e3 }, massic_entropy_t {    0. }) == 1);
        CHECK(region (massic_enthalpy_t {    90e3 }, massic_entropy_t {    0. }) == 1);
        CHECK(region (massic_enthalpy_t {  1500e3 }, massic_entropy_t { 3.4e3 }) == 1);
        CHECK(region (massic_enthalpy_t {  2800e3 }, massic_entropy_t { 6.5e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  2800e3 }, massic_entropy_t { 9.5e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  4100e3 }, massic_entropy_t { 9.5e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  2800e3 }, massic_entropy_t { 6.0e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  3600e3 }, massic_entropy_t { 6.0e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  3600e3 }, massic_entropy_t { 7.0e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  2800e3 }, massic_entropy_t { 5.1e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  2800e3 }, massic_entropy_t { 5.8e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  3400e3 }, massic_entropy_t { 5.8e3 }) == 2);
        CHECK(region (massic_enthalpy_t {  1700e3 }, massic_entropy_t { 3.8e3 }) == 3);
        CHECK(region (massic_enthalpy_t {  2000e3 }, massic_entropy_t { 4.2e3 }) == 3);
        CHECK(region (massic_enthalpy_t {  2100e3 }, massic_entropy_t { 4.3e3 }) == 3);
        CHECK(region (massic_enthalpy_t {  2600e3 }, massic_entropy_t { 5.1e3 }) == 3);
        CHECK(region (massic_enthalpy_t {  2400e3 }, massic_entropy_t { 4.7e3 }) == 3);
        CHECK(region (massic_enthalpy_t {  2700e3 }, massic_entropy_t { 5.0e3 }) == 3);
        CHECK(region (massic_enthalpy_t {  1800e3 }, massic_entropy_t { 5.3e3 }) == 4);
        CHECK(region (massic_enthalpy_t {  2400e3 }, massic_entropy_t { 6.0e3 }) == 4);
        CHECK(region (massic_enthalpy_t {  2500e3 }, massic_entropy_t { 5.5e3 }) == 4);
    }
    SUBCASE("region P-H")
    {
        CHECK(region (pressure_t {     3e6 }, massic_enthalpy_t {  500e3 }) == 1);
        CHECK(region (pressure_t {    80e6 }, massic_enthalpy_t {  500e3 }) == 1);
        CHECK(region (pressure_t {    80e6 }, massic_enthalpy_t { 1500e3 }) == 1);
        CHECK(region (pressure_t { 0.001e6 }, massic_enthalpy_t { 3000e3 }) == 2);
        CHECK(region (pressure_t {     3e6 }, massic_enthalpy_t { 3000e3 }) == 2);
        CHECK(region (pressure_t {     3e6 }, massic_enthalpy_t { 4000e3 }) == 2);
        CHECK(region (pressure_t {     5e6 }, massic_enthalpy_t { 3500e3 }) == 2);
        CHECK(region (pressure_t {     5e6 }, massic_enthalpy_t { 4000e3 }) == 2);
        CHECK(region (pressure_t {    25e6 }, massic_enthalpy_t { 3500e3 }) == 2);
        CHECK(region (pressure_t {    40e6 }, massic_enthalpy_t { 2700e3 }) == 2);
        CHECK(region (pressure_t {    60e6 }, massic_enthalpy_t { 2700e3 }) == 2);
        CHECK(region (pressure_t {    60e6 }, massic_enthalpy_t { 3200e3 }) == 2);
        CHECK(region (pressure_t {    20e6 }, massic_enthalpy_t { 1700e3 }) == 3);
        CHECK(region (pressure_t {    50e6 }, massic_enthalpy_t { 2000e3 }) == 3);
        CHECK(region (pressure_t {   100e6 }, massic_enthalpy_t { 2100e3 }) == 3);
        CHECK(region (pressure_t {    20e6 }, massic_enthalpy_t { 2500e3 }) == 3);
        CHECK(region (pressure_t {    50e6 }, massic_enthalpy_t { 2400e3 }) == 3);
        CHECK(region (pressure_t {   100e6 }, massic_enthalpy_t { 2700e3 }) == 3);
        // 4?
    }
    SUBCASE("region P-S")
    {
        CHECK(region (pressure_t {     3e6 }, massic_entropy_t {  0.5e3 }) == 1);
        CHECK(region (pressure_t {    80e6 }, massic_entropy_t {  0.5e3 }) == 1);
        CHECK(region (pressure_t {    80e6 }, massic_entropy_t {  3.0e3 }) == 1);
        CHECK(region (pressure_t {   0.1e6 }, massic_entropy_t {  7.5e3 }) == 2);
        CHECK(region (pressure_t {   0.1e6 }, massic_entropy_t {    8e3 }) == 2);
        CHECK(region (pressure_t {   2.5e6 }, massic_entropy_t {    8e3 }) == 2);
        CHECK(region (pressure_t {     8e6 }, massic_entropy_t {    6e3 }) == 2);
        CHECK(region (pressure_t {     8e6 }, massic_entropy_t {  7.5e3 }) == 2);
        CHECK(region (pressure_t {    90e6 }, massic_entropy_t {    6e3 }) == 2);
        CHECK(region (pressure_t {    20e6 }, massic_entropy_t { 5.75e3 }) == 2);
        CHECK(region (pressure_t {    80e6 }, massic_entropy_t { 5.25e3 }) == 2);
        CHECK(region (pressure_t {    80e6 }, massic_entropy_t { 5.75e3 }) == 2);
        CHECK(region (pressure_t {    20e6 }, massic_entropy_t {  3.8e3 }) == 3);
        CHECK(region (pressure_t {    50e6 }, massic_entropy_t {  3.6e3 }) == 3);
        CHECK(region (pressure_t {   100e6 }, massic_entropy_t {  4.0e3 }) == 3);
        CHECK(region (pressure_t {    20e6 }, massic_entropy_t {  5.0e3 }) == 3);
        CHECK(region (pressure_t {    50e6 }, massic_entropy_t {  4.5e3 }) == 3);
        CHECK(region (pressure_t {   100e6 }, massic_entropy_t {  5.0e3 }) == 3);
        // 4?
    }
    SUBCASE("main API")
    {
        auto const t = temperature_t { 300. };
        auto const p = pressure_t    { 3e6 };
        auto const h = massic_enthalpy_t { 1e6 };
        auto const s = massic_entropy_t  { 1e3 };
        // pt
        massic_volume                 (p, t);
        density                       (p, t);
        massic_enthalpy               (p, t);
        massic_internal_energy        (p, t);
        massic_entropy                (p, t);
        massic_isobaric_heat_capacity (p, t);
        speed_of_sound                (p, t);
        // hs
        pressure                      (h, s);
        temperature                   (h, s);
        massic_volume                 (h, s);
        density                       (h, s);
        massic_internal_energy        (h, s);
        massic_isobaric_heat_capacity (h, s);
        speed_of_sound                (h, s);
        // ph
        temperature                   (p, h);
        massic_volume                 (p, h);
        density                       (p, h);
        massic_internal_energy        (p, h);
        massic_entropy                (p, h);
        massic_isobaric_heat_capacity (p, h);
        speed_of_sound                (p, h);
        // ps
        temperature                   (p, s);
        massic_volume                 (p, s);
        density                       (p, s);
        massic_enthalpy               (p, s);
        massic_internal_energy        (p, s);
        massic_isobaric_heat_capacity (p, s);
        speed_of_sound                (p, s);

        // tp
        massic_volume                 (t, p);
        density                       (t, p);
        massic_enthalpy               (t, p);
        massic_internal_energy        (t, p);
        massic_entropy                (t, p);
        massic_isobaric_heat_capacity (t, p);
        speed_of_sound                (t, p);
        // sh
        pressure                      (s, h);
        temperature                   (s, h);
        massic_volume                 (s, h);
        density                       (s, h);
        massic_internal_energy        (s, h);
        massic_isobaric_heat_capacity (s, h);
        speed_of_sound                (s, h);
        // hp
        temperature                   (h, p);
        massic_volume                 (h, p);
        density                       (h, p);
        massic_internal_energy        (h, p);
        massic_entropy                (h, p);
        massic_isobaric_heat_capacity (h, p);
        speed_of_sound                (h, p);
        // sp
        temperature                   (s, p);
        massic_volume                 (s, p);
        density                       (s, p);
        massic_enthalpy               (s, p);
        massic_internal_energy        (s, p);
        massic_isobaric_heat_capacity (s, p);
        speed_of_sound                (s, p);

        CHECK_US(massic_enthalpy                  (pressure_t { 0.0035e6 }, temperature_t { 300.   }), massic_enthalpy_t { 0.254991145e4 * 1e3 }, 1e4, eps);
        CHECK_US(massic_enthalpy                  (pressure_t { 0.0035e6 }, temperature_t { 700.   }), massic_enthalpy_t { 0.333568375e4 * 1e3 }, 1e4, eps);
        CHECK_US(massic_enthalpy                  (pressure_t { 0.5e6    }, temperature_t { 1500.  }), massic_enthalpy_t { 0.521976855e4 * 1e3 }, 1e7, eps);
        CHECK_US(massic_enthalpy                  (pressure_t { 3. * 1e6 }, temperature_t { 300.   }), massic_enthalpy_t { 0.115331273e3 * 1e3 }, 1e6, eps);
        CHECK_US(massic_enthalpy                  (pressure_t { 3. * 1e6 }, temperature_t { 500.   }), massic_enthalpy_t { 0.975542239e3 * 1e3 }, 1e6, eps);
        CHECK_US(massic_enthalpy                  (pressure_t { 30.e6    }, temperature_t { 1500.  }), massic_enthalpy_t { 0.516723514e4 * 1e3 }, 1e7, eps);
        CHECK_US(massic_enthalpy                  (pressure_t { 30.e6    }, temperature_t { 2000.  }), massic_enthalpy_t { 0.657122604e4 * 1e3 }, 1e7, eps);
        CHECK_US(massic_enthalpy                  (pressure_t { 30e6     }, temperature_t { 700.   }), massic_enthalpy_t { 0.263149474e4 * 1e3 }, 1e4, eps);
        CHECK_US(massic_enthalpy                  (pressure_t { 80.* 1e6 }, temperature_t { 300.   }), massic_enthalpy_t { 0.184142828e3 * 1e3 }, 1e6, eps);
        CHECK_US(massic_entropy                   (pressure_t { 0.0035e6 }, temperature_t { 300.   }), massic_entropy_t { 0.852238967e1 * 1e3 }, 1e1, eps);
        CHECK_US(massic_entropy                   (pressure_t { 0.0035e6 }, temperature_t { 700.   }), massic_entropy_t { 0.101749996e2 * 1e3 }, 1e2, eps);
        CHECK_US(massic_entropy                   (pressure_t { 0.5e6    }, temperature_t { 1500.  }), massic_entropy_t { 0.965408875e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_entropy                   (pressure_t { 3. * 1e6 }, temperature_t { 300.   }), massic_entropy_t { 0.392294792   * 1e3 }, 1e2, eps);
        CHECK_US(massic_entropy                   (pressure_t { 3. * 1e6 }, temperature_t { 500.   }), massic_entropy_t { 0.258041912e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_entropy                   (pressure_t { 30.e6    }, temperature_t { 1500.  }), massic_entropy_t { 0.772970133e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_entropy                   (pressure_t { 30.e6    }, temperature_t { 2000.  }), massic_entropy_t { 0.853640523e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_entropy                   (pressure_t { 30e6     }, temperature_t { 700.   }), massic_entropy_t { 0.517540298e1 * 1e3 }, 1e1, eps);
        CHECK_US(massic_entropy                   (pressure_t { 80.* 1e6 }, temperature_t { 300.   }), massic_entropy_t { 0.368563852   * 1e3 }, 1e3, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 0.0035e6 }, temperature_t { 300.   }), massic_energy_t { 0.241169160e4 * 1e3 }, 1e4, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 0.0035e6 }, temperature_t { 700.   }), massic_energy_t { 0.301262819e4 * 1e3 }, 1e4, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 0.5e6    }, temperature_t { 1500.  }), massic_energy_t { 0.452749310e4 * 1e3 }, 1e7, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 3. * 1e6 }, temperature_t { 300.   }), massic_energy_t { 0.112324818e3 * 1e3 }, 1e6, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 3. * 1e6 }, temperature_t { 500.   }), massic_energy_t { 0.971934985e3 * 1e3 }, 1e6, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 30.e6    }, temperature_t { 1500.  }), massic_energy_t { 0.447495124e4 * 1e3 }, 1e7, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 30.e6    }, temperature_t { 2000.  }), massic_energy_t { 0.563707038e4 * 1e3 }, 1e7, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 30e6     }, temperature_t { 700.   }), massic_energy_t { 0.246861076e4 * 1e3 }, 1e4, eps);
        CHECK_US(massic_internal_energy           (pressure_t { 80.* 1e6 }, temperature_t { 300.   }), massic_energy_t { 0.106448356e3 * 1e3 }, 1e6, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 0.0035e6 }, temperature_t { 300.   }), massic_heat_capacity_t { 0.191300162e1 * 1e3 }, 1e1, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 0.0035e6 }, temperature_t { 700.   }), massic_heat_capacity_t { 0.208141274e1 * 1e3 }, 1e1, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 0.5e6    }, temperature_t { 1500.  }), massic_heat_capacity_t { 0.261609445e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 3. * 1e6 }, temperature_t { 300.   }), massic_heat_capacity_t { 0.417301218e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 3. * 1e6 }, temperature_t { 500.   }), massic_heat_capacity_t { 0.465580682e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 30.e6    }, temperature_t { 1500.  }), massic_heat_capacity_t { 0.272724317e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 30.e6    }, temperature_t { 2000.  }), massic_heat_capacity_t { 0.288569882e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 30e6     }, temperature_t { 700.   }), massic_heat_capacity_t { 0.103505092e2 * 1e3 }, 1e2, eps);
        CHECK_US(massic_isobaric_heat_capacity    (pressure_t { 80.* 1e6 }, temperature_t { 300.   }), massic_heat_capacity_t { 0.401008987e1 * 1e3 }, 1e4, eps);
        CHECK_US(massic_volume                    (pressure_t { 100e6    }, massic_enthalpy_t { 2100e3 }), massic_volume_t { 1.676229776e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 100e6    }, massic_enthalpy_t { 2700e3 }), massic_volume_t { 2.404234998e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 20e6     }, massic_enthalpy_t { 1700e3 }), massic_volume_t { 1.749903962e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 20e6     }, massic_enthalpy_t { 2500e3 }), massic_volume_t { 6.670547043e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 50e6     }, massic_enthalpy_t { 2000e3 }), massic_volume_t { 1.908139035e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 50e6     }, massic_enthalpy_t { 2400e3 }), massic_volume_t { 2.801244590e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 100e6    }, massic_entropy_t  { 4.0e3  }), massic_volume_t { 1.555893131e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 100e6    }, massic_entropy_t  { 5.0e3  }), massic_volume_t { 2.449610757e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 20e6     }, massic_entropy_t  { 3.8e3  }), massic_volume_t { 1.733791463e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 20e6     }, massic_entropy_t  { 5.0e3  }), massic_volume_t { 6.262101987e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 50e6     }, massic_entropy_t  { 3.6e3  }), massic_volume_t { 1.469680170e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 50e6     }, massic_entropy_t  { 4.5e3  }), massic_volume_t { 2.332634294e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 0.0035e6 }, temperature_t     { 300.   }), massic_volume_t { 0.394913866e2  }, 1e2, eps);
        CHECK_US(massic_volume                    (pressure_t { 0.0035e6 }, temperature_t     { 700.   }), massic_volume_t { 0.923015898e2  }, 1e2, eps);
        CHECK_US(massic_volume                    (pressure_t { 0.5e6    }, temperature_t     { 1500.  }), massic_volume_t { 0.138455090e1  }, 1e1, eps);
        CHECK_US(massic_volume                    (pressure_t { 17.0e6   }, temperature_t     { 626.00 }), massic_volume_t { 8.483262001e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 19.1e6   }, temperature_t     { 635.00 }), massic_volume_t { 1.932829079e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 20.0e6   }, temperature_t     { 630.00 }), massic_volume_t { 1.761696406e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 20.0e6   }, temperature_t     { 638.00 }), massic_volume_t { 1.985387227e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 20.0e6   }, temperature_t     { 640.00 }), massic_volume_t { 6.227528101e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 21.1e6   }, temperature_t     { 640.00 }), massic_volume_t { 1.970999272e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 21.1e6   }, temperature_t     { 644.00 }), massic_volume_t { 5.251009921e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 21.5e6   }, temperature_t     { 644.60 }), massic_volume_t { 2.268366647e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 21.8e6   }, temperature_t     { 643.00 }), massic_volume_t { 2.043919161e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 21.8e6   }, temperature_t     { 648.00 }), massic_volume_t { 5.256844741e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.064e6 }, temperature_t     { 647.05 }), massic_volume_t { 2.717655648e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.064e6 }, temperature_t     { 647.15 }), massic_volume_t { 3.701940010e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.0e6   }, temperature_t     { 646.10 }), massic_volume_t { 2.296350553e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.0e6   }, temperature_t     { 646.84 }), massic_volume_t { 2.698354719e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.0e6   }, temperature_t     { 646.89 }), massic_volume_t { 3.798732962e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.11e6  }, temperature_t     { 648.00 }), massic_volume_t { 4.528072649e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.15e6  }, temperature_t     { 647.50 }), massic_volume_t { 3.694032281e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.3e6   }, temperature_t     { 647.90 }), massic_volume_t { 2.811424405e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.3e6   }, temperature_t     { 648.10 }), massic_volume_t { 3.622226305e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.3e6   }, temperature_t     { 649.00 }), massic_volume_t { 4.556905799e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.5e6   }, temperature_t     { 648.60 }), massic_volume_t { 2.832373260e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.6e6   }, temperature_t     { 646.00 }), massic_volume_t { 2.117860851e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.6e6   }, temperature_t     { 648.60 }), massic_volume_t { 2.533063780e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.6e6   }, temperature_t     { 649.00 }), massic_volume_t { 2.923432711e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.6e6   }, temperature_t     { 649.10 }), massic_volume_t { 3.131208996e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.6e6   }, temperature_t     { 649.40 }), massic_volume_t { 3.715596186e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.8e6   }, temperature_t     { 649.30 }), massic_volume_t { 2.572971781e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.8e6   }, temperature_t     { 649.70 }), massic_volume_t { 2.913311494e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.8e6   }, temperature_t     { 649.90 }), massic_volume_t { 3.221160278e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 22.8e6   }, temperature_t     { 650.20 }), massic_volume_t { 3.664754790e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 23.0e6   }, temperature_t     { 646.00 }), massic_volume_t { 2.062374674e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 23.0e6   }, temperature_t     { 660.00 }), massic_volume_t { 6.109525997e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 23.5e6   }, temperature_t     { 655.00 }), massic_volume_t { 4.545001142e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 23.6e6   }, temperature_t     { 649.00 }), massic_volume_t { 2.163198378e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 23.6e6   }, temperature_t     { 652.00 }), massic_volume_t { 2.651081407e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 23.6e6   }, temperature_t     { 653.00 }), massic_volume_t { 3.273916816e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 24.0e6   }, temperature_t     { 650.00 }), massic_volume_t { 2.166044161e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 24.0e6   }, temperature_t     { 654.00 }), massic_volume_t { 2.967802335e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 24.0e6   }, temperature_t     { 655.00 }), massic_volume_t { 3.550329864e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 24.0e6   }, temperature_t     { 660.00 }), massic_volume_t { 5.100267704e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 24.0e6   }, temperature_t     { 670.00 }), massic_volume_t { 6.427325645e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 26.0e6   }, temperature_t     { 656.00 }), massic_volume_t { 2.245587720e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 26.0e6   }, temperature_t     { 661.00 }), massic_volume_t { 2.970225962e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 26.0e6   }, temperature_t     { 671.00 }), massic_volume_t { 5.019029401e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 3. * 1e6 }, temperature_t     { 300.   }), massic_volume_t { 0.100215168e-2 }, 1e-2, eps);
        CHECK_US(massic_volume                    (pressure_t { 3. * 1e6 }, temperature_t     { 500.   }), massic_volume_t { 0.120241800e-2 }, 1e-2, eps);
        CHECK_US(massic_volume                    (pressure_t { 30.0e6   }, temperature_t     { 650.00 }), massic_volume_t { 1.819560617e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 30.0e6   }, temperature_t     { 670.00 }), massic_volume_t { 2.506897702e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 30.0e6   }, temperature_t     { 675.00 }), massic_volume_t { 3.004627086e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 30.0e6   }, temperature_t     { 690.00 }), massic_volume_t { 4.656470142e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 30.e6    }, temperature_t     { 1500.  }), massic_volume_t { 0.230761299e-1 }, 1e-1, eps);
        CHECK_US(massic_volume                    (pressure_t { 30.e6    }, temperature_t     { 2000.  }), massic_volume_t { 0.311385219e-1 }, 1e-1, eps);
        CHECK_US(massic_volume                    (pressure_t { 30e6     }, temperature_t     { 700.   }), massic_volume_t { 0.542946619e-2 }, 1e-2, eps);
        CHECK_US(massic_volume                    (pressure_t { 50.0e6   }, temperature_t     { 630.00 }), massic_volume_t { 1.470853100e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 50.0e6   }, temperature_t     { 710.00 }), massic_volume_t { 2.204728587e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 80.* 1e6 }, temperature_t     { 300.   }), massic_volume_t { 0.971180894e-3 }, 1e-2, eps);
        CHECK_US(massic_volume                    (pressure_t { 80.0e6   }, temperature_t     { 670.00 }), massic_volume_t { 1.503831359e-3 }, 1e-3, eps);
        CHECK_US(massic_volume                    (pressure_t { 80.0e6   }, temperature_t     { 750.00 }), massic_volume_t { 1.973692940e-3 }, 1e-3, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 0.001e3  }, massic_entropy_t { 0.     }), pressure_t { 9.800980612e-4 }, 1e4, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 1500e3   }, massic_entropy_t { 3.4e3  }), pressure_t { 5.868294423e1  }, 1e1, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 1700e3   }, massic_entropy_t { 3.8e3  }), pressure_t { 2.555703246e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2000e3   }, massic_entropy_t { 4.2e3  }), pressure_t { 4.540873468e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2100e3   }, massic_entropy_t { 4.3e3  }), pressure_t { 6.078123340e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2400e3   }, massic_entropy_t { 4.7e3  }), pressure_t { 6.363924887e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2600e3   }, massic_entropy_t { 5.1e3  }), pressure_t { 3.434999263e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2700e3   }, massic_entropy_t { 5.0e3  }), pressure_t { 8.839043281e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2800e3   }, massic_entropy_t { 5.1e3  }), pressure_t { 9.439202060e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2800e3   }, massic_entropy_t { 5.8e3  }), pressure_t { 8.414574124e6  }, 1e6, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2800e3   }, massic_entropy_t { 6.0e3  }), pressure_t { 4.793911442e6  }, 1e6, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2800e3   }, massic_entropy_t { 6.5e3  }), pressure_t { 1.371012767e6  }, 1e6, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 2800e3   }, massic_entropy_t { 9.5e3  }), pressure_t { 1.879743844e3  }, 1e3, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 3400e3   }, massic_entropy_t { 5.8e3  }), pressure_t { 8.376903879e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 3600e3   }, massic_entropy_t { 6.0e3  }), pressure_t { 8.395519209e7  }, 1e7, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 3600e3   }, massic_entropy_t { 7.0e3  }), pressure_t { 7.527161441e6  }, 1e6, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 4100e3   }, massic_entropy_t { 9.5e3  }), pressure_t { 1.024788997e5  }, 1e5, eps);
        CHECK_US(pressure_hs                      (massic_enthalpy_t { 90e3     }, massic_entropy_t { 0.     }), pressure_t { 9.192954727e1  }, 1e1, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 0.0035e6 }, temperature_t { 300.   }), velocity_t { 0.427920172e3 }, 1e3, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 0.0035e6 }, temperature_t { 700.   }), velocity_t { 0.644289068e3 }, 1e3, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 0.5e6    }, temperature_t { 1500.  }), velocity_t { 0.917068690e3 }, 1e3, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 3. * 1e6 }, temperature_t { 300.   }), velocity_t { 0.150773921e4 }, 1e4, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 3. * 1e6 }, temperature_t { 500.   }), velocity_t { 0.124071337e4 }, 1e4, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 30.e6    }, temperature_t { 1500.  }), velocity_t { 0.928548002e3 }, 1e3, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 30.e6    }, temperature_t { 2000.  }), velocity_t { 0.106736948e4 }, 1e4, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 30e6     }, temperature_t { 700.   }), velocity_t { 0.480386523e3 }, 1e3, eps);
        CHECK_US(speed_of_sound                   (pressure_t { 80.* 1e6 }, temperature_t { 300.   }), velocity_t { 0.163469054e4 }, 1e4, eps);
        CHECK_US(temperature                      (pressure_t { 0.001e6  }, massic_enthalpy_t { 3000e3 }), temperature_t { 0.534433241e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 100e6    }, massic_enthalpy_t { 2100e3 }), temperature_t { 7.336163014e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 100e6    }, massic_enthalpy_t { 2700e3 }), temperature_t { 8.420460876e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 20e6     }, massic_enthalpy_t { 1700e3 }), temperature_t { 6.293083892e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 20e6     }, massic_enthalpy_t { 2500e3 }), temperature_t { 6.418418053e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 25e6     }, massic_enthalpy_t { 3500e3 }), temperature_t { 0.875279054e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 3e6      }, massic_enthalpy_t { 3000e3 }), temperature_t { 0.575373370e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 3e6      }, massic_enthalpy_t { 4000e3 }), temperature_t { 0.101077577e4 }, 1e4, eps);
        CHECK_US(temperature                      (pressure_t { 3e6      }, massic_enthalpy_t { 500e3  }), temperature_t { 0.391798509e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 40e6     }, massic_enthalpy_t { 2700e3 }), temperature_t { 0.743056411e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 50e6     }, massic_enthalpy_t { 2000e3 }), temperature_t { 6.905718338e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 50e6     }, massic_enthalpy_t { 2400e3 }), temperature_t { 7.351848618e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 5e6      }, massic_enthalpy_t { 3500e3 }), temperature_t { 0.801299102e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 5e6      }, massic_enthalpy_t { 4000e3 }), temperature_t { 0.101531583e4 }, 1e4, eps);
        CHECK_US(temperature                      (pressure_t { 60e6     }, massic_enthalpy_t { 2700e3 }), temperature_t { 0.791137067e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 60e6     }, massic_enthalpy_t { 3200e3 }), temperature_t { 0.882756860e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 80e6     }, massic_enthalpy_t { 1500e3 }), temperature_t { 0.611041229e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 80e6     }, massic_enthalpy_t { 500e3  }), temperature_t { 0.378108626e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 0.1e6    }, massic_entropy_t  { 7.5e3  }), temperature_t { 0.399517097e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 0.1e6    }, massic_entropy_t  { 8e3    }), temperature_t { 0.514127081e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 100e6    }, massic_entropy_t  { 4.0e3  }), temperature_t { 7.056880237e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 100e6    }, massic_entropy_t  { 5.0e3  }), temperature_t { 8.474332825e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 2.5e6    }, massic_entropy_t  { 8e3    }), temperature_t { 0.103984917e4 }, 1e4, eps);
        CHECK_US(temperature                      (pressure_t { 20e6     }, massic_entropy_t  { 3.8e3  }), temperature_t { 6.282959869e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 20e6     }, massic_entropy_t  { 5.0e3  }), temperature_t { 6.401176443e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 20e6     }, massic_entropy_t  { 5.75e3 }), temperature_t { 0.697992849e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 3e6      }, massic_entropy_t  { 0.5e3  }), temperature_t { 0.307842258e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 50e6     }, massic_entropy_t  { 3.6e3  }), temperature_t { 6.297158726e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 50e6     }, massic_entropy_t  { 4.5e3  }), temperature_t { 7.163687517e2 }, 1e2, eps);
        CHECK_US(temperature                      (pressure_t { 80e6     }, massic_entropy_t  { 0.5e3  }), temperature_t { 0.309979785e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 80e6     }, massic_entropy_t  { 3.0e3  }), temperature_t { 0.565899909e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 80e6     }, massic_entropy_t  { 5.25e3 }), temperature_t { 0.854011484e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 80e6     }, massic_entropy_t  { 5.75e3 }), temperature_t { 0.949017998e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 8e6      }, massic_entropy_t  { 6e3    }), temperature_t { 0.600484040e3 }, 1e3, eps);
        CHECK_US(temperature                      (pressure_t { 8e6      }, massic_entropy_t  { 7.5e3  }), temperature_t { 0.106495556e4 }, 1e4, eps);
        CHECK_US(temperature                      (pressure_t { 90e6     }, massic_entropy_t  { 6e3    }), temperature_t { 0.103801126e4 }, 1e4, eps);
    }
    SUBCASE("expansion, compressibility, etc.")
    {
            using namespace isto::iapws::r7::detail;
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
                r = region (pressure_t { p }, temperature_t { t });
                INFO ("P = ", p * 1e-5, ", T = ", t - 273.15, ", region = ", r);
                CHECK_US(isobaric_cubic_expansion_coefficient (pressure_t { p }, temperature_t { t }), thermal_expansion_t { a }, fabs (a), 1e-4);
                CHECK_US(isothermal_compressibility           (pressure_t { p }, temperature_t { t }), compressibility_t   { b },       b , 1e-4);
                CHECK_US(relative_pressure_coefficient        (pressure_t { p }, temperature_t { t }), thermal_expansion_t { c }, fabs (c), 1e-4);
                CHECK_US(isothermal_stress_coefficient        (pressure_t { p }, temperature_t { t }), density_t           { d }, fabs (d), 1e-4);
            }
        }
    }
}

