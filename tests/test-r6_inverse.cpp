#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/calculisto/iapws/r6_inverse.hpp"
    using namespace calculisto::iapws::r6_inverse;
#include "../include/calculisto/iapws/detail/data_for_the_tests.hpp"
    using namespace calculisto::iapws::r6::r6_95_2016::detail;

TEST_CASE("r6_inverse.hpp")
{
        using namespace calculisto::iapws;
    for(const auto& e: table_7)
    {
        INFO ("P= ", e.P, ", T= ", e.T);
        CHECK (density_pt (e.P, e.T) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        CHECK (density_tp (e.T, e.P) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        // With initial guess
        CHECK (density_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        CHECK (density_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.D }.scale (1e3).epsilon (1e-6));
        // And info
        CHECK (density_pt (e.P, e.T, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.D }.scale (1e3).epsilon (1e-6));
        CHECK (density_tp (e.T, e.P, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.D }.scale (1e3).epsilon (1e-6));
    }
    {
            const auto
        [ r, i ] = density_tp (300.0,  0.992418352e-1 * 1e6, info::convergence);
        CHECK(i.convergence.size () > 1);
        /*
        for (auto&& [ v, f, df ]: i.convergence)
        {
            MESSAGE (v, ", ", f, ", ", df);
        }
        */
    }
    for(const auto& e: table_7)
    {
        INFO ("P= ", e.P, ", T= ", e.T);
        CHECK (massic_isochoric_heat_capacity_pt (e.P, e.T) == Approx { e.Cv }.scale (1e3).epsilon (1e-6));
        CHECK (massic_isochoric_heat_capacity_tp (e.T, e.P) == Approx { e.Cv }.scale (1e3).epsilon (1e-6));
        CHECK (massic_isochoric_heat_capacity_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.Cv }.scale (1e3).epsilon (1e-6));
        CHECK (massic_isochoric_heat_capacity_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.Cv }.scale (1e3).epsilon (1e-6));
        CHECK (massic_isochoric_heat_capacity_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.Cv }.scale (1e3).epsilon (1e-6));
        CHECK (massic_isochoric_heat_capacity_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.Cv }.scale (1e3).epsilon (1e-6));
        CHECK (massic_isochoric_heat_capacity_pt (e.P, e.T, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.Cv }.scale (1e3).epsilon (1e-6));
        CHECK (massic_isochoric_heat_capacity_tp (e.T, e.P, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.Cv }.scale (1e3).epsilon (1e-6));

        CHECK (speed_of_sound_pt (e.P, e.T) == Approx { e.W }.scale (1e3).epsilon (1e-6));
        CHECK (speed_of_sound_tp (e.T, e.P) == Approx { e.W }.scale (1e3).epsilon (1e-6));
        CHECK (speed_of_sound_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.W }.scale (1e3).epsilon (1e-6));
        CHECK (speed_of_sound_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.W }.scale (1e3).epsilon (1e-6));
        CHECK (speed_of_sound_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.W }.scale (1e3).epsilon (1e-6));
        CHECK (speed_of_sound_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.W }.scale (1e3).epsilon (1e-6));
        CHECK (speed_of_sound_pt (e.P, e.T, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.W }.scale (1e3).epsilon (1e-6));
        CHECK (speed_of_sound_tp (e.T, e.P, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.W }.scale (1e3).epsilon (1e-6));

        CHECK (massic_entropy_pt (e.P, e.T) == Approx { e.S }.scale (1e3).epsilon (1e-6));
        CHECK (massic_entropy_tp (e.T, e.P) == Approx { e.S }.scale (1e3).epsilon (1e-6));
        CHECK (massic_entropy_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.S }.scale (1e3).epsilon (1e-6));
        CHECK (massic_entropy_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.S }.scale (1e3).epsilon (1e-6));
        CHECK (massic_entropy_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.S }.scale (1e3).epsilon (1e-6));
        CHECK (massic_entropy_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.S }.scale (1e3).epsilon (1e-6));
        CHECK (massic_entropy_pt (e.P, e.T, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.S }.scale (1e3).epsilon (1e-6));
        CHECK (massic_entropy_tp (e.T, e.P, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.S }.scale (1e3).epsilon (1e-6));

    }
    SUBCASE ("Saturation")
    {
        for (auto const& entry: r7::detail::table_1)
        {
                const auto
            temperature = entry.temperature_kelvin;
            INFO("Temperature = ", temperature);
                const auto
            expected_pressure = entry.saturation_pressure * 1e5;
                const auto
            expected_density_liquid = 1 / entry.massic_volume_liquid;
                const auto
            expected_density_gas = 1 / entry.massic_volume_gas;
                const auto
            pressure_inital_guess = r7::saturation_pressure_t (temperature);
                const auto
            density_liquid_inital_guess = r7::r1::density_pt (
                  pressure_inital_guess
                , temperature
            );
                const auto
            density_gas_inital_guess = r7::r2::density_pt (
                  pressure_inital_guess
                , temperature
            );

            try
            {
                    const auto
                [ p_s, d_l, d_g, info_data ] = saturation_pressure_t (
                      temperature
                    , pressure_inital_guess
                    , density_liquid_inital_guess
                    , density_gas_inital_guess
                    , {
                          .max_iter = 100
                        , .converged = [](
                              Eigen::Matrix <double, 3, 1> const& current
                            , Eigen::Matrix <double, 3, 1> const& past
                            , Eigen::Matrix <double, 3, 1> const& result
                          ){
                            return
                                   (past - current).norm () / current.norm ()
                                   < 1e-8
                                || result.squaredNorm () == 0
                            ;
                            }
                      }
                    , info::convergence
                );
                CHECK(info_data.converged == true);
                CHECK(p_s == Approx { expected_pressure }.scale (expected_pressure).epsilon (1e-3));
                CHECK(d_l == Approx { expected_density_liquid }.scale (expected_density_liquid).epsilon (1e-3));
                CHECK(d_g == Approx { expected_density_gas }.scale (expected_density_gas).epsilon (1e-3));
            }
            catch (...)
            {
                MESSAGE("Exception at T = ", temperature);
            }
        }
        for (auto const& entry: r7::detail::table_1)
        {
                const auto
            temperature = entry.temperature_kelvin;
            INFO("Temperature = ", temperature);
                const auto
            expected_pressure = entry.saturation_pressure * 1e5;
                const auto
            expected_density_liquid = 1 / entry.massic_volume_liquid;
                const auto
            expected_density_gas = 1 / entry.massic_volume_gas;
            try
            {
                    const auto
                [ p_s, d_l, d_g ] = saturation_pressure_t (
                      temperature
                );
                CHECK(p_s == Approx { expected_pressure }.scale (expected_pressure).epsilon (1e-3));
                CHECK(d_l == Approx { expected_density_liquid }.scale (expected_density_liquid).epsilon (1e-3));
                CHECK(d_g == Approx { expected_density_gas }.scale (expected_density_gas).epsilon (1e-3));
            }
            catch (...)
            {
                MESSAGE("Exception at T = ", temperature);
            }
        }
    }
} // TEST_CASE("r6_inverse.hpp")
