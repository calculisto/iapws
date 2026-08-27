#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/calculisto/iapws/r6_inverse.hpp"
    using namespace calculisto::thermodynamics::iapws::r6_inverse;
#include "../include/calculisto/iapws/detail/data_for_the_tests.hpp"
    using namespace calculisto::thermodynamics::iapws::r6;
#include <calculisto/finite_difference/finite_difference.hpp>
    using calculisto::finite_difference::central_finite_difference;

TEST_CASE("r6_inverse.hpp")
{
        using namespace calculisto::thermodynamics::iapws;
    for(const auto& e: r6::detail::table_7)
    {
        INFO ("P= ", e.P, ", T= ", e.T, ", D= ", e.D);

        // Density
            const auto
        D_pt = density_pt (e.P, e.T);
            const auto
        D_tp = density_tp (e.T, e.P);
        CHECK (D_tp == Approx { D_pt });
        CHECK (D_tp == Approx { e.D });
        // With initial guess
        CHECK (density_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.D });
        CHECK (density_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.D });
        // And info
        CHECK (density_pt (e.P, e.T, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.D });
        CHECK (density_tp (e.T, e.P, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.D });

        // Temperature
            const auto
        T_pd = temperature_pd (e.P, e.D);
            const auto
        T_dp = temperature_dp (e.D, e.P);
        CHECK (T_dp == Approx { T_pd });
        CHECK (T_dp == Approx { e.T });
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
    for(const auto& e: r6::detail::table_7)
    {
        INFO ("P= ", e.P, ", T= ", e.T);
        CHECK (massic_isochoric_heat_capacity_pt (e.P, e.T) == Approx { e.Cv });
        CHECK (massic_isochoric_heat_capacity_tp (e.T, e.P) == Approx { e.Cv });
        CHECK (massic_isochoric_heat_capacity_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.Cv });
        CHECK (massic_isochoric_heat_capacity_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.Cv });
        CHECK (massic_isochoric_heat_capacity_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.Cv });
        CHECK (massic_isochoric_heat_capacity_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.Cv });
        CHECK (massic_isochoric_heat_capacity_pt (e.P, e.T, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.Cv });
        CHECK (massic_isochoric_heat_capacity_tp (e.T, e.P, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.Cv });

        CHECK (speed_of_sound_pt (e.P, e.T) == Approx { e.W });
        CHECK (speed_of_sound_tp (e.T, e.P) == Approx { e.W });
        CHECK (speed_of_sound_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.W });
        CHECK (speed_of_sound_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.W });
        CHECK (speed_of_sound_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.W });
        CHECK (speed_of_sound_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.W });
        CHECK (speed_of_sound_pt (e.P, e.T, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.W });
        CHECK (speed_of_sound_tp (e.T, e.P, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.W });

        CHECK (massic_entropy_pt (e.P, e.T) == Approx { e.S });
        CHECK (massic_entropy_tp (e.T, e.P) == Approx { e.S });
        CHECK (massic_entropy_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.S });
        CHECK (massic_entropy_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.S });
        CHECK (massic_entropy_pt (e.P, e.T, r7::density_pt (e.P, e.T)) == Approx { e.S });
        CHECK (massic_entropy_tp (e.T, e.P, r7::density_pt (e.P, e.T)) == Approx { e.S });
        CHECK (massic_entropy_pt (e.P, e.T, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.S });
        CHECK (massic_entropy_tp (e.T, e.P, r7::density_pt (e.P, e.T), info::iterations).first == Approx { e.S });

    }
    SUBCASE ("Saturation")
    {
            using namespace r6::detail;

        for (auto i = 0u; i < table_13_1_liquid.size (); ++i)
        {
                const auto
            temperature = table_13_1_liquid.at (i).T;
            INFO("Temperature = ", temperature);
                const auto
            expected_pressure = table_13_1_liquid[i].P * 1e6;
                const auto
            expected_density_liquid = table_13_1_liquid[i].D;
                const auto
            expected_density_gas = table_13_1_gas.at (i).D;
            try
            {
                    const auto
                [ p_s, d_l, d_g ] = saturation_pressure_t (temperature);
                CHECK(p_s == Approx { expected_pressure }.scale (expected_pressure).epsilon (1e-3));
                CHECK(d_l == Approx { expected_density_liquid }.scale (expected_density_liquid).epsilon (1e-3));
                CHECK(d_g == Approx { expected_density_gas }.scale (expected_density_gas).epsilon (1e-3));
            }
            catch (...)
            {
                MESSAGE(
                      "Exception at T = "
                    , temperature
                    , ", we are probably too close to the critical point"
                );
            }
        }
        // FIXME: we need to test this against table 13.1 of Wagner et Pruss,
        // 2002 instead.
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
                // TODO: epsilons is probably too high here
                CHECK(p_s == Approx { expected_pressure }.scale (expected_pressure).epsilon (1e-0));
                CHECK(d_l == Approx { expected_density_liquid }.scale (expected_density_liquid).epsilon (1e-0));
                CHECK(d_g == Approx { expected_density_gas }.scale (expected_density_gas).epsilon (1e-0));
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
    SUBCASE("R6: Derivative of density w.r.t. temperature")
    {
        // This is here because we need r6_inverse to test that
        for(const auto& e: r6::detail::table_7)
        {
            INFO("D= ", e.D, ", T= ", e.T, ", P= ", e.P);
                const auto
            d_D_d_T = d_density_d_temperature_dt (e.D, e.T);
                const auto
            d_D_d_T_fd = central_finite_difference <0> (
                  r6_inverse::density_tp <double, double>
                , 1e-5
                , e.T
                , e.P
            );
            CHECK(d_D_d_T == Approx { d_D_d_T_fd }.scale (fabs (d_D_d_T_fd)).epsilon (1e-4));
        }
    }
} // TEST_CASE("r6_inverse.hpp")
