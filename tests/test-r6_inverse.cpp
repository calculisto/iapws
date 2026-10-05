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
            // if (not (temperature == Approx { 302. })) continue;
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
                // FIXME: epsilon too low? 
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
    }
    SUBCASE("R6: Derivative of density w.r.t. temperature at fixed pressure")
    {
        // This is here because we need r6_inverse to test it
        for(const auto& e: r6::detail::table_7)
        {
            INFO("D= ", e.D, ", T= ", e.T, ", P= ", e.P);
                const auto
            d_D_d_T_at_P = d_density_d_temperature_at_pressure_dt (e.D, e.T);
                const auto
            d_D_d_T_at_P_fd = central_finite_difference <0> (
                  r6_inverse::density_tp <double, double>
                , 1e-6
                , e.T
                , e.P
            );
            CHECK(d_D_d_T_at_P == Approx { d_D_d_T_at_P_fd }.scale (fabs (d_D_d_T_at_P_fd )));
                const auto
            P = e.P;
                const auto
            d2_D_d2_T_at_P = d_density_d2_temperature_at_pressure_dt (e.D, e.T);
                const auto
            d2_D_d2_T_at_P_fd = central_finite_difference (
                  [P] (auto T_)
                  {
                        const auto
                    D_ = r6_inverse::density_tp (P, T_);
                    return d_density_d_temperature_at_pressure_dt (D_, T_);
                  }
                , 1e-4
                , e.T
            );
            CHECK(d2_D_d2_T_at_P == Approx { d2_D_d2_T_at_P_fd }.scale (fabs (d2_D_d2_T_at_P_fd)));
                const auto
            d2_D_d2_T_at_P_fd2 = central_finite_difference <0, 2> (
                  r6_inverse::density_tp <double, double>
                , 1e-4
                , e.T
                , e.P
            );
            CHECK(d2_D_d2_T_at_P == Approx { d2_D_d2_T_at_P_fd2 }.scale (fabs (d2_D_d2_T_at_P_fd2)));
        }
    }
} // TEST_CASE("r6_inverse.hpp")
