#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/isto/iapws/r6_inverse.hpp"
    using namespace isto::iapws::r6_inverse;
#include "../include/isto/iapws/detail/data_for_the_tests.hpp"
    using namespace isto::iapws::r6::r6_95_2016::detail;

TEST_CASE("r6_inverse.hpp")
{
        using namespace isto::iapws;
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
} // TEST_CASE("r6_inverse.hpp")
