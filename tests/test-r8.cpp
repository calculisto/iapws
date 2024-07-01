#include <doctest/doctest.h>
    using doctest::Approx;
#include "test.hpp"
#include "../include/calculisto/iapws/r8.hpp"
#include "../include/calculisto/iapws/r6_inverse.hpp"
    using namespace calculisto::iapws::r8;
    using namespace calculisto::iapws::r6_inverse;

TEST_CASE("r8.hpp")
{
SUBCASE ("details")
{
        using calculisto::iapws::r8::molar_mass_of_water;
    CHECK (density_pt (   0.101325e6, 240.) == Approx { 54.33701e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK (density_pt (   0.101325e6, 300.) == Approx { 55.31735e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK (density_pt (  10e6       , 300.) == Approx { 55.56148e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK (density_pt (1000e6       , 300.) == Approx { 68.69265e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK (density_pt (  10e6       , 650.) == Approx {  2.24692e3 * molar_mass_of_water }.scale (1e1).epsilon (1e-7));
    CHECK (density_pt ( 100e6       , 650.) == Approx { 40.31090e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK (density_pt ( 500e6       , 650.) == Approx { 52.58636e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
    CHECK (density_pt (  10e6       , 870.) == Approx {  1.45275e3 * molar_mass_of_water }.scale (1e1).epsilon (1e-5));
    CHECK (density_pt ( 100e6       , 870.) == Approx { 20.98927e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-6));
    CHECK (density_pt ( 500e6       , 870.) == Approx { 45.01376e3 * molar_mass_of_water }.scale (1e2).epsilon (1e-7));
}
SUBCASE ("main API")
{
    CHECK (relative_permittivity_dt (density_pt (   0.101325e6, 240.), 240.) == Approx { 104.34982 }.scale (1e2).epsilon (1e-7));
    CHECK (relative_permittivity_dt (density_pt (   0.101325e6, 300.), 300.) == Approx {  77.74735 }.scale (1e1).epsilon (1e-7));
    CHECK (relative_permittivity_dt (density_pt (  10e6       , 300.), 300.) == Approx {  78.11269 }.scale (1e1).epsilon (1e-7));
    CHECK (relative_permittivity_dt (density_pt (1000e6       , 300.), 300.) == Approx { 103.69632 }.scale (1e2).epsilon (1e-7));
    CHECK (relative_permittivity_dt (density_pt (  10e6       , 650.), 650.) == Approx {   1.26715 }.scale (1e0).epsilon (1e-5));
    CHECK (relative_permittivity_dt (density_pt ( 100e6       , 650.), 650.) == Approx {  17.71733 }.scale (1e1).epsilon (1e-6));
    CHECK (relative_permittivity_dt (density_pt ( 500e6       , 650.), 650.) == Approx {  26.62132 }.scale (1e1).epsilon (1e-7));
    CHECK (relative_permittivity_dt (density_pt (  10e6       , 870.), 870.) == Approx {   1.12721 }.scale (1e0).epsilon (1e-6));
    CHECK (relative_permittivity_dt (density_pt ( 100e6       , 870.), 870.) == Approx {   4.98281 }.scale (1e0).epsilon (1e-6));
    CHECK (relative_permittivity_dt (density_pt ( 500e6       , 870.), 870.) == Approx {  15.09746 }.scale (1e1).epsilon (1e-7));
} // SUBCASE ("main API")
SUBCASE ("Derivatives")
{
        struct
    data_t
    {
            double
        T;
            double
        p;
            double
        rho;
            double
        epsilon;
            double
        dedp;
            double
        dedt;
            double
        dedpp;
            double
        dedtt;
            double
        dedpt;
    };
    // https://doi.org/10.1063/1.555997
        constexpr auto
    data = std::array <data_t, 41>
    {{
          { 270     , 0.101325e6 , 55.4827e3 , 89.1821 , 0.0426805e-6  , -0.409375  , -0.56745e4 , 0.22655e-2 , -0.28474e3 }
        , { 300     , 0.101325e6 , 55.3174e3 , 77.7474 , 0.0371860e-6  , -0.355908  , -0.57134e4 , 0.15732e-2 , -0.11007e3 }
        , { 300     , 10.0e6     , 55.5615e3 , 78.1127 , 0.0366343e-6  , -0.357011  , -0.54389e4 , 0.16084e-2 , -0.11278e3 }
        , { 300     , 100.0e6    , 57.5729e3 , 81.2159 , 0.0325748e-6  , -0.367852  , -0.37492e4 , 0.19041e-2 , -0.12492e3 }
        , { 300     , 1000.0e6   , 68.6927e3 , 103.696 , 0.0212872e-6  , -0.481815  , -0.50341e5 , 0.35760e-2 , -0.12120e3 }
        , { 350     , 0.101325e6 , 54.0502e3 , 61.7889 , 0.0348273e-6  , -0.284834  , -0.77769e4 , 0.12805e-2 , -0.24820e5 }
        , { 350     , 10.0e6     , 54.2922e3 , 62.1299 , 0.0340834e-6  , -0.284884  , -0.72639e4 , 0.12933e-2 , -0.75556e5 }
        , { 350     , 100.0e6    , 56.2620e3 , 64.9510 , 0.0290383e-6  , -0.286956  , -0.43440e4 , 0.13909e-2 , -0.33972e4 }
        , { 350     , 1000.0e6   , 67.2951e3 , 83.6084 , 0.0168932e-6  , -0.334865  , -0.44415e5 , 0.21576e-2 , -0.57491e4 }
        , { 400     , 10.0e6     , 52.3122e3 , 49.3850 , 0.0350839e-6  , -0.227249  , -0.10647e3 , 0.10124e-2 ,  0.46103e4 }
        , { 400     , 100.0e6    , 54.4995e3 , 52.2009 , 0.0282478e-6  , -0.225663  , -0.54619e4 , 0.10752e-2 , -0.15232e5 }
        , { 400     , 1000.0e6   , 65.9422e3 , 69.1249 , 0.0147872e-6  , -0.251380  , -0.42334e5 , 0.13147e-2 , -0.32072e4 }
        , { 450     , 10.0e6     , 49.7447e3 , 39.1716 , 0.0389336e-6  , -0.183607  , -0.17795e3 , 0.73430e-3 ,  0.11397e3 }
        , { 450     , 100.0e6    , 52.3729e3 , 42.1495 , 0.0287410e-6  , -0.178549  , -0.73129e4 , 0.81662e-3 ,  0.21191e4 }
        , { 450     , 1000.0e6   , 64.5983e3 , 58.0200 , 0.0134011e-6  , -0.195865  , -0.40222e5 , 0.94387e-3 , -0.24640e4 }
        , { 500     , 10.0e6     , 46.5175e3 , 30.7941 , 0.0476630e-6  , -0.153818  , -0.36298e3 , 0.45362e-3 ,  0.25663e3 }
        , { 500     , 100.0e6    , 49.9140e3 , 34.1490 , 0.0304254e-6  , -0.143247  , -0.10411e3 , 0.60342e-3 ,  0.47211e4 }
        , { 500     , 1000.0e6   , 63.2530e3 , 49.3017 , 0.0122648e-6  , -0.154787  , -0.39098e5 , 0.71269e-3 , -0.21035e4 }
        , { 550     , 10.0e6     , 42.2875e3 , 23.5308 , 0.0695410e-6  , -0.139993  , -0.10916e2 , 0.48765e-4 ,  0.73366e3 }
        , { 550     , 100.0e6    , 47.1146e3 , 27.6672 , 0.0335966e-6  , -0.117403  , -0.15734e3 , 0.43823e-3 ,  0.81280e4 }
        , { 550     , 1000.0e6   , 61.9090e3 , 42.3775 , 0.0112878e-6  , -0.123594  , -0.39018e5 , 0.54304e-3 , -0.18089e4 }
        , { 600     , 100.0e6    , 43.9346e3 , 22.2903 , 0.0387450e-6  , -0.0986724 , -0.25140e3 , 0.31810e-3 ,  0.12680e3 }
        , { 600     , 1000.0e6   , 60.5715e3 , 36.8195 , 0.0104519e-6  , -0.0997915 , -0.39700e5 , 0.41510e-3 , -0.15389e4 }
        , { 650     , 100.0e6    , 40.3109e3 , 17.7173 , 0.0464881e-6  , -0.0848987 , -0.41946e3 , 0.23996e-3 ,  0.18454e3 }
        , { 650     , 1000.0e6   , 59.2461e3 , 32.3058 , 0.00974377e-6 , -0.0815521 , -0.40842e5 , 0.31914e-3 , -0.12995e4 }
        , { 700     , 100.0e6    , 36.1790e3 , 13.7544 , 0.0571346e-6  , -0.0738557 , -0.69943e3 , 0.21178e-3 ,  0.23622e3 }
        , { 700     , 1000.0e6   , 57.9372e3 , 28.5951 , 0.00914626e-6 , -0.0674716 , -0.42191e5 , 0.24759e-3 , -0.10969e4 }
        , { 750     , 100.0e6    , 31.5575e3 , 10.3360 , 0.0686834e-6  , -0.0625236 , -0.10220e2 , 0.25447e-3 ,  0.19846e3 }
        , { 750     , 1000.0e6   , 56.6489e3 , 25.5072 , 0.00864069e-6 , -0.0564897 , -0.43557e5 , 0.19428e-3 , -0.93141e5 }
        , { 800     , 100.0e6    , 26.7676e3 , 7.56225 , 0.0734450e-6  , -0.0478270 , -0.98076e3 , 0.32280e-3 , -0.30196e4 }
        , { 800     , 1000.0e6   , 55.3842e3 , 22.9077 , 0.00820919e-6 , -0.0478210 , -0.44801e5 , 0.15437e-3 , -0.79981e5 }
        , { 273.150 , 0.101325e6 , 55.4998e3 , 87.9035 , 0.0418297e-6  , -0.402570  , -0.55072e4 , 0.20674e-2 , -0.25617e3 }
        , { 273.150 , 100.0e6    , 58.0218e3 , 91.8380 , 0.0371765e-6  , -0.426331  , -0.40152e4 , 0.26183e-2 , -0.21982e3 }
        , { 298.144 , 0.101325e6 , 55.3447e3 , 78.4106 , 0.0373966e-6  , -0.358840  , -0.56606e4 , 0.15870e-2 , -0.11689e3 }
        , { 298.144 , 50.0e6     , 56.5321e3 , 80.2114 , 0.0348736e-6  , -0.364953  , -0.45377e4 , 0.17668e-2 , -0.12688e3 }
        , { 298.144 , 200.0e6    , 59.5001e3 , 85.0175 , 0.0296950e-6  , -0.384517  , -0.25672e4 , 0.22195e-2 , -0.13070e3 }
        , { 373.124 , 0.101325e6 , 53.1975e3 , 55.5333 , 0.0350994e-6  , -0.256700  , -0.92518e4 , 0.11525e-2 ,  0.25266e4 }
        , { 373.124 , 100.0e6    , 55.4962e3 , 58.6727 , 0.0284747e-6  , -0.256655  , -0.47909e4 , 0.12341e-2 , -0.16042e4 }
        , { 473.110 , 100.0e6    , 51.2774e3 , 38.2317 , 0.0293588e-6  , -0.160903  , -0.85449e4 , 0.71214e-3 ,  0.32458e4 }
        , { 673.102 , 100.0e6    , 38.4676e3 , 15.8180 , 0.0510736e-6  , -0.0796130 , -0.53399e3 , 0.21952e-3 ,  0.21205e3 }
        , { 773.071 , 100.0e6    , 29.3314e3 , 8.96472 , 0.0723386e-6  , -0.0562017 , -0.10733e2 , 0.29368e-3 ,  0.11028e3 }
    }};
    for (auto [T, p, rho, epsilon, dedp, dedt, dedpp, dedtt, dedpt]: data)
    {
            const auto
        d = density_pt (p, T);
        INFO("T= ", T, ", p= ", p, ", rho= ", rho * molar_mass_of_water);
        CHECK (density_pt (p, T) == Approx { rho * molar_mass_of_water }.scale (rho * molar_mass_of_water).epsilon (1e-6));
        CHECK (relative_permittivity_dt (d, T) == Approx { epsilon }.scale (epsilon).epsilon (1e-5));
        CHECK (d_relative_permittivity_d_p_dt (density_pt (p, T), T) == Approx { dedp }.scale (fabs (dedp)).epsilon (1e-5));
        CHECK (d_relative_permittivity_d_t_dt (density_pt (p, T), T) == Approx { dedt }.scale (fabs (dedt)).epsilon (1e-5));
        // The values in the reference are wrong, see below.
        // CHECK (d_relative_permittivity_dtt_dt (density_pt (p, T), T) == Approx { dedtt }.scale (fabs (dedtt)).epsilon (1e-5));
    }
    SUBCASE ("Our derivatives w.r.t. temperature and density are good")
    {
            auto
        dtr = [](auto delta, auto f, double r, double t, auto... args)
        {
            return (
                  f (r + delta, t + delta, args...)
                - f (r - delta, t + delta, args...)
                - f (r + delta, t - delta, args...)
                + f (r - delta, t - delta, args...)
            ) / 4 / delta / delta;
        };

        for (auto [T, p, rho, epsilon, dedp, dedt, dedpp, dedtt, dedpt]: data)
        {
            INFO("T= ", T, ", p= ", p, ", rho= ", rho * molar_mass_of_water);
                const auto
            expected = dtr (1e-1, e <double, double>, rho, T);
                const auto
            value = dedtr (rho, T);
            CHECK (value == Approx { expected }.scale (fabs (expected)).epsilon (1e-4));
        }
    }
} // SUBCASE ("Derivatives")
} // TEST_CASE("r10.hpp")
