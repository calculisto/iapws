#pragma once
#include "detail/common.hpp"
#include "detail/template_pow.hpp"
    using calculisto::template_pow::pow;

// TODO long double coeffs ?

    namespace
calculisto::iapws::r7
{
    inline namespace
r7_97_2012
{
    constexpr auto
massic_gas_constant = 0.461526e3;

    constexpr auto
critical_temperature = 647.096;

    constexpr auto
critical_pressure = 22.064e6;

    constexpr auto
critical_density = 322.0;

// §5  Equations for Region 1.
    namespace
r1
{
    constexpr array_t
I //{{{
{
      0.0  // 1
    , 0.0  // 2
    , 0.0  // 3
    , 0.0  // 4
    , 0.0  // 5
    , 0.0  // 6
    , 0.0  // 7
    , 0.0  // 8
    , 1.0  // 9
    , 1.0  // 10
    , 1.0  // 11
    , 1.0  // 12
    , 1.0  // 13
    , 1.0  // 14
    , 2.0  // 15
    , 2.0  // 16
    , 2.0  // 17
    , 2.0  // 18
    , 2.0  // 19
    , 3.0  // 20
    , 3.0  // 21
    , 3.0  // 22
    , 4.0  // 23
    , 4.0  // 24
    , 4.0  // 25
    , 5.0  // 26
    , 8.0  // 27
    , 8.0  // 28
    , 21.0 // 29
    , 23.0 // 30
    , 29.0 // 31
    , 30.0 // 32
    , 31.0 // 33
    , 32.0 // 34
};//}}}

    constexpr array_t
J //{{{
{
      -2.0   // 1
    , -1.0   // 2
    ,  0.0   // 3
    ,  1.0   // 4
    ,  2.0   // 5
    ,  3.0   // 6
    ,  4.0   // 7
    ,  5.0   // 8
    , -9.0   // 9
    , -7.0   // 10
    , -1.0   // 11
    ,  0.0   // 12
    ,  1.0   // 13
    ,  3.0   // 14
    , -3.0   // 15
    ,  0.0   // 16
    ,  1.0   // 17
    ,  3.0   // 18
    ,  17.0  // 19
    , -4.0   // 20
    ,  0.0   // 21
    ,  6.0   // 22
    , -5.0   // 23
    , -2.0   // 24
    ,  10.0  // 25
    , -8.0   // 26
    , -11.0  // 27
    , -6.0   // 28
    , -29.0  // 29
    , -31.0  // 30
    , -38.0  // 31
    , -39.0  // 32
    , -40.0  // 33
    , -41.0  // 34
}; //}}}

    constexpr array_t
n //{{{
{
       0.14632971213167      // 1
    , -0.84548187169114      // 2
    , -0.37563603672040e1    // 3
    ,  0.33855169168385e1    // 4
    , -0.95791963387872      // 5
    ,  0.15772038513228      // 6
    , -0.16616417199501e-1  // 7
    ,  0.81214629983568e-3   // 8
    ,  0.28319080123804e-3   // 9
    , -0.60706301565874e-3   // 10
    , -0.18990068218419e-1   // 11
    , -0.32529748770505e-1   // 12
    , -0.21841717175414e-1   // 13
    , -0.52838357969930e-4   // 14
    , -0.47184321073267e-3   // 15
    , -0.30001780793026e-3   // 16
    ,  0.47661393906987e-4   // 17
    , -0.44141845330846e-5   // 18
    , -0.72694996297594e-15  // 19
    , -0.31679644845054e-4   // 20
    , -0.28270797985312e-5   // 21
    , -0.85205128120103e-9   // 22
    , -0.22425281908000e-5   // 23
    , -0.65171222895601e-6   // 24
    , -0.14341729937924e-12  // 25
    , -0.40516996860117e-6   // 26
    , -0.12734301741641e-8   // 27
    , -0.17424871230634e-9   // 28
    , -0.68762131295531e-18  // 29
    ,  0.14478307828521e-19  // 30
    ,  0.26335781662795e-22  // 31
    , -0.11947622640071e-22  // 32
    ,  0.18228094581404e-23  // 33
    , -0.93537087292458e-25  // 34
}; //}}}

// low-level dimensionless {{{
    constexpr auto
gamma (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n
        * pow (7.1 - pi, I)
        * pow (tau - 1.222, J)
    );
}
    constexpr auto
gamma_p (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
        - n
        * I * pow (7.1 - pi, I - 1.)
        * pow (tau - 1.222, J)
    );
}
    constexpr auto
gamma_t (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n
        * pow (7.1 - pi, I)
        * J * pow (tau - 1.222, J - 1.)
    );
}

    constexpr auto
gamma_pp (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n * I
        * (I - 1.) * pow (7.1 - pi, I - 2.)
        * pow (tau - 1.222, J)
    );
}
    constexpr auto
gamma_tt (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n
        * pow (7.1 - pi, I)
        * J * (J - 1.) * pow (tau - 1.222, J - 2.)
    );
}
    constexpr auto
gamma_pt (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
        - n * I
        * pow (7.1 - pi, I - 1.)
        * J * pow (tau - 1.222, J - 1.)
    );
}
// }}}

    constexpr auto
massic_volume_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = (1386. * (1)) / temperature;
    return pi * gamma_p (pi, tau) * massic_gas_constant * temperature / pressure;
}
    constexpr auto
massic_enthalpy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return tau * gamma_t (pi, tau) * massic_gas_constant * temperature;
}
    constexpr auto
massic_internal_energy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return (tau * gamma_t (pi, tau) - pi * gamma_p (pi, tau)) 
        * massic_gas_constant * temperature;
}
    constexpr auto
massic_entropy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return (tau * gamma_t (pi, tau) - gamma (pi, tau)) * massic_gas_constant;
}
    constexpr auto
massic_isobaric_heat_capacity_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return - tau * tau * gamma_tt (pi, tau) * massic_gas_constant;
}
    constexpr auto
massic_isochoric_heat_capacity_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return (- tau * tau * gamma_tt (pi, tau) + pow <2> (gamma_p (pi, tau) 
        - tau * gamma_pt (pi, tau)) / gamma_pp (pi, tau)) * massic_gas_constant;
}
    constexpr auto
speed_of_sound_pt (auto const& pressure, auto const& temperature)
{
        using std::pow;
        using std::sqrt;
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return sqrt (
          pow (gamma_p (pi, tau), 2)
        / (
              pow (gamma_p (pi, tau) - tau * gamma_pt (pi, tau), 2)
            / (tau * tau * gamma_tt (pi, tau))
            - gamma_pp (pi, tau)
          )
        * massic_gas_constant * temperature
    );
}
    constexpr auto
density_pt (auto const& pressure, auto const& temperature)
{
    return 1 / (massic_volume_pt (pressure, temperature));
}
    constexpr auto
isobaric_cubic_expansion_coefficient_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return (1. - tau * gamma_pt (pi, tau) / gamma_p (pi, tau)) / temperature;
}
    constexpr auto
isothermal_compressibility_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return (- pi * gamma_pp (pi, tau) / gamma_p (pi, tau)) / pressure;
}
    constexpr auto
relative_pressure_coefficient_pt (auto const& pressure, auto const& temperature)
{
    return isobaric_cubic_expansion_coefficient_pt (pressure, temperature) 
        / pressure / isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
isothermal_stress_coefficient_pt (auto const& pressure, auto const& temperature)
{
    return 1. / pressure / massic_volume_pt (pressure, temperature) 
        / isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
massic_gibbs_free_energy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = (1386. * (1)) / temperature;
    return gamma (pi, tau) * massic_gas_constant * temperature;
}
    constexpr auto
massic_helmholtz_free_energy_pt (auto const& pressure, auto const& temperature)
{
    return massic_internal_energy_pt (pressure, temperature) 
        - temperature * massic_entropy_pt(pressure, temperature);
}
    constexpr auto
temperature_ph (auto const& pressure, auto const& enthalpy)
{
    constexpr auto
I = array_t //{{{
{
      0. // 1
    , 0. // 2
    , 0. // 3
    , 0. // 4
    , 0. // 5
    , 0. // 6
    , 1. // 7
    , 1. // 8
    , 1. // 9
    , 1. // 10
    , 1. // 11
    , 1. // 12
    , 1. // 13
    , 2. // 14
    , 2. // 15
    , 3. // 16
    , 3. // 17
    , 4. // 18
    , 5. // 19
    , 6. // 20
}; //}}}
    constexpr auto
J = array_t //{{{
{
      0.  // 1
    , 1.  // 2
    , 2.  // 3
    , 6.  // 4
    , 22. // 5
    , 32. // 6
    , 0.  // 7
    , 1.  // 8
    , 2.  // 9
    , 3.  // 10
    , 4.  // 11
    , 10. // 12
    , 32. // 13
    , 10. // 14
    , 32. // 15
    , 10. // 16
    , 32. // 17
    , 32. // 18
    , 32. // 19
    , 32. // 20
}; //}}}
    constexpr auto
n = array_t  //{{{
{
      -0.23872489924521e3   // 1
    ,  0.40421188637945e3   // 2
    ,  0.11349746881718e3   // 3
    , -0.58457616048039e1   // 4
    , -0.15285482413140e-3  // 5
    , -0.10866707695377e-5  // 6
    , -0.13391744872602e2   // 7
    ,  0.43211039183559e2   // 8
    , -0.54010067170506e2   // 9
    ,  0.30535892203916e2   // 10
    , -0.65964749423638e1   // 11
    ,  0.93965400878363e-2  // 12
    ,  0.11573647505340e-6  // 13
    , -0.25858641282073e-4  // 14
    , -0.40644363084799e-8  // 15
    ,  0.66456186191635e-7  // 16
    ,  0.80670734103027e-10 // 17
    , -0.93477771213947e-12 // 18
    ,  0.58265442020601e-14 // 19
    , -0.15020185953503e-16 // 20
}; //}}}
        auto const
    pi = pressure / 1e6;
        auto const
    eta = enthalpy / (2500e3);
    return sum (n * pow (pi, I) * pow (eta + 1., J)) * (1);
}
    constexpr auto
temperature_ps (auto const& pressure, auto const& entropy)
{
    constexpr auto
I = array_t //{{{
{
      0. // 1
    , 0. // 2
    , 0. // 3
    , 0. // 4
    , 0. // 5
    , 0. // 6
    , 1. // 7
    , 1. // 8
    , 1. // 9
    , 1. // 10
    , 1. // 11
    , 1. // 12
    , 2. // 13
    , 2. // 14
    , 2. // 15
    , 2. // 16
    , 2. // 17
    , 3. // 18
    , 3. // 19
    , 4. // 20
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0. // 1
    ,  1. // 2
    ,  2. // 3
    ,  3. // 4
    , 11. // 5
    , 31. // 6
    ,  0. // 7
    ,  1. // 8
    ,  2. // 9
    ,  3. // 10
    , 12. // 11
    , 31. // 12
    ,  0. // 13
    ,  1. // 14
    ,  2. // 15
    ,  9. // 16
    , 31. // 17
    , 10. // 18
    , 32. // 19
    , 32. // 20
};  //}}}
    constexpr auto
n = array_t //{{{
{
       0.17478268058307e3    // 1
    ,  0.34806930892873e2    // 2
    ,  0.65292584978455e1    // 3
    ,  0.33039981775489      // 4
    , -0.19281382923196e-6   // 5
    , -0.24909197244573e-22  // 6
    , -0.26107636489332      // 7
    ,  0.22592965981586      // 8
    , -0.64256463395226e-1   // 9
    ,  0.78876289270526e-2   // 10
    ,  0.35672110607366e-9   // 11
    ,  0.17332496994895e-23  // 12
    ,  0.56608900654837e-3   // 13
    , -0.32635483139717e-3   // 14
    ,  0.44778286690632e-4   // 15
    , -0.51322156908507e-9   // 16
    , -0.42522657042207e-25  // 17
    ,  0.26400441360689e-12  // 18
    ,  0.78124600459723e-28  // 19
    , -0.30732199903668e-30  // 20
}; //}}}
        auto const
    pi = pressure / 1e6;
        auto const
    sigma = entropy / (1e3);
    return sum (n * pow (pi, I) * pow (sigma + 2., J)) * (1);
}
    constexpr auto
pressure_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
        constexpr auto
    I = array_t //{{{
    {
          0.
        , 0.
        , 0.
        , 0.
        , 0.
        , 0.
        , 0.
        , 0.
        , 1.
        , 1.
        , 1.
        , 1.
        , 2.
        , 2.
        , 2.
        , 3.
        , 4.
        , 4.
        , 5.
    }; //}}}
        constexpr auto
    J = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        ,  4.
        ,  5.
        ,  6.
        ,  8.
        , 14.
        ,  0.
        ,  1.
        ,  4.
        ,  6.
        ,  0.
        ,  1.
        , 10.
        ,  4.
        ,  1.
        ,  4.
        ,  0.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
          -0.691997014660582
        , -0.183612548787560e2
        , -0.928332409297335e1
        ,  0.659639569909906e2
        , -0.162060388912024e2
        ,  0.450620017338667e3
        ,  0.854680678224170e3
        ,  0.607523214001162e4
        ,  0.326487682621856e2
        , -0.269408844582931e2
        , -0.319947848334300e3
        , -0.928354307043320e3
        ,  0.303634537455249e2
        , -0.650540422444146e2
        , -0.430991316516130e4
        , -0.747512324096068e3
        ,  0.730000345529245e3
        ,  0.114284032569021e4
        , -0.436407041874559e3
    }; //}}}
        using std::pow;
        auto const
    eta = massic_enthalpy / (3400e3);
        auto const
    sigma = massic_entropy / (7.6e3);
    return sum (n * pow (eta + 0.05, I) * pow (sigma + 0.05, J)) * (100.);
}
    constexpr auto
temperature_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    return temperature_ph (pressure_hs (massic_enthalpy, massic_entropy), massic_enthalpy);
}
} // namespace r1
// §5  Equations for Region 2.
    namespace
r2
{
    constexpr array_t
J_0 // {{{
{
       0.0 // 1
    ,  1.0 // 2
    , -5.0 // 3
    , -4.0 // 4
    , -3.0 // 5
    , -2.0 // 6
    , -1.0 // 7
    ,  2.0 // 8
    ,  3.0 // 9
}; //}}}
    constexpr array_t
n_0 //{{{
{
      -0.96927686500217e1   // 1
    ,  0.10086655968018e2   // 2
    , -0.56087911283020e-2  // 3
    ,  0.71452738081455e-1  // 4
    , -0.40710498223928     // 5
    ,  0.14240819171444e1   // 6
    , -0.43839511319450e1   // 7
    , -0.28408632460772     // 8
    ,  0.21268463753307e-1  // 9
}; //}}}
    constexpr array_t
I_r //{{{
{
       1.0 // 1
    ,  1.0 // 2
    ,  1.0 // 3
    ,  1.0 // 4
    ,  1.0 // 5
    ,  2.0 // 6
    ,  2.0 // 7
    ,  2.0 // 8
    ,  2.0 // 9
    ,  2.0 // 10
    ,  3.0 // 11
    ,  3.0 // 12
    ,  3.0 // 13
    ,  3.0 // 14
    ,  3.0 // 15
    ,  4.0 // 16
    ,  4.0 // 17
    ,  4.0 // 18
    ,  5.0 // 19
    ,  6.0 // 20
    ,  6.0 // 21
    ,  6.0 // 22
    ,  7.0 // 23
    ,  7.0 // 24
    ,  7.0 // 25
    ,  8.0 // 26
    ,  8.0 // 27
    ,  9.0 // 28
    , 10.0 // 29
    , 10.0 // 30
    , 10.0 // 31
    , 16.0 // 32
    , 16.0 // 33
    , 18.0 // 34
    , 20.0 // 35
    , 20.0 // 36
    , 20.0 // 37
    , 21.0 // 38
    , 22.0 // 39
    , 23.0 // 40
    , 24.0 // 41
    , 24.0 // 42
    , 24.0 // 43
}; //}}}
    constexpr array_t
J_r //{{{
{
        0.0  // 1
     ,  1.0  // 2
     ,  2.0  // 3
     ,  3.0  // 4
     ,  6.0  // 5
     ,  1.0  // 6
     ,  2.0  // 7
     ,  4.0  // 8
     ,  7.0  // 9
     , 36.0  // 10
     ,  0.0  // 11
     ,  1.0  // 12
     ,  3.0  // 13
     ,  6.0  // 14
     , 35.0  // 15
     ,  1.0  // 16
     ,  2.0  // 17
     ,  3.0  // 18
     ,  7.0  // 19
     ,  3.0  // 20
     , 16.0  // 21
     , 35.0  // 22
     ,  0.0  // 23
     , 11.0  // 24
     , 25.0  // 25
     ,  8.0  // 26
     , 36.0  // 27
     , 13.0  // 28
     ,  4.0  // 29
     , 10.0  // 30
     , 14.0  // 31
     , 29.0  // 32
     , 50.0  // 33
     , 57.0  // 34
     , 20.0  // 35
     , 35.0  // 36
     , 48.0  // 37
     , 21.0  // 38
     , 53.0  // 39
     , 39.0  // 40
     , 26.0  // 41
     , 40.0  // 42
     , 58.0  // 43
}; //}}}
    constexpr array_t
n_r //{{{
{
      -0.17731742473213e-2     // 1
    , -0.17834862292358e-1     // 2
    , -0.45996013696365e-1     // 3
    , -0.57581259083432e-1     // 4
    , -0.50325278727930e-1     // 5
    , -0.33032641670203e-4     // 6
    , -0.18948987516315e-3     // 7
    , -0.39392777243355e-2     // 8
    , -0.43797295650573e-1     // 9
    , -0.26674547914087e-4     // 10
    ,  0.20481737692309e-7     // 11
    ,  0.43870667284435e-6     // 12
    , -0.32277677238570e-4     // 13
    , -0.15033924542148e-2     // 14
    , -0.40668253562649e-1     // 15
    , -0.78847309559367e-9     // 16
    ,  0.12790717852285e-7     // 17
    ,  0.48225372718507e-6     // 18
    ,  0.22922076337661e-5     // 19
    , -0.16714766451061e-10    // 20
    , -0.21171472321355e-2     // 21
    , -0.23895741934104e2      // 22
    , -0.59059564324270e-17    // 23
    , -0.12621808899101e-5     // 24
    , -0.38946842435739e-1     // 25
    ,  0.11256211360459e-10    // 26
    , -0.82311340897998e1      // 27
    ,  0.19809712802088e-7     // 28
    ,  0.10406965210174e-18    // 29
    , -0.10234747095929e-12    // 30
    , -0.10018179379511e-8     // 31
    , -0.80882908646985e-10    // 32
    ,  0.10693031879409        // 33
    , -0.33662250574171        // 34
    ,  0.89185845355421e-24    // 35
    ,  0.30629316876232e-12    // 36
    , -0.42002467698208e-5     // 37
    , -0.59056029685639e-25    // 38
    ,  0.37826947613457e-5     // 39
    , -0.12768608934681e-14    // 40
    ,  0.73087610595061e-28    // 41
    ,  0.55414715350778e-16    // 42
    , -0.94369707241210e-6     // 43
}; //}}}


// low-level dimensionless {{{
    constexpr auto
gamma_0 (auto const& pi, auto const& tau)
{
        using std::pow;
        using std::log;
    return log (pi) + sum (
          n_0
        * pow (tau, J_0)
    );
}
    constexpr auto
gamma_0_p (auto const& pi, auto const& /*tau*/)
{
    return 1 / pi;
}
    constexpr auto
gamma_0_pp (auto const& pi, auto const& /*tau*/)
{
    return - 1 / pi / pi;
}
    constexpr auto
gamma_0_t (auto const& /*pi*/, auto const& tau)
{
        using std::pow;
    return sum (
          n_0
        * J_0
        * pow (tau, J_0 - 1.)
    );
}
    constexpr auto
gamma_0_tt (auto const& /*pi*/, auto const& tau)
{
        using std::pow;
    return sum (
          n_0
        * J_0
        * (J_0 - 1.)
        * pow (tau, J_0 - 2.)
    );
}
    constexpr auto
gamma_r (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * pow (tau - 0.5, J_r)
    );
}
    constexpr auto
gamma_r_p (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * I_r
        * pow (pi, I_r - 1)
        * pow (tau - 0.5, J_r)
    );
}
    constexpr auto
gamma_r_pp (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * I_r
        * (I_r - 1.)
        * pow (pi, I_r - 2)
        * pow (tau - 0.5, J_r)
    );
}
    constexpr auto
gamma_r_t (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * J_r
        * pow (tau - 0.5, J_r - 1.)
    );
}
    constexpr auto
gamma_r_tt (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * J_r
        * (J_r - 1.)
        * pow (tau - 0.5, J_r - 2.)
    );
}
    constexpr auto
gamma_r_pt (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * I_r
        * pow (pi, I_r - 1.)
        * J_r
        * pow (tau - 0.5, J_r - 1.)
    );
}

    constexpr auto
gamma (auto const& pi, auto const& tau)
{
    return gamma_0 (pi, tau) + gamma_r (pi, tau);
}
    constexpr auto
gamma_p (auto const& pi, auto const& tau)
{
    return gamma_0_p (pi, tau) + gamma_r_p (pi, tau);
}
    constexpr auto
gamma_t (auto const& pi, auto const& tau)
{
    return gamma_0_t (pi, tau) + gamma_r_t (pi, tau);
}
    constexpr auto
gamma_pp (auto const& pi, auto const& tau)
{
    return gamma_0_pp (pi, tau) + gamma_r_pp (pi, tau);
}
    constexpr auto
gamma_tt (auto const& pi, auto const& tau)
{
    return gamma_0_tt (pi, tau) + gamma_r_tt (pi, tau);
}
    constexpr auto
gamma_pt (auto const& pi, auto const& tau)
{
    return /*gamma_0_pt (pi, tau) +*/ gamma_r_pt (pi, tau);
}
// }}}
    constexpr auto
massic_volume_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return pi * gamma_p (pi, tau) * massic_gas_constant * temperature / pressure;
}
    constexpr auto
massic_enthalpy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return tau * gamma_t (pi, tau) * massic_gas_constant * temperature;
}
    constexpr auto
massic_internal_energy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return (tau * gamma_t (pi, tau) - pi * gamma_p (pi, tau)) 
        * massic_gas_constant * temperature;
}
    constexpr auto
massic_entropy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return (tau * gamma_t (pi, tau) - gamma (pi, tau)) * massic_gas_constant;
}
    constexpr auto
massic_isobaric_heat_capacity_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return - tau * tau * gamma_tt (pi, tau) * massic_gas_constant;
}
    constexpr auto
massic_isochoric_heat_capacity_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return (- tau * tau * gamma_tt (pi, tau) + pow <2> (gamma_p (pi, tau) 
        - tau * gamma_pt (pi, tau)) / gamma_pp (pi, tau)) * massic_gas_constant;
}
    constexpr auto
speed_of_sound_pt (auto const& pressure, auto const& temperature)
{
        using std::pow;
        using std::sqrt;
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return sqrt (
          pow (1. + pi * gamma_r_p (pi, tau), 2)
        / (
              (1. - pi * pi * gamma_r_pp (pi, tau))
            + pow (1. + pi * gamma_r_p (pi, tau) - tau * pi * gamma_r_pt (pi, tau), 2)
            / (tau * tau * gamma_tt (pi, tau))
          )
        * massic_gas_constant * temperature
    );
}
    constexpr auto
density_pt (auto const& pressure, auto const& temperature)
{
    return 1 / massic_volume_pt (pressure, temperature);
}
    constexpr auto
isobaric_cubic_expansion_coefficient_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = (540.) / temperature;
    return (1. - tau * gamma_pt (pi, tau) / gamma_p (pi, tau)) / temperature;
}
    constexpr auto
isothermal_compressibility_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = (540.) / temperature;
    return (- pi * gamma_pp (pi, tau) / gamma_p (pi, tau)) / pressure;
}
    constexpr auto
relative_pressure_coefficient_pt (auto const& pressure, auto const& temperature)
{
    return isobaric_cubic_expansion_coefficient_pt (pressure, temperature) 
        / pressure / isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
isothermal_stress_coefficient_pt (auto const& pressure, auto const& temperature)
{
    return 1. / pressure / massic_volume_pt (pressure, temperature) 
        / isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
massic_gibbs_free_energy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = (540.) / temperature;
    return gamma (pi, tau) * massic_gas_constant * temperature;
}
    constexpr auto
massic_helmholtz_free_energy_pt (auto const& pressure, auto const& temperature)
{
    return massic_internal_energy_pt (pressure, temperature) 
        - temperature * massic_entropy_pt(pressure, temperature);
}
    constexpr auto
b2bc_h (auto const& massic_enthalpy)
{
        auto const
    eta = massic_enthalpy / (1e3);
    return (0.90584278514723e3 - 0.67955786399241 * eta + 0.12809002730136e-3 * eta * eta) * 1e6;
}
    constexpr auto
b2bc_p (auto const& pressure)
{
        using std::sqrt;
        auto const
    pi = pressure / 1e6;
    return (0.26526571908428e4 + sqrt ((pi - 0.45257578905948e1 ) / 0.12809002730136e-3)) * (1e3);
}
    constexpr auto
subregion_ph (auto const& pressure, auto const& massic_enthalpy)
{
    if (pressure < 4. * 1e6) return 1;
    if (pressure < b2bc_h (massic_enthalpy)) return 2;
    return 3;
}
    constexpr auto
subregion_ps (auto const& pressure, auto const& massic_entropy)
{
    if (pressure < 4. * 1e6) return 1;
    if (massic_entropy >= 5.85 * (1e3)) return 2;
    return 3;
}
    namespace
t_ph
{
    namespace
a
{
    constexpr auto
I = array_t //{{{
{
      0.
    , 0.
    , 0.
    , 0.
    , 0.
    , 0.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 2.
    , 2.
    , 2.
    , 2.
    , 2.
    , 2.
    , 2.
    , 2.
    , 3.
    , 3.
    , 4.
    , 4.
    , 4.
    , 5.
    , 5.
    , 5.
    , 6.
    , 6.
    , 7.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    ,  2.
    ,  3.
    ,  7.
    , 20.
    ,  0.
    ,  1.
    ,  2.
    ,  3.
    ,  7.
    ,  9.
    , 11.
    , 18.
    , 44.
    ,  0.
    ,  2.
    ,  7.
    , 36.
    , 38.
    , 40.
    , 42.
    , 44.
    , 24.
    , 44.
    , 12.
    , 32.
    , 44.
    , 32.
    , 36.
    , 42.
    , 34.
    , 44.
    , 28.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.10898952318288e4
    ,  0.84951654495535e3
    , -0.10781748091826e3
    ,  0.33153654801263e2
    , -0.74232016790248e1
    ,  0.11765048724356e2
    ,  0.18445749355790e1
    , -0.41792700549624e1
    ,  0.62478196935812e1
    , -0.17344563108114e2
    , -0.20058176862096e3
    ,  0.27196065473796e3
    , -0.45511318285818e3
    ,  0.30919688604755e4
    ,  0.25226640357872e6
    , -0.61707422868339e-2
    , -0.31078046629583
    ,  0.11670873077107e2
    ,  0.12812798404046e9
    , -0.98554909623276e9
    ,  0.28224546973002e10
    , -0.35948971410703e10
    ,  0.17227349913197e10
    , -0.13551334240775e5
    ,  0.12848734664650e8
    ,  0.13865724283226e1
    ,  0.23598832556514e6
    , -0.13105236545054e8
    ,  0.73999835474766e4
    , -0.55196697030060e6
    ,  0.37154085996233e7
    ,  0.19127729239660e5
    , -0.41535164835634e6
    , -0.62459855192507e2
}; //}}}
} // namespace a
    namespace
b
{
    constexpr auto
I = array_t //{{{
{
      0.
    , 0.
    , 0.
    , 0.
    , 0.
    , 0.
    , 0.
    , 0.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 2.
    , 2.
    , 2.
    , 2.
    , 3.
    , 3.
    , 3.
    , 3.
    , 4.
    , 4.
    , 4.
    , 4.
    , 4.
    , 4.
    , 5.
    , 5.
    , 5.
    , 6.
    , 7.
    , 7.
    , 9.
    , 9.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    ,  2.
    , 12.
    , 18.
    , 24.
    , 28.
    , 40.
    ,  0.
    ,  2.
    ,  6.
    , 12.
    , 18.
    , 24.
    , 28.
    , 40.
    ,  2.
    ,  8.
    , 18.
    , 40.
    ,  1.
    ,  2.
    , 12.
    , 24.
    ,  2.
    , 12.
    , 18.
    , 24.
    , 28.
    , 40.
    , 18.
    , 24.
    , 40.
    , 28.
    ,  2.
    , 28.
    ,  1.
    , 40.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.14895041079516e4
    ,  0.74307798314034e3
    , -0.97708318797837e2
    ,  0.24742464705674e1
    , -0.63281320016026
    ,  0.11385952129658e1
    , -0.47811863648625
    ,  0.85208123431544e-2
    ,  0.93747147377932
    ,  0.33593118604916e1
    ,  0.33809355601454e1
    ,  0.16844539671904
    ,  0.73875745236695
    , -0.47128737436186
    ,  0.15020273139707
    , -0.21764114219750e-2
    , -0.21810755324761e-1
    , -0.10829784403677
    , -0.46333324635812e-1
    ,  0.71280351959551e-4
    ,  0.11032831789999e-3
    ,  0.18955248387902e-3
    ,  0.30891541160537e-2
    ,  0.13555504554949e-2
    ,  0.28640237477456e-6
    , -0.10779857357512e-4
    , -0.76462712454814e-4
    ,  0.14052392818316e-4
    , -0.31083814331434e-4
    , -0.10302738212103e-5
    ,  0.28217281635040e-6
    ,  0.12704902271945e-5
    ,  0.73803353468292e-7
    , -0.11030139238909e-7
    , -0.81456365207833e-13
    , -0.25180545682962e-10
    , -0.17565233969407e-17
    ,  0.86934156344163e-14
}; //}}}
} // namespace b
    namespace
c
{
    constexpr auto
I = array_t //{{{
{
      -7.
    , -7.
    , -6.
    , -6.
    , -5.
    , -5.
    , -2.
    , -2.
    , -1.
    , -1.
    ,  0.
    ,  0.
    ,  1.
    ,  1.
    ,  2.
    ,  6.
    ,  6.
    ,  6.
    ,  6.
    ,  6.
    ,  6.
    ,  6.
    ,  6.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  4.
    ,  0.
    ,  2.
    ,  0.
    ,  2.
    ,  0.
    ,  1.
    ,  0.
    ,  2.
    ,  0.
    ,  1.
    ,  4.
    ,  8.
    ,  4.
    ,  0.
    ,  1.
    ,  4.
    , 10.
    , 12.
    , 16.
    , 20.
    , 22.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.32368398555242e13
    ,  0.73263350902181e13
    ,  0.35825089945447e12
    , -0.58340131851590e12
    , -0.10783068217470e11
    ,  0.20825544563171e11
    ,  0.61074783564516e6
    ,  0.85977722535580e6
    , -0.25745723604170e5
    ,  0.31081088422714e5
    ,  0.12082315865936e4
    ,  0.48219755109255e3
    ,  0.37966001272486e1
    , -0.10842984880077e2
    , -0.45364172676660e-1
    ,  0.14559115658698e-12
    ,  0.11261597407230e-11
    , -0.17804982240686e-10
    ,  0.12324579690832e-6
    , -0.11606921130984e-5
    ,  0.27846367088554e-4
    , -0.59270038474176e-3
    ,  0.12918582991878e-2
}; //}}}
} // namespace cc
} // namespace h
    namespace
t_ps
{
    namespace
a
{
    constexpr auto
I = array_t //{{{
{
      -1.5
    , -1.5
    , -1.5
    , -1.5
    , -1.5
    , -1.5
    , -1.25
    , -1.25
    , -1.25
    , -1.0
    , -1.0
    , -1.0
    , -1.0
    , -1.0
    , -1.0
    , -0.75
    , -0.75
    , -0.5
    , -0.5
    , -0.5
    , -0.5
    , -0.25
    , -0.25
    , -0.25
    , -0.25
    ,  0.25
    ,  0.25
    ,  0.25
    ,  0.25
    ,  0.5
    ,  0.5
    ,  0.5
    ,  0.5
    ,  0.5
    ,  0.5
    ,  0.5
    ,  0.75
    ,  0.75
    ,  0.75
    ,  0.75
    ,  1.0
    ,  1.0
    ,  1.25
    ,  1.25
    ,  1.5
    ,  1.5
}; //}}}
    constexpr auto
J = array_t //{{{
{
      -24.
    , -23.
    , -19.
    , -13.
    , -11.
    , -10.
    , -19.
    , -15.
    ,  -6.
    , -26.
    , -21.
    , -17.
    , -16.
    ,  -9.
    ,  -8.
    , -15.
    , -14.
    , -26.
    , -13.
    ,  -9.
    ,  -7.
    , -27.
    , -25.
    , -11.
    ,  -6.
    ,   1.
    ,   4.
    ,   8.
    ,  11.
    ,   0.
    ,   1.
    ,   5.
    ,   6.
    ,  10.
    ,  14.
    ,  16.
    ,   0.
    ,   4.
    ,   9.
    ,  17.
    ,   7.
    ,  18.
    ,   3.
    ,  15.
    ,   5.
    ,  18.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.39235983861984e6
    ,  0.51526573827270e6
    ,  0.40482443161048e5
    , -0.32193790923902e3
    ,  0.96961424218694e2
    , -0.22867846371773e2
    , -0.44942914124357e6
    , -0.50118336020166e4
    ,  0.35684463560015
    ,  0.44235335848190e5
    , -0.13673388811708e5
    ,  0.42163260207864e6
    ,  0.22516925837475e5
    ,  0.47442144865646e3
    , -0.14931130797647e3
    , -0.19781126320452e6
    , -0.23554399470760e5
    , -0.19070616302076e5
    ,  0.55375669883164e5
    ,  0.38293691437363e4
    , -0.60391860580567e3
    ,  0.19363102620331e4
    ,  0.42660643698610e4
    , -0.59780638872718e4
    , -0.70401463926862e3
    ,  0.33836784107553e3
    ,  0.20862786635187e2
    ,  0.33834172656196e-1
    , -0.43124428414893e-4
    ,  0.16653791356412e3
    , -0.13986292055898e3
    , -0.78849547999872
    ,  0.72132411753872e-1
    , -0.59754839398283e-2
    , -0.12141358953904e-4
    ,  0.23227096733871e-6
    , -0.10538463566194e2
    ,  0.20718925496502e1
    , -0.72193155260427e-1
    ,  0.20749887081120e-6
    , -0.18340657911379e-1
    ,  0.29036272348696e-6
    ,  0.21037527893619
    ,  0.25681239729999e-3
    , -0.12799002933781e-1
    , -0.82198102652018e-5
}; //}}}
} // namespace a
    namespace
b
{
    constexpr auto
I = array_t //{{{
{
      -6.
    , -6.
    , -5.
    , -5.
    , -4.
    , -4.
    , -4.
    , -3.
    , -3.
    , -3.
    , -3.
    , -2.
    , -2.
    , -2.
    , -2.
    , -1.
    , -1.
    , -1.
    , -1.
    , -1.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  2.
    ,  2.
    ,  2.
    ,  3.
    ,  3.
    ,  3.
    ,  4.
    ,  4.
    ,  5.
    ,  5.
    ,  5.
}; //}}}
    constexpr auto
J = array_t //{{{
{
        0.
    ,  11.
    ,   0.
    ,  11.
    ,   0.
    ,   1.
    ,  11.
    ,   0.
    ,   1.
    ,  11.
    ,  12.
    ,   0.
    ,   1.
    ,   6.
    ,  10.
    ,   0.
    ,   1.
    ,   5.
    ,   8.
    ,   9.
    ,   0.
    ,   1.
    ,   2.
    ,   4.
    ,   5.
    ,   6.
    ,   9.
    ,   0.
    ,   1.
    ,   2.
    ,   3.
    ,   7.
    ,   8.
    ,   0.
    ,   1.
    ,   5.
    ,   0.
    ,   1.
    ,   3.
    ,   0.
    ,   1.
    ,   0.
    ,   1.
    ,   2.
}; //}}}
    constexpr auto
n = array_t //{{{
{
        0.31687665083497e6
    ,   0.20864175881858e2
    ,  -0.39859399803599e6
    ,  -0.21816058518877e2
    ,   0.22369785194242e6
    ,  -0.27841703445817e4
    ,   0.99207436071480e1
    ,  -0.75197512299157e5
    ,   0.29708605951158e4
    ,  -0.34406878548526e1
    ,   0.3881556424911
    ,   0.17511295085750e5
    ,  -0.14237112854449e4
    ,   0.10943803364167e1
    ,   0.89971619308495
    ,  -0.33759740098958e4
    ,   0.47162885818355e3
    ,  -0.19188241993679e1
    ,   0.41078580492196
    ,  -0.33465378172097
    ,   0.13870034777505e4
    ,  -0.40663326195838e3
    ,   0.41727347159610e2
    ,   0.21932549434532e1
    ,  -0.10320050009077e1
    ,   0.35882943516703
    ,   0.52511453726066e-2
    ,   0.12838916450705e2
    ,  -0.28642437219381e1
    ,   0.56912683664855
    ,  -0.99962954584931e-1
    ,  -0.32632037778459e-2
    ,   0.23320922576723e-3
    ,  -0.15334809857450
    ,   0.29072288239902e-1
    ,   0.37534702741167e-3
    ,   0.17296691702411e-2
    ,  -0.38556050844504e-3
    ,  -0.35017712292608e-4
    ,  -0.14566393631492e-4
    ,   0.56420857267269e-5
    ,   0.41286150074605e-7
    ,  -0.20684671118824e-7
    ,   0.16409393674725e-8
}; //}}}
} // namespace b
    namespace
c
{
    constexpr auto
I = array_t //{{{
{
      -2.
    , -2.
    , -1.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  2.
    ,  2.
    ,  2.
    ,  3.
    ,  3.
    ,  3.
    ,  4.
    ,  4.
    ,  4.
    ,  5.
    ,  5.
    ,  5.
    ,  6.
    ,  6.
    ,  7.
    ,  7.
    ,  7.
    ,  7.
    ,  7.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      0.
    , 1.
    , 0.
    , 0.
    , 1.
    , 2.
    , 3.
    , 0.
    , 1.
    , 3.
    , 4.
    , 0.
    , 1.
    , 2.
    , 0.
    , 1.
    , 5.
    , 0.
    , 1.
    , 4.
    , 0.
    , 1.
    , 2.
    , 0.
    , 1.
    , 0.
    , 1.
    , 3.
    , 4.
    , 5.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.90968501005365e3
    ,  0.24045667088420e4
    , -0.59162326387130e3
    ,  0.54145404128074e3
    , -0.27098308411192e3
    ,  0.97976525097926e3
    , -0.46966772959435e3
    ,  0.14399274604723e2
    , -0.19104204230429e2
    ,  0.53299167111971e1
    , -0.21252975375934e2
    , -0.31147334413760
    ,  0.60334840894623
    , -0.42764839702509e-1
    ,  0.58185597255259e-2
    , -0.14597008284753e-1
    ,  0.56631175631027e-2
    , -0.76155864584577e-4
    ,  0.22440342919332e-3
    , -0.12561095013413e-4
    ,  0.63323132660934e-6
    , -0.20541989675375e-5
    ,  0.36405370390082e-7
    , -0.29759897789215e-8
    ,  0.10136618529763e-7
    ,  0.59925719692351e-11
    , -0.20677870105164e-10
    , -0.20874278181886e-10
    ,  0.10162166825089e-9
    , -0.16429828281347e-9
}; //}}}
} // namespace cc
} // namespace s

    constexpr auto
temperature_ph (auto const& pressure, auto const& massic_enthalpy)
{
        auto const
    pi = pressure / 1e6;
        auto const
    eta = massic_enthalpy / (2000e3);
        using namespace t_ph;
    switch (subregion_ph (pressure, massic_enthalpy))
    {
        case 1: return sum (a::n * pow (pi,       a::I) * pow (eta - 2.1, a::J));
        case 2: return sum (b::n * pow (pi - 2.,  b::I) * pow (eta - 2.6, b::J));
        case 3: return sum (c::n * pow (pi + 25., c::I) * pow (eta - 1.8, c::J));
    }
    throw internal_error_e {};
}
    constexpr auto
temperature_ps (auto const& pressure, auto const& massic_entropy)
{
        auto const
    pi = pressure / 1e6;
        using namespace t_ps;
    switch (subregion_ps (pressure, massic_entropy))
    {
        case 1:
        {
                auto const
            sigma = massic_entropy / (2e3);
            return sum (a::n * pow (pi, a::I) * pow (sigma - 2., a::J)) ;
        }
        case 2:
        {
                auto const
            sigma = massic_entropy / (0.7853e3);
            return sum (b::n * pow (pi, b::I) * pow (10. - sigma, b::J));
        }
        case 3:
        {
                auto const
            sigma = massic_entropy / (2.9251e3);
            return sum (c::n * pow (pi, c::I) * pow (2. - sigma, c::J)) ;
        }
    }
    throw internal_error_e {};
}
    constexpr auto
b2ab_s (auto const& massic_entropy)
{
        auto const
    sigma = massic_entropy / (1e3);
    return (
          -0.349898083432139e4
        +  0.257560716905876e4 * sigma
        + -0.421073558227969e3 * sigma * sigma
        +  0.276349063799944e2 * sigma * sigma * sigma
    ) * (1e3);
}
    constexpr auto
subregion_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    if (massic_enthalpy < b2ab_s (massic_entropy)) return 1;
    if (massic_entropy >= 5.85 * (1e3)) return 2;
    return 3;
}

    namespace
p_hs
{
    namespace
a
{
    constexpr auto
I = array_t ///{{{
{
      0.
    , 0.
    , 0.
    , 0.
    , 0.
    , 0.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 1.
    , 2.
    , 2.
    , 2.
    , 3.
    , 3.
    , 3.
    , 3.
    , 3.
    , 4.
    , 5.
    , 5.
    , 6.
    , 7.
}; //}}}

    constexpr auto
J = array_t ///{{{
{
       1.
    ,  3.
    ,  6.
    , 16.
    , 20.
    , 22.
    ,  0.
    ,  1.
    ,  2.
    ,  3.
    ,  5.
    ,  6.
    , 10.
    , 16.
    , 20.
    , 22.
    ,  3.
    , 16.
    , 20.
    ,  0.
    ,  2.
    ,  3.
    ,  6.
    , 16.
    , 16.
    ,  3.
    , 16.
    ,  3.
    ,  1.
}; //}}}

    constexpr auto
n = array_t ///{{{
{
      -0.182575361923032e-1
    , -0.125229548799536
    ,  0.592290437320145
    ,  0.604769706185122e1
    ,  0.238624965444474e3
    , -0.298639090222922e3
    ,  0.512250813040750e-1
    , -0.437266515606486
    ,  0.413336902999504
    , -0.516468254574773e1
    , -0.557014838445711e1
    ,  0.128555037824478e2
    ,  0.114144108953290e2
    , -0.119504225652714e3
    , -0.284777985961560e4
    ,  0.431757846408006e4
    ,  0.112894040802650e1
    ,  0.197409186206319e4
    ,  0.151612444706087e4
    ,  0.141324451421235e-1
    ,  0.585501282219601
    , -0.297258075863012e1
    ,  0.594567314847319e1
    , -0.623656565798905e4
    ,  0.965986235133332e4
    ,  0.681500934948134e1
    , -0.633207286824489e4
    , -0.558919224465760e1
    ,  0.400645798472063e-1
}; //}}}
} // namespace a
    namespace
b
{
    constexpr auto
I = array_t ///{{{
{
       0.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  2.
    ,  2.
    ,  2.
    ,  3.
    ,  3.
    ,  3.
    ,  3.
    ,  4.
    ,  4.
    ,  5.
    ,  5.
    ,  6.
    ,  6.
    ,  6.
    ,  7.
    ,  7.
    ,  8.
    ,  8.
    ,  8.
    ,  8.
    , 12.
    , 14.
}; //}}}

    constexpr auto
J = array_t ///{{{
{
       0.
    ,  1.
    ,  2.
    ,  4.
    ,  8.
    ,  0.
    ,  1.
    ,  2.
    ,  3.
    ,  5.
    , 12.
    ,  1.
    ,  6.
    , 18.
    ,  0.
    ,  1.
    ,  7.
    , 12.
    ,  1.
    , 16.
    ,  1.
    , 12.
    ,  1.
    ,  8.
    , 18.
    ,  1.
    , 16.
    ,  1.
    ,  3.
    , 14.
    , 18.
    , 10.
    , 16.
}; //}}}

    constexpr auto
n = array_t ///{{{
{
       0.801496989929495e-1
    , -0.543862807146111
    ,  0.337455597421283
    ,  0.890555451157450e1
    ,  0.313840736431485e3
    ,  0.797367065977789
    , -0.121616973556240e1
    ,  0.872803386937477e1
    , -0.169769781757602e2
    , -0.186552827328416e3
    ,  0.951159274344237e5
    , -0.189168510120494e2
    , -0.433407037194840e4
    ,  0.543212633012715e9
    ,  0.144793408386013
    ,  0.128024559637516e3
    , -0.672309534071268e5
    ,  0.336972380095287e8
    , -0.586634196762720e3
    , -0.221403224769889e11
    ,  0.171606668708389e4
    , -0.570817595806302e9
    , -0.312109693178482e4
    , -0.207841384633010e7
    ,  0.305605946157786e13
    ,  0.322157004314333e4
    ,  0.326810259797295e12
    , -0.144104158934487e4
    ,  0.410694867802691e3
    ,  0.109077066873024e12
    , -0.247964654258893e14
    ,  0.188801906865134e10
    , -0.123651009018773e15
}; //}}}
} // namespace b
    namespace
c
{
    constexpr auto
I = array_t ///{{{
{
       0.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  2.
    ,  2.
    ,  2.
    ,  2.
    ,  2.
    ,  3.
    ,  3.
    ,  3.
    ,  3.
    ,  3.
    ,  4.
    ,  5.
    ,  5.
    ,  5.
    ,  5.
    ,  6.
    ,  6.
    , 10.
    , 12.
    , 16.
}; //}}}

    constexpr auto
J = array_t ///{{{
{
       0.
    ,  1.
    ,  2.
    ,  3.
    ,  4.
    ,  8.
    ,  0.
    ,  2.
    ,  5.
    ,  8.
    , 14.
    ,  2.
    ,  3.
    ,  7.
    , 10.
    , 18.
    ,  0.
    ,  5.
    ,  8.
    , 16.
    , 18.
    , 18.
    ,  1.
    ,  4.
    ,  6.
    , 14.
    ,  8.
    , 18.
    ,  7.
    ,  7.
    , 10.
}; //}}}

    constexpr auto
n = array_t ///{{{
{
       0.112225607199012
    , -0.339005953606712e1
    , -0.320503911730094e2
    , -0.197597305104900e3
    , -0.407693861553446e3
    ,  0.132943775222331e5
    ,  0.170846839774007e1
    ,  0.373694198142245e2
    ,  0.358144365815434e4
    ,  0.423014446424664e6
    , -0.751071025760063e9
    ,  0.523446127607898e2
    , -0.228351290812417e3
    , -0.960652417056937e6
    , -0.807059292526074e8
    ,  0.162698017225669e13
    ,  0.772465073604171
    ,  0.463929973837746e5
    , -0.137317885134128e8
    ,  0.170470392630512e13
    , -0.251104628187308e14
    ,  0.317748830835520e14
    ,  0.538685623675312e2
    , -0.553089094625169e5
    , -0.102861522421405e7
    ,  0.204249418756234e13
    ,  0.273918446626977e9
    , -0.263963146312685e16
    , -0.107890854108088e10
    , -0.296492620980124e11
    , -0.111754907323424e16
}; //}}}
} // namespace c
} // namespace p_hs
    constexpr auto
pressure_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
        using namespace p_hs;
    switch (subregion_hs (massic_enthalpy, massic_entropy))
    {
        case 1:
        {
                auto const
            eta = massic_enthalpy / (4200e3);
                auto const
            sigma = massic_entropy / (12e3);
            return pow (sum (a::n * pow (eta - 0.5, a::I) * pow (sigma - 1.2, a::J)), 4) * (4e6);
        }
        case 2:
        {
                auto const
            eta = massic_enthalpy / (4100e3);
                auto const
            sigma = massic_entropy / (7.9e3);
            return pow (sum (b::n * pow (eta - 0.6, b::I) * pow (sigma - 1.01, b::J)), 4) * (100e6);
        }
        case 3:
        {
                auto const
            eta = massic_enthalpy / (3500e3);
                auto const
            sigma = massic_entropy / (5.9e3);
            return pow (sum (c::n * pow (eta - 0.7, c::I) * pow (sigma - 1.1, c::J)), 4) * (100e6);
        }
    }
    throw internal_error_e {};
}
    constexpr auto
temperature_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    return temperature_ph (pressure_hs (massic_enthalpy, massic_entropy), massic_enthalpy);
}
} //namespace r2

    namespace
r2_metastable_vapor
{
    constexpr array_t
n_0 //{{{
{
      -0.96937268393049e1   // 1
    ,  0.10087275970006e2   // 2
    , -0.56087911283020e-2  // 3
    ,  0.71452738081455e-1  // 4
    , -0.40710498223928     // 5
    ,  0.14240819171444e1   // 6
    , -0.43839511319450e1   // 7
    , -0.28408632460772     // 8
    ,  0.21268463753307e-1  // 9
}; //}}}
    constexpr array_t
I_r //{{{
{
      1.
    , 1.
    , 1.
    , 1.
    , 2.
    , 2.
    , 2.
    , 3.
    , 3.
    , 4.
    , 4.
    , 5.
    , 5.
}; //}}}
    constexpr array_t
J_r //{{{
{
       0.
    ,  2.
    ,  5.
    , 11.
    ,  1.
    ,  7.
    , 16.
    ,  4.
    , 16.
    ,  7.
    , 10.
    ,  9.
    , 10.
}; //}}}
    constexpr array_t
n_r //{{{
{
      -0.73362260186506e-2
    , -0.88223831943146e-1
    , -0.72334555213245e-1
    , -0.40813178534455e-2
    ,  0.20097803380207e-2
    , -0.53045921898642e-1
    , -0.76190409086970e-2
    , -0.63498037657313e-2
    , -0.86043093028588e-1
    ,  0.75321581522770e-2
    , -0.79238375446139e-2
    , -0.22888160778447e-3
    , -0.26456501482810e-2
}; //}}}
// low-level dimensionless {{{

    using r2::J_0;
    constexpr auto
gamma_0 (auto const& pi, auto const& tau)
{
        using std::pow;
        using std::log;
    return log (pi) + sum (
          n_0
        * pow (tau, J_0)
    );
}
    constexpr auto
gamma_0_p (auto const& pi, auto const& /*tau*/)
{
    return 1 / pi;
}
    constexpr auto
gamma_0_pp (auto const& pi, auto const& /*tau*/)
{
    return - 1 / pi / pi;
}
    constexpr auto
gamma_0_t (auto const& /*pi*/, auto const& tau)
{
        using std::pow;
    return sum (
          n_0
        * J_0
        * pow (tau, J_0 - 1.)
    );
}
    constexpr auto
gamma_0_tt (auto const& /*pi*/, auto const& tau)
{
        using std::pow;
    return sum (
          n_0
        * J_0
        * (J_0 - 1.)
        * pow (tau, J_0 - 2.)
    );
}
    constexpr auto
gamma_r (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * pow (tau - 0.5, J_r)
    );
}
    constexpr auto
gamma_r_p (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * I_r
        * pow (pi, I_r - 1)
        * pow (tau - 0.5, J_r)
    );
}
    constexpr auto
gamma_r_pp (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * I_r
        * (I_r - 1.)
        * pow (pi, I_r - 2)
        * pow (tau - 0.5, J_r)
    );
}
    constexpr auto
gamma_r_t (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * J_r
        * pow (tau - 0.5, J_r - 1.)
    );
}
    constexpr auto
gamma_r_tt (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * J_r
        * (J_r - 1.)
        * pow (tau - 0.5, J_r - 2.)
    );
}
    constexpr auto
gamma_r_pt (auto const& pi, auto const& tau)
{
        using std::pow;
    return sum (
          n_r
        * I_r
        * pow (pi, I_r - 1.)
        * J_r
        * pow (tau - 0.5, J_r - 1.)
    );
}

    constexpr auto
gamma (auto const& pi, auto const& tau)
{
    return gamma_0 (pi, tau) + gamma_r (pi, tau);
}
    constexpr auto
gamma_p (auto const& pi, auto const& tau)
{
    return gamma_0_p (pi, tau) + gamma_r_p (pi, tau);
}
    constexpr auto
gamma_t (auto const& pi, auto const& tau)
{
    return gamma_0_t (pi, tau) + gamma_r_t (pi, tau);
}
    constexpr auto
gamma_pp (auto const& pi, auto const& tau)
{
    return gamma_0_pp (pi, tau) + gamma_r_pp (pi, tau);
}
    constexpr auto
gamma_tt (auto const& pi, auto const& tau)
{
    return gamma_0_tt (pi, tau) + gamma_r_tt (pi, tau);
}
    constexpr auto
gamma_pt (auto const& pi, auto const& tau)
{
    return /*gamma_0_pt (pi, tau) +*/ gamma_r_pt (pi, tau);
}
// }}}
    constexpr auto
massic_volume_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return pi * gamma_p (pi, tau) * massic_gas_constant * temperature / pressure;
}
    constexpr auto
massic_enthalpy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return tau * gamma_t (pi, tau) * massic_gas_constant * temperature;
}
    constexpr auto
massic_internal_energy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return (tau * gamma_t (pi, tau) - pi * gamma_p (pi, tau)) * massic_gas_constant * temperature;
}
    constexpr auto
massic_entropy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return (tau * gamma_t (pi, tau) - gamma (pi, tau)) * massic_gas_constant;
}
    constexpr auto
massic_isobaric_heat_capacity_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return - tau * tau * gamma_tt (pi, tau) * massic_gas_constant;
}
    constexpr auto
massic_isochoric_heat_capacity_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return (- tau * tau * gamma_tt (pi, tau) + pow <2> (gamma_p (pi, tau) - tau * gamma_pt (pi, tau)) / gamma_pp (pi, tau)) * massic_gas_constant;
}
    constexpr auto
speed_of_sound_pt (auto const& pressure, auto const& temperature)
{
        using std::pow;
        using std::sqrt;
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return sqrt (
          pow (1. + pi * gamma_r_p (pi, tau), 2)
        / (
              (1. - pi * pi * gamma_r_pp (pi, tau))
            + pow (1. + pi * gamma_r_p (pi, tau) - tau * pi * gamma_r_pt (pi, tau), 2)
            / (tau * tau * gamma_tt (pi, tau))
          )
        * massic_gas_constant * temperature
    );
}
    constexpr auto
massic_gibbs_free_energy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 540 * (1) / temperature;
    return gamma (pi, tau) * massic_gas_constant * temperature;
}
} // namespace r2_metastable_vapor

// Saturated region
    namespace
pt_sat
{
    constexpr auto
n = array_t //{{{
{
       0.11670521452767e4
    , -0.72421316703206e6
    , -0.17073846940092e2
    ,  0.12020824702470e5
    , -0.32325550322333e7
    ,  0.14915108613530e2
    , -0.48232657361591e4
    ,  0.40511340542057e6
    , -0.23855557567849
    ,  0.65017534844798e3
}; //}}}
} // namespace detail
    auto
saturation_pressure_t (auto const& temperature)
{
        using std::pow;
        using std::sqrt;
        using namespace pt_sat;
        auto const
    t = temperature / (1);
        auto const
    theta = t + n[8] / (t - n[9]);
        auto const
    A = theta * theta + n[0] * theta + n[1];
        auto const
    B = n[2] * theta * theta + n[3] * theta + n[4];
        auto const
    C = n[5] * theta * theta + n[6] * theta + n[7];
    return pow (2. * C / (-B + sqrt (B * B - 4. * A * C)), 4) * 1e6;
}

    constexpr auto
saturation_temperature_p (auto const& pressure)
{
        using std::pow;
        using std::sqrt;
        using namespace pt_sat;
        auto const
    beta = pow (pressure / 1e6, 0.25);
        auto const
    E = beta * beta + n[2] * beta + n[5];
        auto const
    F = n[0] * beta * beta + n[3] * beta + n[6];
        auto const
    G = n[1] * beta * beta + n[4] * beta + n[7];
        auto const
    D = 2. * G / (-F - sqrt (F * F - 4. * E * G));
    return (n[9] + D - sqrt (pow (n[9] + D, 2) - 4. * (n[8] + n[9] * D))) / 2.;
}
    namespace
t_sat_hs
{
    constexpr auto
I = array_t //{{{
{
        0.
    ,   0.
    ,   0.
    ,   1.
    ,   1.
    ,   1.
    ,   1.
    ,   2.
    ,   2.
    ,   2.
    ,   3.
    ,   3.
    ,   3.
    ,   3.
    ,   4.
    ,   4.
    ,   5.
    ,   5.
    ,   5.
    ,   5.
    ,   6.
    ,   6.
    ,   6.
    ,   8.
    ,  10.
    ,  10.
    ,  12.
    ,  14.
    ,  14.
    ,  16.
    ,  16.
    ,  18.
    ,  18.
    ,  18.
    ,  20.
    ,  28.
}; //}}}

    constexpr auto
J = array_t //{{{
{
       0.
    ,  3.
    , 12.
    ,  0.
    ,  1.
    ,  2.
    ,  5.
    ,  0.
    ,  5.
    ,  8.
    ,  0.
    ,  2.
    ,  3.
    ,  4.
    ,  0.
    ,  1.
    ,  1.
    ,  2.
    ,  4.
    , 16.
    ,  6.
    ,  8.
    , 22.
    ,  1.
    , 20.
    , 36.
    , 24.
    ,  1.
    , 28.
    , 12.
    , 32.
    , 14.
    , 22.
    , 36.
    , 24.
    , 36.
}; //}}}

    constexpr auto
n = array_t //{{{
{
       0.179882673606601
    , -0.267507455199603
    ,  0.116276722612600e1
    ,  0.147545428713616
    , -0.512871635973248
    ,  0.421333567697984
    ,  0.563749522189870
    ,  0.429274443819153
    , -0.335704552142140e1
    ,  0.108890916499278e2
    , -0.248483390456012
    ,  0.304153221906390
    , -0.494819763939905
    ,  0.107551674933261e1
    ,  0.733888415457688e-1
    ,  0.140170545411085e-1
    , -0.106110975998808
    ,  0.168324361811875e-1
    ,  0.125028363714877e1
    ,  0.101316840309509e4
    , -0.151791558000712e1
    ,  0.524277865990866e2
    ,  0.230495545563912e5
    ,  0.249459806365456e-1
    ,  0.210796467412137e7
    ,  0.366836848613065e9
    , -0.144814105365163e9
    , -0.179276373003590e-2
    ,  0.489955602100459e10
    ,  0.471262212070518e3
    , -0.829294390198652e11
    , -0.171545662263191e4
    ,  0.355777682973575e7
    ,  0.586062760258436e12
    , -0.129887635078195e8
    ,  0.317247449371057e11
}; //}}}
} // namespace t_sat_hs
    constexpr auto
saturation_temperature_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
        using namespace t_sat_hs;
        auto const
    eta = massic_enthalpy / (2800e3);
        auto const
    sigma = massic_entropy / (9.2e3);
    return sum (n * pow (eta - 0.119, I) * pow (sigma - 1.07, J)) * (550);
}
    namespace
p_sat_h
{
    constexpr auto
I = array_t //{{{
{
       0.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  5.
    ,  7.
    ,  8.
    , 14.
    , 20.
    , 22.
    , 24.
    , 28.
    , 36.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    ,  3.
    ,  4.
    , 36.
    ,  3.
    ,  0.
    , 24.
    , 16.
    , 16.
    ,  3.
    , 18.
    ,  8.
    , 24.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.600073641753024
    , -0.936203654849857e1
    ,  0.246590798594147e2
    , -0.107014222858224e3
    , -0.915821315805768e14
    , -0.862332011700662e4
    , -0.235837344740032e2
    ,  0.252304969384128e18
    , -0.389718771997719e19
    , -0.333775713645296e23
    ,  0.356499469636328e11
    , -0.148547544720641e27
    ,  0.330611514838798e19
    ,  0.813641294467829e38
}; //}}}
} // namespace p_sat_h
    constexpr auto
saturation_pressure_h (auto const& massic_enthalpy)
{
        using namespace p_sat_h;
        auto const
    eta = massic_enthalpy / (2600e3);
    return sum (n * pow (eta - 1.02, I) * pow (eta - 0.608, J)) * (22e6);
}
    namespace
p_sat_s
{
    constexpr auto
I = array_t //{{{
{
       0.
    ,  1.
    ,  1.
    ,  4.
    , 12.
    , 12.
    , 16.
    , 24.
    , 28.
    , 32.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    , 32.
    ,  7.
    ,  4.
    , 14.
    , 36.
    , 10.
    ,  0.
    , 18.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.639767553612785
    , -0.129727445396014e2
    , -0.224595125848403e16
    ,  0.177466741801846e7
    ,  0.717079349571538e10
    , -0.378829107169011e18
    , -0.955586736431328e35
    ,  0.187269814676188e24
    ,  0.119254746466473e12
    ,  0.110649277244882e37
}; //}}}
} //namespace p_sat_s
    constexpr auto
saturation_pressure_s (auto const& massic_entropy)
{
        using namespace p_sat_s;
        auto const
    sigma = massic_entropy / (5.2e3);
    return sum (n * pow (sigma - 1.03, I) * pow (sigma - 0.699, J)) * (22e6);
}

// §5  Equations for Region 3.
    namespace
r3
{
    namespace
direct
{
    constexpr array_t
I
{ //{{{
       0.0 //  2 XXX WARNING starts at 2 XXX
    ,  0.0 //  3
    ,  0.0 //  4
    ,  0.0 //  5
    ,  0.0 //  6
    ,  0.0 //  7
    ,  0.0 //  8
    ,  1.0 //  9
    ,  1.0 // 10
    ,  1.0 // 11
    ,  1.0 // 12
    ,  2.0 // 13
    ,  2.0 // 14
    ,  2.0 // 15
    ,  2.0 // 16
    ,  2.0 // 17
    ,  2.0 // 18
    ,  3.0 // 19
    ,  3.0 // 20
    ,  3.0 // 21
    ,  3.0 // 22
    ,  3.0 // 23
    ,  4.0 // 24
    ,  4.0 // 25
    ,  4.0 // 26
    ,  4.0 // 27
    ,  5.0 // 28
    ,  5.0 // 29
    ,  5.0 // 30
    ,  6.0 // 31
    ,  6.0 // 32
    ,  6.0 // 33
    ,  7.0 // 34
    ,  8.0 // 35
    ,  9.0 // 36
    ,  9.0 // 37
    , 10.0 // 38
    , 10.0 // 39
    , 11.0 // 40
}; //}}}
    constexpr array_t
J //{{{
{
       0.0  //  2 XXX WARNING XXX
    ,  1.0  //  3
    ,  2.0  //  4
    ,  7.0  //  5
    , 10.0  //  6
    , 12.0  //  7
    , 23.0  //  8
    ,  2.0  //  9
    ,  6.0  // 10
    , 15.0  // 11
    , 17.0  // 12
    ,  0.0  // 13
    ,  2.0  // 14
    ,  6.0  // 15
    ,  7.0  // 16
    , 22.0  // 17
    , 26.0  // 18
    ,  0.0  // 19
    ,  2.0  // 20
    ,  4.0  // 21
    , 16.0  // 22
    , 26.0  // 23
    ,  0.0  // 24
    ,  2.0  // 25
    ,  4.0  // 26
    , 26.0  // 27
    ,  1.0  // 28
    ,  3.0  // 29
    , 26.0  // 30
    ,  0.0  // 31
    ,  2.0  // 32
    , 26.0  // 33
    ,  2.0  // 34
    , 26.0  // 35
    ,  2.0  // 36
    , 26.0  // 37
    ,  0.0  // 38
    ,  1.0  // 39
    , 26.0  // 40
}; //}}}

    constexpr double
n1 = 0.10658070028513e1; //  1 XXX WARNING XXX

    constexpr array_t
n //{{{
{
      -0.15732845290239e2   //  2 XXX WARNING XXX -0.15732845290239e2
    ,  0.20944396974307e2   //  3                  0.20944396974307e2
    , -0.76867707878716e1   //  4                 -0.76867707878716e1
    ,  0.26185947787954e1   //  5                  0.26185947787954e1
    , -0.28080781148620e1   //  6                 -0.28080781148620e1
    ,  0.12053369696517e1   //  7                  0.12053369696517e1
    , -0.84566812812502e-2  //  8                 -0.84566812812502e-2
    , -0.12654315477714e1   //  9                 -0.12654315477714e1
    , -0.11524407806681e1   // 10                 -0.11524407806681e1
    ,  0.88521043984318     // 11                  0.88521043984318
    , -0.64207765181607     // 12                 -0.64207765181607
    ,  0.38493460186671     // 13                  0.38493460186671
    , -0.85214708824206     // 14                 -0.85214708824206
    ,  0.48972281541877e1   // 15                  0.48972281541877e1
    , -0.30502617256965e1   // 16                 -0.30502617256965e1
    ,  0.39420536879154e-1  // 17                  0.39420536879154e-1
    ,  0.12558408424308     // 18                  0.12558408424308
    , -0.27999329698710     // 19                 -0.27999329698710
    ,  0.13899799569460e1   // 20                  0.13899799569460e1
    , -0.20189915023570e1   // 21                 -0.20189915023570e1
    , -0.82147637173963e-2  // 22                 -0.82147637173963e-2
    , -0.47596035734923     // 23                 -0.47596035734923
    ,  0.43984074473500e-1  // 24                  0.43984074473500e-1
    , -0.44476435428739     // 25                 -0.44476435428739
    ,  0.90572070719733     // 26                  0.90572070719733
    ,  0.70522450087967     // 27                  0.70522450087967
    ,  0.10770512626332     // 28                  0.10770512626332
    , -0.32913623258954     // 29                 -0.32913623258954
    , -0.50871062041158     // 30                 -0.50871062041158
    , -0.22175400873096e-1  // 31                 -0.22175400873096e-1
    ,  0.94260751665092e-1  // 32                  0.94260751665092e-1
    ,  0.16436278447961     // 33                  0.16436278447961
    , -0.13503372241348e-1  // 34                 -0.13503372241348e-1
    , -0.14834345352472e-1  // 35                 -0.14834345352472e-1
    ,  0.57922953628084e-3  // 36                  0.57922953628084e-3
    ,  0.32308904703711e-2  // 37                  0.32308904703711e-2
    ,  0.80964802996215e-4  // 38                  0.80964802996215e-4
    , -0.16557679795037e-3  // 39                 -0.16557679795037e-3
    , -0.44923899061815e-4  // 40                 -0.44923899061815e-4
}; //}}}
} // namespace direct
// low-level dimensionless {{{
    constexpr auto
phi (auto const& delta, auto const& tau)
{
        using namespace direct;
        using std::log;
        using std::pow;
    return
          n1
        * log (delta)
        + sum (
              n
            * pow (delta, I)
            * pow (tau, J)
          )
    ;
}
    constexpr auto
phi_d (auto const& delta, auto const& tau)
{
        using namespace direct;
        using std::pow;
    return
          n1
        / delta
        + sum (
              n
            * I
            * pow (delta, I - 1.)
            * pow (tau, J)
          )
    ;
}
    constexpr auto
phi_dd (auto const& delta, auto const& tau)
{
        using namespace direct;
        using std::pow;
    return
        - n1
        / delta
        / delta
        + sum (
              n
            * I
            * (I - 1.)
            * pow (delta, I - 2.)
            * pow (tau, J)
          )
    ;
}
    constexpr auto
phi_t (auto const& delta, auto const& tau)
{
        using namespace direct;
        using std::pow;
    return
          sum (
              n
            * pow (delta, I)
            * J
            * pow (tau, J - 1.)
          )
    ;
}
    constexpr auto
phi_tt (auto const& delta, auto const& tau)
{
        using namespace direct;
        using std::pow;
    return
          sum (
              n
            * pow (delta, I)
            * J
            * (J - 1.)
            * pow (tau, J - 2.)
          )
    ;
}
    constexpr auto
phi_dt (auto const& delta, auto const& tau)
{
        using namespace direct;
        using std::pow;
    return
          sum (
              n
            * I
            * pow (delta, I - 1.)
            * J
            * pow (tau, J - 1.)
          )
    ;
}
// }}}
    constexpr auto
pressure_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (delta * phi_d (delta, tau)) * density * massic_gas_constant * temperature;
}
    constexpr auto
massic_internal_energy_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (tau * phi_t (delta, tau)) * massic_gas_constant * temperature;
}
    constexpr auto
massic_entropy_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (tau * phi_t (delta, tau) - phi (delta, tau)) * massic_gas_constant;
}
    constexpr auto
massic_enthalpy_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (tau * phi_t (delta, tau) + delta * phi_d (delta, tau)) * massic_gas_constant * temperature;
}
    constexpr auto
massic_isochoric_heat_capacity_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (- tau * tau * phi_tt (delta, tau)) * massic_gas_constant;
}
    constexpr auto
massic_isobaric_heat_capacity_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (- tau * tau * phi_tt (delta, tau) + pow (delta * phi_d (delta, tau) - delta * tau * phi_dt (delta, tau), 2) / (2. * delta * phi_d (delta, tau) + delta * delta * phi_dd (delta, tau))) * massic_gas_constant;
}
    constexpr auto
speed_of_sound_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return sqrt ((
          2. * delta * phi_d (delta, tau)
        + delta * delta * phi_dd (delta, tau)
        - pow (delta * phi_d (delta, tau) - delta * tau * phi_dt (delta, tau), 2)
        / (tau * tau * phi_tt (delta, tau))
    ) * massic_gas_constant * temperature);
}
    constexpr auto
isobaric_cubic_expansion_coefficient_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return ((phi_d (delta, tau) - tau * phi_dt (delta, tau)) / (2. * phi_d (delta, tau) + delta * phi_dd (delta, tau))) / temperature;
}
    constexpr auto
isothermal_compressibility_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (1. / (2. * delta * phi_d (delta, tau) + delta * delta * phi_dd (delta, tau))) / density / massic_gas_constant / temperature;
}
    constexpr auto
relative_pressure_coefficient_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (1. - tau * phi_dt (delta, tau) / phi_d (delta, tau)) / temperature;
}
    constexpr auto
isothermal_stress_coefficient_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (2. + delta * phi_dd (delta, tau) / phi_d (delta, tau)) * density;
}
    constexpr auto
massic_helmholtz_free_energy_dt (auto const& density, auto const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return phi (delta, tau) * massic_gas_constant * temperature;
}
    constexpr auto
massic_gibbs_free_energy_dt (auto const& density, auto const& temperature)
{
    return massic_helmholtz_free_energy_dt (density, temperature) 
        - temperature * massic_entropy_dt(density, temperature);
}
    constexpr auto
b3ab_p (auto const& pressure)
{
        auto const
    pi = pressure / 1e6;
    return (
           0.201464004206875e4
        +  0.374696550136983e1  * pi
        + -0.219921901054187e-1 * pi * pi
        +  0.875131686009950e-4 * pi * pi * pi
    ) * ( 1e3);
}
    namespace
t_ph
{
    namespace
a
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    , -12.
    , -12.
    , -12.
    , -12.
    , -12.
    , -12.
    , -10.
    , -10.
    , -10.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -5.
    ,  -3.
    ,  -2.
    ,  -2.
    ,  -2.
    ,  -1.
    ,  -1.
    ,   0.
    ,   0.
    ,   1.
    ,   3.
    ,   3.
    ,   4.
    ,   4.
    ,  10.
    ,  12.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    ,  2.
    ,  6.
    , 14.
    , 16.
    , 20.
    , 22.
    ,  1.
    ,  5.
    , 12.
    ,  0.
    ,  2.
    ,  4.
    , 10.
    ,  2.
    ,  0.
    ,  1.
    ,  3.
    ,  4.
    ,  0.
    ,  2.
    ,  0.
    ,  1.
    ,  1.
    ,  0.
    ,  1.
    ,  0.
    ,  3.
    ,  4.
    ,  5.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.133645667811215e-6
    ,  0.455912656802978e-5
    , -0.146294640700979e-4
    ,  0.639341312970080e-2
    ,  0.372783927268847e3
    , -0.718654377460447e4
    ,  0.573494752103400e6
    , -0.267569329111439e7
    , -0.334066283302614e-4
    , -0.245479214069597e-1
    ,  0.478087847764996e2
    ,  0.764664131818904e-5
    ,  0.128350627676972e-2
    ,  0.171219081377331e-1
    , -0.851007304583213e1
    , -0.136513461629781e-1
    , -0.384460997596657e-5
    ,  0.337423807911655e-2
    , -0.551624873066791
    ,  0.729202277107470
    , -0.992522757376041e-2
    , -0.119308831407288
    ,  0.793929190615421
    ,  0.454270731799386
    ,  0.209998591259910
    , -0.642109823904738e-2
    , -0.235155868604540e-1
    ,  0.252233108341612e-2
    , -0.764885133368119e-2
    ,  0.136176427574291e-1
    , -0.133027883575669e-1
}; //}}}
} // namespace a
    namespace
b
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    , -10.
    , -10.
    , -10.
    , -10.
    , -10.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -4.
    ,  -4.
    ,  -3.
    ,  -2.
    ,  -2.
    ,  -1.
    ,  -1.
    ,  -1.
    ,  -1.
    ,  -1.
    ,  -1.
    ,   0.
    ,   0.
    ,   1.
    ,   3.
    ,   5.
    ,   6.
    ,   8.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    ,  0.
    ,  1.
    ,  5.
    , 10.
    , 12.
    ,  0.
    ,  1.
    ,  2.
    ,  4.
    , 10.
    ,  0.
    ,  1.
    ,  2.
    ,  0.
    ,  1.
    ,  5.
    ,  0.
    ,  4.
    ,  2.
    ,  4.
    ,  6.
    , 10.
    , 14.
    , 16.
    ,  0.
    ,  2.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.323254573644920e-4
    , -0.127575556587181e-3
    , -0.475851877356068e-3
    ,  0.156183014181602e-2
    ,  0.105724860113781
    , -0.858514221132534e2
    ,  0.724140095480911e3
    ,  0.296475810273257e-2
    , -0.592721983365988e-2
    , -0.126305422818666e-1
    , -0.115716196364853
    ,  0.849000969739595e2
    , -0.108602260086615e-1
    ,  0.154304475328851e-1
    ,  0.750455441524466e-1
    ,  0.252520973612982e-1
    , -0.602507901232996e-1
    , -0.307622221350501e1
    , -0.574011959864879e-1
    ,  0.503471360939849e1
    , -0.925081888584834
    ,  0.391733882917546e1
    , -0.773146007130190e2
    ,  0.949308762098587e4
    , -0.141043719679409e7
    ,  0.849166230819026e7
    ,  0.861095729446704
    ,  0.323346442811720
    ,  0.873281936020439
    , -0.436653048526683
    ,  0.286596714529479
    , -0.131778331276228
    ,  0.676682064330275e-2
}; //}}}
} // namespace b
} // namespace t_ph
    constexpr auto
temperature_ph (auto const& pressure, auto const& massic_enthalpy)
{
        using namespace t_ph;
        auto const
    pi = pressure / (100e6);
    if (massic_enthalpy < b3ab_p (pressure))
    {
            auto const
        eta = massic_enthalpy / (2300e3);
        return sum (a::n * pow (pi + 0.240, a::I) * pow (eta - 0.615, a::J)) * (760);
    }
        auto const
    eta = massic_enthalpy / (2800e3);
    return sum (b::n * pow (pi + 0.298, b::I) * pow (eta - 0.720, b::J)) * (860);
}
    namespace
v_ph
{
    namespace
a
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    , -12.
    , -12.
    , -10.
    , -10.
    , -10.
    ,  -8.
    ,  -8.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -4.
    ,  -4.
    ,  -3.
    ,  -2.
    ,  -2.
    ,  -1.
    ,  -1.
    ,  -1.
    ,  -1.
    ,   0.
    ,   0.
    ,   1.
    ,   1.
    ,   1.
    ,   2.
    ,   2.
    ,   3.
    ,   4.
    ,   5.
    ,   8.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       6.
    ,  8.
    , 12.
    , 18.
    ,  4.
    ,  7.
    , 10.
    ,  5.
    , 12.
    ,  3.
    ,  4.
    , 22.
    ,  2.
    ,  3.
    ,  7.
    ,  3.
    , 16.
    ,  0.
    ,  1.
    ,  2.
    ,  3.
    ,  0.
    ,  1.
    ,  0.
    ,  1.
    ,  2.
    ,  0.
    ,  2.
    ,  0.
    ,  2.
    ,  2.
    ,  2.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.529944062966028e-2
    , -0.170099690234461
    ,  0.111323814312927e2
    , -0.217898123145125e4
    , -0.506061827980875e-3
    ,  0.556495239685324
    , -0.943672726094016e1
    , -0.297856807561527
    ,  0.939353943717186e2
    ,  0.192944939465981e-1
    ,  0.421740664704763
    , -0.368914126282330e7
    , -0.737566847600639e-2
    , -0.354753242424366
    , -0.199768169338727e1
    ,  0.115456297059049e1
    ,  0.568366875815960e4
    ,  0.808169540124668e-2
    ,  0.172416341519307
    ,  0.104270175292927e1
    , -0.297691372792847
    ,  0.560394465163593
    ,  0.275234661176914
    , -0.148347894866012
    , -0.651142513478515e-1
    , -0.292468715386302e1
    ,  0.664876096952665e-1
    ,  0.352335014263844e1
    , -0.146340792313332e-1
    , -0.224503486668184e1
    ,  0.110533464706142e1
    , -0.408757344495612e-1
}; //}}}
} //namespace a
    namespace
b
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -4.
    ,  -4.
    ,  -4.
    ,  -3.
    ,  -3.
    ,  -2.
    ,  -2.
    ,  -1.
    ,  -1.
    ,  -1.
    ,  -1.
    ,   0.
    ,   1.
    ,   1.
    ,   2.
    ,   2.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    ,  0.
    ,  1.
    ,  3.
    ,  6.
    ,  7.
    ,  8.
    ,  0.
    ,  1.
    ,  2.
    ,  5.
    ,  6.
    , 10.
    ,  3.
    ,  6.
    , 10.
    ,  0.
    ,  2.
    ,  1.
    ,  2.
    ,  0.
    ,  1.
    ,  4.
    ,  5.
    ,  0.
    ,  0.
    ,  1.
    ,  2.
    ,  6.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.225196934336318e-8
    ,  0.140674363313486e-7
    ,  0.233784085280560e-5
    , -0.331833715229001e-4
    ,  0.107956778514318e-2
    , -0.271382067378863
    ,  0.107202262490333e1
    , -0.853821329075382
    , -0.215214194340526e-4
    ,  0.769656088222730e-3
    , -0.431136580433864e-2
    ,  0.453342167309331
    , -0.507749535873652
    , -0.100475154528389e3
    , -0.219201924648793
    , -0.321087965668917e1
    ,  0.607567815637771e3
    ,  0.557686450685932e-3
    ,  0.187499040029550
    ,  0.905368030448107e-2
    ,  0.285417173048685
    ,  0.329924030996098e-1
    ,  0.239897419685483
    ,  0.482754995951394e1
    , -0.118035753702231e2
    ,  0.169490044091791
    , -0.179967222507787e-1
    ,  0.371810116332674e-1
    , -0.536288335065096e-1
    ,  0.160697101092520e1
}; //}}}
} // namespace b
} //namespace v_ph
    constexpr auto
massic_volume_ph (auto const& pressure, auto const& massic_enthalpy)
{
        using namespace v_ph;
        auto const
    pi = pressure / (100e6);
    if (massic_enthalpy < b3ab_p (pressure))
    {
            auto const
        eta = massic_enthalpy / (2100e3);
        return sum (a::n * pow (pi + 0.128, a::I) * pow (eta - 0.727, a::J)) * (0.0028);
    }
        auto const
    eta = massic_enthalpy / (2800e3);
    return sum (b::n * pow (pi + 0.0661, b::I) * pow (eta - 0.720, b::J)) * (0.0088);
}
    namespace
t_ps
{
    namespace
a
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    , -10.
    , -10.
    , -10.
    , -10.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -5.
    ,  -5.
    ,  -5.
    ,  -4.
    ,  -4.
    ,  -4.
    ,  -2.
    ,  -2.
    ,  -1.
    ,  -1.
    ,   0.
    ,   0.
    ,   0.
    ,   1.
    ,   2.
    ,   2.
    ,   3.
    ,   8.
    ,   8.
    ,  10.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      28.
    , 32.
    ,  4.
    , 10.
    , 12.
    , 14.
    ,  5.
    ,  7.
    ,  8.
    , 28.
    ,  2.
    ,  6.
    , 32.
    ,  0.
    , 14.
    , 32.
    ,  6.
    , 10.
    , 36.
    ,  1.
    ,  4.
    ,  1.
    ,  6.
    ,  0.
    ,  1.
    ,  4.
    ,  0.
    ,  0.
    ,  3.
    ,  2.
    ,  0.
    ,  1.
    ,  2.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.150042008263875e10
    , -0.159397258480424e12
    ,  0.502181140217975e-3
    , -0.672057767855466e2
    ,  0.145058545404456e4
    , -0.823889534888890e4
    , -0.154852214233853
    ,  0.112305046746695e2
    , -0.297000213482822e2
    ,  0.438565132635495e11
    ,  0.137837838635464e-2
    , -0.297478527157462e1
    ,  0.971777947349413e13
    , -0.571527767052398e-4
    ,  0.288307949778420e5
    , -0.744428289262703e14
    ,  0.128017324848921e2
    , -0.368275545889071e3
    ,  0.664768904779177e16
    ,  0.449359251958880e-1
    , -0.422897836099655e1
    , -0.240614376434179
    , -0.474341365254924e1
    ,  0.724093999126110
    ,  0.923874349695897
    ,  0.399043655281015e1
    ,  0.384066651868009e-1
    , -0.359344365571848e-2
    , -0.735196448821653
    ,  0.188367048396131
    ,  0.141064266818704e-3
    , -0.257418501496337e-2
    ,  0.123220024851555e-2
}; //}}}
} // namespace a
    namespace
b
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    , -12.
    , -12.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -5.
    ,  -5.
    ,  -5.
    ,  -5.
    ,  -5.
    ,  -4.
    ,  -3.
    ,  -3.
    ,  -2.
    ,   0.
    ,   2.
    ,   3.
    ,   4.
    ,   5.
    ,   6.
    ,   8.
    ,  12.
    ,  14.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       1.
    ,  3.
    ,  4.
    ,  7.
    ,  0.
    ,  1.
    ,  3.
    ,  0.
    ,  2.
    ,  4.
    ,  0.
    ,  1.
    ,  2.
    ,  4.
    ,  6.
    , 12.
    ,  1.
    ,  6.
    ,  2.
    ,  0.
    ,  1.
    ,  1.
    ,  0.
    , 24.
    ,  0.
    ,  3.
    ,  1.
    ,  2.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.527111701601660
    , -0.401317830052742e2
    ,  0.153020073134484e3
    , -0.224799398218827e4
    , -0.193993484669048
    , -0.140467557893768e1
    ,  0.426799878114024e2
    ,  0.752810643416743
    ,  0.226657238616417e2
    , -0.622873556909932e3
    , -0.660823667935396
    ,  0.841267087271658
    , -0.253717501764397e2
    ,  0.485708963532948e3
    ,  0.880531517490555e3
    ,  0.265015592794626e7
    , -0.359287150025783
    , -0.656991567673753e3
    ,  0.241768149185367e1
    ,  0.856873461222588
    ,  0.655143675313458
    , -0.213535213206406
    ,  0.562974957606348e-2
    , -0.316955725450471e15
    , -0.699997000152457e-3
    ,  0.119845803210767e-1
    ,  0.193848122022095e-4
    , -0.215095749182309e-4
}; //}}}
} // namespace b
} // namespace t_ps
    constexpr auto
temperature_ps (auto const& pressure, auto const& massic_entropy)
{
        using namespace t_ps;
        auto const
    pi = pressure / (100e6);
    if (massic_entropy < 4.41202148223476e3)
    {
            auto const
        sigma = massic_entropy / (4.4e3);
        return sum (a::n * pow (pi + 0.240, a::I) * pow (sigma - 0.703, a::J)) * (760);
    }
        auto const
    sigma = massic_entropy / (5.3e3);
    return sum (b::n * pow (pi + 0.760, b::I) * pow (sigma - 0.818, b::J)) * (860);
}
    namespace
v_ps
{
    namespace
a
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    , -12.
    , -10.
    , -10.
    , -10.
    , -10.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -8.
    ,  -6.
    ,  -5.
    ,  -4.
    ,  -3.
    ,  -3.
    ,  -2.
    ,  -2.
    ,  -1.
    ,  -1.
    ,   0.
    ,   0.
    ,   0.
    ,   1.
    ,   2.
    ,   4.
    ,   5.
    ,   6.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      10.
    , 12.
    , 14.
    ,  4.
    ,  8.
    , 10.
    , 20.
    ,  5.
    ,  6.
    , 14.
    , 16.
    , 28.
    ,  1.
    ,  5.
    ,  2.
    ,  4.
    ,  3.
    ,  8.
    ,  1.
    ,  2.
    ,  0.
    ,  1.
    ,  3.
    ,  0.
    ,  0.
    ,  2.
    ,  2.
    ,  0.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.795544074093975e2
    , -0.238261242984590e4
    ,  0.176813100617787e5
    , -0.110524727080379e-2
    , -0.153213833655326e2
    ,  0.297544599376982e3
    , -0.350315206871242e8
    ,  0.277513761062119
    , -0.523964271036888
    , -0.148011182995403e6
    ,  0.160014899374266e7
    ,  0.170802322663427e13
    ,  0.246866996006494e-3
    ,  0.165326084797980e1
    , -0.118008384666987
    ,  0.253798642355900e1
    ,  0.965127704669424
    , -0.282172420532826e2
    ,  0.203224612353823
    ,  0.110648186063513e1
    ,  0.526127948451280
    ,  0.277000018736321
    ,  0.108153340501132e1
    , -0.744127885357893e-1
    ,  0.164094443541384e-1
    , -0.680468275301065e-1
    ,  0.257988576101640e-1
    , -0.145749861944416e-3
}; //}}}
} // namespace a
    namespace
b
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    , -12.
    , -12.
    , -12.
    , -12.
    , -10.
    , -10.
    , -10.
    , -10.
    ,  -8.
    ,  -5.
    ,  -5.
    ,  -5.
    ,  -4.
    ,  -4.
    ,  -4.
    ,  -4.
    ,  -3.
    ,  -2.
    ,  -2.
    ,  -2.
    ,  -2.
    ,  -2.
    ,  -2.
    ,   0.
    ,   0.
    ,   0.
    ,   1.
    ,   1.
    ,   2.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    ,  2.
    ,  3.
    ,  5.
    ,  6.
    ,  0.
    ,  1.
    ,  2.
    ,  4.
    ,  0.
    ,  1.
    ,  2.
    ,  3.
    ,  0.
    ,  1.
    ,  2.
    ,  3.
    ,  1.
    ,  0.
    ,  1.
    ,  2.
    ,  3.
    ,  4.
    , 12.
    ,  0.
    ,  1.
    ,  2.
    ,  0.
    ,  2.
    ,  2.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.591599780322238e-4
    , -0.185465997137856e-2
    ,  0.104190510480013e-1
    ,  0.598647302038590e-2
    , -0.771391189901699
    ,  0.172549765557036e1
    , -0.467076079846526e-3
    ,  0.134533823384439e-1
    , -0.808094336805495e-1
    ,  0.508139374365767
    ,  0.128584643361683e-2
    , -0.163899353915435e1
    ,  0.586938199318063e1
    , -0.292466667918613e1
    , -0.614076301499537e-2
    ,  0.576199014049172e1
    , -0.121613320606788e2
    ,  0.167637540957944e1
    , -0.744135838773463e1
    ,  0.378168091437659e-1
    ,  0.401432203027688e1
    ,  0.160279837479185e2
    ,  0.317848779347728e1
    , -0.358362310304853e1
    , -0.115995260446827e7
    ,  0.199256573577909
    , -0.122270624794624
    , -0.191449143716586e2
    , -0.150448002905284e-1
    ,  0.146407900162154e2
    , -0.327477787188230e1
}; //}}}
} // namespace b
} // namespace v_ps
    constexpr auto
massic_volume_ps (auto const& pressure, auto const& massic_entropy)
{
        using namespace v_ps;
        auto const
    pi = pressure / (100e6);
    if (massic_entropy < 4.41202148223476e3)
    {
            auto const
        sigma = massic_entropy / (4.4e3);
        return sum (a::n * pow (pi + 0.187, a::I) * pow (sigma - 0.755, a::J)) * (0.0028);
    }
        auto const
    sigma = massic_entropy / (5.3e3);
    return sum (b::n * pow (pi + 0.298, b::I) * pow (sigma - 0.816, b::J)) * (0.0088);
}
    namespace
p_hs
{
    namespace
a
{
    constexpr auto
I = array_t //{{{
{
        0.
    ,   0.
    ,   0.
    ,   1.
    ,   1.
    ,   1.
    ,   1.
    ,   1.
    ,   2.
    ,   2.
    ,   3.
    ,   3.
    ,   3.
    ,   4.
    ,   4.
    ,   4.
    ,   4.
    ,   5.
    ,   6.
    ,   7.
    ,   8.
    ,  10.
    ,  10.
    ,  14.
    ,  18.
    ,  20.
    ,  22.
    ,  22.
    ,  24.
    ,  28.
    ,  28.
    ,  32.
    ,  32.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
    ,  1.
    ,  5.
    ,  0.
    ,  3.
    ,  4.
    ,  8.
    , 14.
    ,  6.
    , 16.
    ,  0.
    ,  2.
    ,  3.
    ,  0.
    ,  1.
    ,  4.
    ,  5.
    , 28.
    , 28.
    , 24.
    ,  1.
    , 32.
    , 36.
    , 22.
    , 28.
    , 36.
    , 16.
    , 28.
    , 36.
    , 16.
    , 36.
    , 10.
    , 28.
}; //}}}

    constexpr auto
n = array_t //{{{
{
       0.770889828326934e1
    , -0.260835009128688e2
    ,  0.267416218930389e3
    ,  0.172221089496844e2
    , -0.293542332145970e3
    ,  0.614135601882478e3
    , -0.610562757725674e5
    , -0.651272251118219e8
    ,  0.735919313521937e5
    , -0.116646505914191e11
    ,  0.355267086434461e2
    , -0.596144543825955e3
    , -0.475842430145708e3
    ,  0.696781965359503e2
    ,  0.335674250377312e3
    ,  0.250526809130882e5
    ,  0.146997380630766e6
    ,  0.538069315091534e20
    ,  0.143619827291346e22
    ,  0.364985866165994e20
    , -0.254741561156775e4
    ,  0.240120197096563e28
    , -0.393847464679496e30
    ,  0.147073407024852e25
    , -0.426391250432059e32
    ,  0.194509340621077e39
    ,  0.666212132114896e24
    ,  0.706777016552858e34
    ,  0.175563621975576e42
    ,  0.108408607429124e29
    ,  0.730872705175151e44
    ,  0.159145847398870e25
    ,  0.377121605943324e41
}; //}}}
} // namespace a
    namespace
b
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -12.
    , -12.
    , -12.
    , -12.
    , -10.
    , -10.
    , -10.
    , -10.
    ,  -8.
    ,  -8.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -6.
    ,  -5.
    ,  -4.
    ,  -4.
    ,  -4.
    ,  -3.
    ,  -3.
    ,  -3.
    ,  -3.
    ,  -2.
    ,  -2.
    ,  -1.
    ,   0.
    ,   2.
    ,   2.
    ,   5.
    ,   6.
    ,   8.
    ,  10.
    ,  14.
    ,  14.
}; //}}}

    constexpr auto
J = array_t //{{{
{
       2.
    , 10.
    , 12.
    , 14.
    , 20.
    ,  2.
    , 10.
    , 14.
    , 18.
    ,  2.
    ,  8.
    ,  2.
    ,  6.
    ,  7.
    ,  8.
    , 10.
    ,  4.
    ,  5.
    ,  8.
    ,  1.
    ,  3.
    ,  5.
    ,  6.
    ,  0.
    ,  1.
    ,  0.
    ,  3.
    ,  0.
    ,  1.
    ,  0.
    ,  1.
    ,  1.
    ,  1.
    ,  3.
    ,  7.
}; //}}}

    constexpr auto
n = array_t //{{{
{
       0.125244360717979e-12
    , -0.126599322553713e-1
    ,  0.506878030140626e1
    ,  0.317847171154202e2
    , -0.391041161399932e6
    , -0.975733406392044e-10
    , -0.186312419488279e2
    ,  0.510973543414101e3
    ,  0.373847005822362e6
    ,  0.299804024666572e-7
    ,  0.200544393820342e2
    , -0.498030487662829e-5
    , -0.102301806360030e2
    ,  0.552819126990325e2
    , -0.206211367510878e3
    , -0.794012232324823e4
    ,  0.782248472028153e1
    , -0.586544326902468e2
    ,  0.355073647696481e4
    , -0.115303107290162e-3
    , -0.175092403171802e1
    ,  0.257981687748160e3
    , -0.727048374179467e3
    ,  0.121644822609198e-3
    ,  0.393137871762692e-1
    ,  0.704181005909296e-2
    , -0.829108200698110e2
    , -0.265178818131250
    ,  0.137531682453991e2
    , -0.522394090753046e2
    ,  0.240556298941048e4
    , -0.227361631268929e5
    ,  0.890746343932567e5
    , -0.239234565822486e8
    ,  0.568795808129714e10
}; //}}}
} // namespace b
} // namespace p_hs
    constexpr auto
pressure_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
        using namespace p_hs;
    if (massic_entropy <= (4.41202148223476e3))
    {
            auto const
        eta = massic_enthalpy / (2300e3);
            auto const
        sigma = massic_entropy / (4.4e3);
        return sum (a::n * pow (eta - 1.01, a::I) * pow (sigma - 0.750, a::J)) * (99e6);
    }
        auto const
    eta = massic_enthalpy / (2800e3);
        auto const
    sigma = massic_entropy / (5.3e3);
    return (16.6e6) / sum (b::n * pow (eta - 0.681, b::I) * pow (sigma - 0.792, b::J));
}
    constexpr auto
temperature_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    return temperature_ph (pressure_hs (massic_enthalpy, massic_entropy), massic_enthalpy);
}
    constexpr auto
massic_volume_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    return massic_volume_ps (pressure_hs (massic_enthalpy, massic_entropy), massic_entropy);
}
#define ISTO_IAPWS_GEN_T3_1(X)                         \
    constexpr auto                                     \
t_3##X##_p (auto const& pressure)            \
{                                                      \
        using namespace t_3##X;                        \
        auto const                                     \
    pi = pressure / 1e6;              \
        using std::pow;                                \
    return sum (n * pow (pi, I)) * (1); \
}
#define ISTO_IAPWS_GEN_T3_2(X)                               \
    constexpr auto                                           \
t_3##X##_p (auto const& pressure)                    \
{                                                            \
        using namespace t_3##X;                              \
        auto const                                           \
    pi = pressure / 1e6;                    \
        using std::pow, std::log;                            \
    return sum (n * pow (log (pi), I)) * (1); \
}
    namespace
t_3ab
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        , -1.
        , -2.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.154793642129415e4
        , -0.187661219490113e3
        ,  0.213144632222113e2
        , -0.191887498864292e4
        ,  0.918419702359447e3
    }; //}}}
} // namespace t_3ab
ISTO_IAPWS_GEN_T3_2(ab);
    namespace
t_3cd
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        ,  3.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.585276966696349e3
        ,  0.278233532206915e1
        , -0.127283549295878e-1
        ,  0.159090746562729e-3
    }; //}}}
} // namespace t_3cd
ISTO_IAPWS_GEN_T3_1(cd);
    constexpr auto
t_3ef_p (auto const& pressure)
{
        auto const
    pi = pressure / 1e6;
    return (3.727888004 * (pi - 22.064) + 647.096) * (1);
}
    namespace
t_3gh
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        ,  3.
        ,  4.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
          -0.249284240900418e5
        ,  0.428143584791546e4
        , -0.269029173140130e3
        ,  0.751608051114157e1
        , -0.787105249910383e-1
    }; //}}}
} // namespace t_3gh
ISTO_IAPWS_GEN_T3_1(gh);
    namespace
t_3ij
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        ,  3.
        ,  4.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.584814781649163e3
        , -0.616179320924617
        ,  0.260763050899562
        , -0.587071076864459e-2
        ,  0.515308185433082e-4
    }; //}}}
} // namespace t_3ij
ISTO_IAPWS_GEN_T3_1(ij);
    namespace
t_3jk
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        ,  3.
        ,  4.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.617229772068439e3
        , -0.770600270141675e1
        ,  0.697072596851896
        , -0.157391839848015e-1
        ,  0.137897492684194e-3
    }; //}}}
} // namespace t_3jk
ISTO_IAPWS_GEN_T3_1(jk);
    namespace
t_3mn
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        ,  3.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.535339483742384e3
        ,  0.761978122720128e1
        , -0.158365725441648
        ,  0.192871054508108e-2
    }; //}}}
} // namespace t_3mn
ISTO_IAPWS_GEN_T3_1(mn);
    namespace
t_3op
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        , -1.
        , -2.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.969461372400213e3
        , -0.332500170441278e3
        ,  0.642859598466067e2
        ,  0.773845935768222e3
        , -0.152313732937084e4
    }; //}}}
} // namespace t_3op
ISTO_IAPWS_GEN_T3_2(op);
    namespace
t_3qu
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        ,  3.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.565603648239126e3
        ,  0.529062258221222e1
        , -0.102020639611016
        ,  0.122240301070145e-2
    }; //}}}
} // namespace t_3qu
ISTO_IAPWS_GEN_T3_1(qu);
    namespace
t_3rx
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        ,  3.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.584561202520006e3
        , -0.102961025163669e1
        ,  0.243293362700452
        , -0.294905044740799e-2
    }; //}}}
} // namespace t_3rx
ISTO_IAPWS_GEN_T3_1(rx);
    namespace
t_3uv
{
        constexpr auto
    I = array_t //{{{
    {
          0.
        , 1.
        , 2.
        , 3.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
           0.528199646263062e3
        ,  0.890579602135307e1
        , -0.222814134903755
        ,  0.286791682263697e-2
    }; //}}}
} // namespace t_3uv
ISTO_IAPWS_GEN_T3_1(uv);
    namespace
t_3wx
{
        constexpr auto
    I = array_t //{{{
    {
           0.
        ,  1.
        ,  2.
        , -1.
        , -2.
    }; //}}}
        constexpr auto
    n = array_t //{{{
    {
          0.728052609145380e1
        , 0.973505869861952e2
        , 0.147370491183191e2
        , 0.329196213998375e3
        , 0.873371668682417e3
    }; //}}}
} // namespace t_3wx
ISTO_IAPWS_GEN_T3_2(wx);
#undef ISTO_IAPWS_GEN_T3_1
#undef ISTO_IAPWS_GEN_T3_2

    namespace
v_3a
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -12.
	, -10.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -6.
	, -5.
	, -5.
	, -5.
	, -4.
	, -3.
	, -3.
	, -3.
	, -3.
	, -2.
	, -2.
	, -2.
	, -1.
	, -1.
	, -1.
	,  0.
	,  0.
	,  1.
	,  1.
	,  2.
	,  2.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       5.
	,  10.
	,  12.
	,  5.
	,  10.
	,  12.
	,  5.
	,  8.
	,  10.
	,  1.
	,  1.
	,  5.
	,  10.
	,  8.
	,  0.
	,  1.
	,  3.
	,  6.
	,  0.
	,  2.
	,  3.
	,  0.
	,  1.
	,  2.
	,  0.
	,  1.
	,  0.
	,  2.
	,  0.
	,  2.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.110879558823853e-2
	,  0.572616740810616e3
	, -0.767051948380852e5
	, -0.253321069529674e-1
	,  0.628008049345689e4
	,  0.234105654131876e6
	,  0.216867826045856
	, -0.156237904341963e3
	, -0.269893956176613e5
	, -0.180407100085505e-3
	,  0.116732227668261e-2
	,  0.266987040856040e2
	,  0.282776617243286e5
	, -0.242431520029523e4
	,  0.435217323022733e-3
	, -0.122494831387441e-1
	,  0.179357604019989e1
	,  0.442729521058314e2
	, -0.593223489018342e-2
	,  0.453186261685774
	,  0.135825703129140e1
	,  0.408748415856745e-1
	,  0.474686397863312
	,  0.118646814997915e1
	,  0.546987265727549
	,  0.195266770452643
	, -0.502268790869663e-1
	, -0.369645308193377
	,  0.633828037528420e-2
	,  0.797441793901017e-1
}; //}}}
    constexpr auto v_s = 0.0024;
    constexpr auto p_s = 100e6 ;
    constexpr auto t_s = 760.  ;
    constexpr auto N   = 30.;
    constexpr auto a   = 0.085;
    constexpr auto b   = 0.817;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3a
    namespace
v_3b
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -10.
	, -10.
	, -8.
	, -6.
	, -6.
	, -6.
	, -5.
	, -5.
	, -5.
	, -4.
	, -4.
	, -4.
	, -3.
	, -3.
	, -3.
	, -3.
	, -3.
	, -2.
	, -2.
	, -2.
	, -1.
	, -1.
	,  0.
	,  0.
	,  1.
	,  1.
	,  2.
	,  3.
	,  4.
	,  4.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       10.
	,  12.
	,  8.
	,  14.
	,  8.
	,  5.
	,  6.
	,  8.
	,  5.
	,  8.
	,  10.
	,  2.
	,  4.
	,  5.
	,  0.
	,  1.
	,  2.
	,  3.
	,  5.
	,  0.
	,  2.
	,  5.
	,  0.
	,  2.
	,  0.
	,  1.
	,  0.
	,  2.
	,  0.
	,  2.
	,  0.
	,  1.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.827670470003621e-1
	,  0.416887126010565e2
	,  0.483651982197059e-1
	, -0.291032084950276e5
	, -0.111422582236948e3
	, -0.202300083904014e-1
	,  0.294002509338515e3
	,  0.140244997609658e3
	, -0.344384158811459e3
	,  0.361182452612149e3
	, -0.140699677420738e4
	, -0.202023902676481e-2
	,  0.171346792457471e3
	, -0.425597804058632e1
	,  0.691346085000334e-5
	,  0.151140509678925e-2
	, -0.416375290166236e-1
	, -0.413754957011042e2
	, -0.506673295721637e2
	, -0.572212965569023e-3
	,  0.608817368401785e1
	,  0.239600660256161e2
	,  0.122261479925384e-1
	,  0.216356057692938e1
	,  0.398198903368642
	, -0.116892827834085
	, -0.102845919373532
	, -0.492676637589284
	,  0.655540456406790e-1
	, -0.240462535078530
	, -0.269798180310075e-1
	,  0.128369435967012
}; //}}}
    constexpr auto v_s = 0.0041;
    constexpr auto p_s = 100e6 ;
    constexpr auto t_s = 860.  ;
    constexpr auto N   = 32.;
    constexpr auto a   = 0.280;
    constexpr auto b   = 0.779;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3b
    namespace
v_3c
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -12.
	, -10.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -6.
	, -5.
	, -5.
	, -5.
	, -4.
	, -4.
	, -3.
	, -3.
	, -2.
	, -2.
	, -2.
	, -1.
	, -1.
	, -1.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  2.
	,  2.
	,  2.
	,  2.
	,  3.
	,  3.
	,  8.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       6.
	,  8.
	,  10.
	,  6.
	,  8.
	,  10.
	,  5.
	,  6.
	,  7.
	,  8.
	,  1.
	,  4.
	,  7.
	,  2.
	,  8.
	,  0.
	,  3.
	,  0.
	,  4.
	,  5.
	,  0.
	,  1.
	,  2.
	,  0.
	,  1.
	,  2.
	,  0.
	,  2.
	,  0.
	,  1.
	,  3.
	,  7.
	,  0.
	,  7.
	,  1.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.311967788763030e1
	,  0.276713458847564e5
	,  0.322583103403269e8
	, -0.342416065095363e3
	, -0.899732529907377e6
	, -0.793892049821251e8
	,  0.953193003217388e2
	,  0.229784742345072e4
	,  0.175336675322499e6
	,  0.791214365222792e7
	,  0.319933345844209e-4
	, -0.659508863555767e2
	, -0.833426563212851e6
	,  0.645734680583292e-1
	, -0.382031020570813e7
	,  0.406398848470079e-4
	,  0.310327498492008e2
	, -0.892996718483724e-3
	,  0.234604891591616e3
	,  0.377515668966951e4
	,  0.158646812591361e-1
	,  0.707906336241843
	,  0.126016225146570e2
	,  0.736143655772152
	,  0.676544268999101
	, -0.178100588189137e2
	, -0.156531975531713
	,  0.117707430048158e2
	,  0.840143653860447e-1
	, -0.186442467471949
	, -0.440170203949645e2
	,  0.123290423502494e7
	, -0.240650039730845e-1
	, -0.107077716660869e7
	,  0.438319858566475e-1
}; //}}}
    constexpr auto v_s = 0.0022;
    constexpr auto p_s = 40e6  ;
    constexpr auto t_s = 690.  ;
    constexpr auto N   = 35.;
    constexpr auto a   = 0.259;
    constexpr auto b   = 0.903;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3c
    namespace
v_3d
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -10.
	, -10.
	, -10.
	, -10.
	, -10.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -8.
	, -6.
	, -6.
	, -5.
	, -5.
	, -5.
	, -5.
	, -4.
	, -4.
	, -4.
	, -3.
	, -3.
	, -2.
	, -2.
	, -1.
	, -1.
	, -1.
	,  0.
	,  0.
	,  1.
	,  1.
	,  3.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       4.
	,  6.
	,  7.
	,  10.
	,  12.
	,  16.
	,  0.
	,  2.
	,  4.
	,  6.
	,  8.
	,  10.
	,  14.
	,  3.
	,  7.
	,  8.
	,  10.
	,  6.
	,  8.
	,  1.
	,  2.
	,  5.
	,  7.
	,  0.
	,  1.
	,  7.
	,  2.
	,  4.
	,  0.
	,  1.
	,  0.
	,  1.
	,  5.
	,  0.
	,  2.
	,  0.
	,  6.
	,  0.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.452484847171645e-9
	,  0.315210389538801e-4
	, -0.214991352047545e-2
	,  0.508058874808345e3
	, -0.127123036845932e8
	,  0.115371133120497e13
	, -0.197805728776273e-15
	,  0.241554806033972e-10
	, -0.156481703640525e-5
	,  0.277211346836625e-2
	, -0.203578994462286e2
	,  0.144369489909053e7
	, -0.411254217946539e11
	,  0.623449786243773e-5
	, -0.221774281146038e2
	, -0.689315087933158e5
	, -0.195419525060713e8
	,  0.316373510564015e4
	,  0.224040754426988e7
	, -0.436701347922356e-5
	, -0.404213852833996e-3
	, -0.348153203414663e3
	, -0.385294213555289e6
	,  0.135203700099403e-6
	,  0.134648383271089e-3
	,  0.125031835351736e6
	,  0.968123678455841e-1
	,  0.225660517512438e3
	, -0.190102435341872e-3
	, -0.299628410819229e-1
	,  0.500833915372121e-2
	,  0.387842482998411
	, -0.138535367777182e4
	,  0.870745245971773
	,  0.171946252068742e1
	, -0.326650121426383e-1
	,  0.498044171727877e4
	,  0.551478022765087e-2
}; //}}}
    constexpr auto v_s = 0.0029;
    constexpr auto p_s = 40e6  ;
    constexpr auto t_s = 690.  ;
    constexpr auto N   = 38.;
    constexpr auto a   = 0.559;
    constexpr auto b   = 0.939;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3d
    namespace
v_3e
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -10.
	, -10.
	, -10.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -6.
	, -5.
	, -4.
	, -4.
	, -3.
	, -3.
	, -3.
	, -2.
	, -2.
	, -2.
	, -2.
	, -1.
	,  0.
	,  0.
	,  1.
	,  1.
	,  1.
	,  2.
	,  2.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       14.
	,  16.
	,  3.
	,  6.
	,  10.
	,  14.
	,  16.
	,  7.
	,  8.
	,  10.
	,  6.
	,  6.
	,  2.
	,  4.
	,  2.
	,  6.
	,  7.
	,  0.
	,  1.
	,  3.
	,  4.
	,  0.
	,  0.
	,  1.
	,  0.
	,  4.
	,  6.
	,  0.
	,  2.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.715815808404721e9
	, -0.114328360753449e12
	,  0.376531002015720e-11
	, -0.903983668691157e-4
	,  0.665695908836252e6
	,  0.535364174960127e10
	,  0.794977402335603e11
	,  0.922230563421437e2
	, -0.142586073991215e6
	, -0.111796381424162e7
	,  0.896121629640760e4
	, -0.669989239070491e4
	,  0.451242538486834e-2
	, -0.339731325977713e2
	, -0.120523111552278e1
	,  0.475992667717124e5
	, -0.266627750390341e6
	, -0.153314954386524e-3
	,  0.305638404828265
	,  0.123654999499486e3
	, -0.104390794213011e4
	, -0.157496516174308e-1
	,  0.685331118940253
	,  0.178373462873903e1
	, -0.544674124878910
	,  0.204529931318843e4
	, -0.228342359328752e5
	,  0.413197481515899
	, -0.341931835910405e2
}; //}}}
    constexpr auto v_s = 0.0032;
    constexpr auto p_s = 40e6  ;
    constexpr auto t_s = 710.  ;
    constexpr auto N   = 29.;
    constexpr auto a   = 0.587;
    constexpr auto b   = 0.918;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3e
    namespace
v_3f
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  0.
	,  0.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  1.
	,  1.
	,  2.
	,  2.
	,  3.
	,  3.
	,  3.
	,  4.
	,  5.
	,  5.
	,  6.
	,  7.
	,  7.
	,  10.
	,  12.
	,  12.
	,  12.
	,  14.
	,  14.
	,  14.
	,  14.
	,  14.
	,  16.
	,  16.
	,  18.
	,  18.
	,  20.
	,  20.
	,  20.
	,  22.
	,  24.
	,  24.
	,  28.
	,  32.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      -3.
	, -2.
	, -1.
	,  0.
	,  1.
	,  2.
	, -1.
	,  1.
	,  2.
	,  3.
	,  0.
	,  1.
	, -5.
	, -2.
	,  0.
	, -3.
	, -8.
	,  1.
	, -6.
	, -4.
	,  1.
	, -6.
	, -10.
	, -8.
	, -4.
	, -12.
	, -10.
	, -8.
	, -6.
	, -4.
	, -10.
	, -8.
	, -12.
	, -10.
	, -12.
	, -10.
	, -6.
	, -12.
	, -12.
	, -4.
	, -12.
	, -12.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.251756547792325e-7
	,  0.601307193668763e-5
	, -0.100615977450049e-2
	,  0.999969140252192
	,  0.214107759236486e1
	, -0.165175571959086e2
	, -0.141987303638727e-2
	,  0.269251915156554e1
	,  0.349741815858722e2
	, -0.300208695771783e2
	, -0.131546288252539e1
	, -0.839091277286169e1
	,  0.181545608337015e-9
	, -0.591099206478909e-3
	,  0.152115067087106e1
	,  0.252956470663225e-4
	,  0.100726265203786e-14
	, -0.149774533860650e1
	, -0.793940970562969e-9
	, -0.150290891264717e-3
	,  0.151205531275133e1
	,  0.470942606221652e-5
	,  0.195049710391712e-12
	, -0.911627886266077e-8
	,  0.604374640201265e-3
	, -0.225132933900136e-15
	,  0.610916973582981e-11
	, -0.303063908043404e-6
	, -0.137796070798409e-4
	, -0.919296736666106e-3
	,  0.639288223132545e-9
	,  0.753259479898699e-6
	, -0.400321478682929e-12
	,  0.756140294351614e-8
	, -0.912082054034891e-11
	, -0.237612381140539e-7
	,  0.269586010591874e-4
	, -0.732828135157839e-10
	,  0.241995578306660e-9
	, -0.405735532730322e-3
	,  0.189424143498011e-9
	, -0.486632965074563e-9
}; //}}}
    constexpr auto v_s = 0.0064;
    constexpr auto p_s = 40e6  ;
    constexpr auto t_s = 730.  ;
    constexpr auto N   = 42.;
    constexpr auto a   = 0.587;
    constexpr auto b   = 0.891;
    constexpr auto c   = 0.5;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3f
    namespace
v_3g
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -10.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -8.
	, -6.
	, -6.
	, -5.
	, -5.
	, -4.
	, -3.
	, -2.
	, -2.
	, -2.
	, -2.
	, -1.
	, -1.
	, -1.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  1.
	,  3.
	,  5.
	,  6.
	,  8.
	,  10.
	,  10.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       7.
	,  12.
	,  14.
	,  18.
	,  22.
	,  24.
	,  14.
	,  20.
	,  24.
	,  7.
	,  8.
	,  10.
	,  12.
	,  8.
	,  22.
	,  7.
	,  20.
	,  22.
	,  7.
	,  3.
	,  5.
	,  14.
	,  24.
	,  2.
	,  8.
	,  18.
	,  0.
	,  1.
	,  2.
	,  0.
	,  1.
	,  3.
	,  24.
	,  22.
	,  12.
	,  3.
	,  0.
	,  6.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.412209020652996e-4
	, -0.114987238280587e7
	,  0.948180885032080e10
	, -0.195788865718971e18
	,  0.496250704871300e25
	, -0.105549884548496e29
	, -0.758642165988278e12
	, -0.922172769596101e23
	,  0.725379072059348e30
	, -0.617718249205859e2
	,  0.107555033344858e5
	, -0.379545802336487e8
	,  0.228646846221831e12
	, -0.499741093010619e7
	, -0.280214310054101e31
	,  0.104915406769586e7
	,  0.613754229168619e28
	,  0.802056715528378e32
	, -0.298617819828065e8
	, -0.910782540134681e2
	,  0.135033227281565e6
	, -0.712949383408211e19
	, -0.104578785289542e37
	,  0.304331584444093e2
	,  0.593250797959445e10
	, -0.364174062110798e28
	,  0.921791403532461
	, -0.337693609657471
	, -0.724644143758508e2
	, -0.110480239272601
	,  0.536516031875059e1
	, -0.291441872156205e4
	,  0.616338176535305e40
	, -0.120889175861180e39
	,  0.818396024524612e23
	,  0.940781944835829e9
	, -0.367279669545448e5
	, -0.837513931798655e16
}; //}}}
    constexpr auto v_s = 0.0027;
    constexpr auto p_s = 25e6  ;
    constexpr auto t_s = 660.  ;
    constexpr auto N   = 38.;
    constexpr auto a   = 0.872;
    constexpr auto b   = 0.971;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3g
    namespace
v_3h
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -10.
	, -10.
	, -10.
	, -10.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -8.
	, -8.
	, -6.
	, -6.
	, -6.
	, -5.
	, -5.
	, -5.
	, -4.
	, -4.
	, -3.
	, -3.
	, -2.
	, -1.
	, -1.
	,  0.
	,  1.
	,  1.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       8.
	,  12.
	,  4.
	,  6.
	,  8.
	,  10.
	,  14.
	,  16.
	,  0.
	,  1.
	,  6.
	,  7.
	,  8.
	,  4.
	,  6.
	,  8.
	,  2.
	,  3.
	,  4.
	,  2.
	,  4.
	,  1.
	,  2.
	,  0.
	,  0.
	,  2.
	,  0.
	,  0.
	,  2.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.561379678887577e-1
	,  0.774135421587083e10
	,  0.111482975877938e-8
	, -0.143987128208183e-2
	,  0.193696558764920e4
	, -0.605971823585005e9
	,  0.171951568124337e14
	, -0.185461154985145e17
	,  0.387851168078010e-16
	, -0.395464327846105e-13
	, -0.170875935679023e3
	, -0.212010620701220e4
	,  0.177683337348191e8
	,  0.110177443629575e2
	, -0.234396091693313e6
	, -0.656174421999594e7
	,  0.156362212977396e-4
	, -0.212946257021400e1
	,  0.135249306374858e2
	,  0.177189164145813
	,  0.139499167345464e4
	, -0.703670932036388e-2
	, -0.152011044389648
	,  0.981916922991113e-4
	,  0.147199658618076e-2
	,  0.202618487025578e2
	,  0.899345518944240
	, -0.211346402240858
	,  0.249971752957491e2
}; //}}}
    constexpr auto v_s = 0.0032;
    constexpr auto p_s = 25e6  ;
    constexpr auto t_s = 660.  ;
    constexpr auto N   = 29.;
    constexpr auto a   = 0.898;
    constexpr auto b   = 0.983;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3h
    namespace
v_3i
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  1.
	,  1.
	,  2.
	,  3.
	,  3.
	,  4.
	,  4.
	,  4.
	,  5.
	,  5.
	,  5.
	,  7.
	,  7.
	,  8.
	,  8.
	,  10.
	,  12.
	,  12.
	,  12.
	,  14.
	,  14.
	,  14.
	,  14.
	,  18.
	,  18.
	,  18.
	,  18.
	,  18.
	,  20.
	,  20.
	,  22.
	,  24.
	,  24.
	,  32.
	,  32.
	,  36.
	,  36.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
	,  1.
	,  10.
	, -4.
	, -2.
	, -1.
	,  0.
	,  0.
	, -5.
	,  0.
	, -3.
	, -2.
	, -1.
	, -6.
	, -1.
	,  12.
	, -4.
	, -3.
	, -6.
	,  10.
	, -8.
	, -12.
	, -6.
	, -4.
	, -10.
	, -8.
	, -4.
	,  5.
	, -12.
	, -10.
	, -8.
	, -6.
	,  2.
	, -12.
	, -10.
	, -12.
	, -12.
	, -8.
	, -10.
	, -5.
	, -10.
	, -8.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.106905684359136e1
	, -0.148620857922333e1
	,  0.259862256980408e15
	, -0.446352055678749e-11
	, -0.566620757170032e-6
	, -0.235302885736849e-2
	, -0.269226321968839
	,  0.922024992944392e1
	,  0.357633505503772e-11
	, -0.173942565562222e2
	,  0.700681785556229e-5
	, -0.267050351075768e-3
	, -0.231779669675624e1
	, -0.753533046979752e-12
	,  0.481337131452891e1
	, -0.223286270422356e22
	, -0.118746004987383e-4
	,  0.646412934136496e-2
	, -0.410588536330937e-9
	,  0.422739537057241e20
	,  0.313698180473812e-12
	,  0.164395334345040e-23
	, -0.339823323754373e-5
	, -0.135268639905021e-1
	, -0.723252514211625e-14
	,  0.184386437538366e-8
	, -0.463959533752385e-1
	, -0.992263100376750e14
	,  0.688169154439335e-16
	, -0.222620998452197e-10
	, -0.540843018624083e-7
	,  0.345570606200257e-2
	,  0.422275800304086e11
	, -0.126974478770487e-14
	,  0.927237985153679e-9
	,  0.612670812016489e-13
	, -0.722693924063497e-11
	, -0.383669502636822e-3
	,  0.374684572410204e-3
	, -0.931976897511086e5
	, -0.247690616026922e-1
	,  0.658110546759474e2
}; //}}}
    constexpr auto v_s = 0.0041;
    constexpr auto p_s = 25e6  ;
    constexpr auto t_s = 660.  ;
    constexpr auto N   = 42.;
    constexpr auto a   = 0.910;
    constexpr auto b   = 0.984;
    constexpr auto c   = 0.5;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3i
    namespace
v_3j
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  1.
	,  2.
	,  2.
	,  3.
	,  4.
	,  4.
	,  5.
	,  5.
	,  5.
	,  6.
	,  10.
	,  12.
	,  12.
	,  14.
	,  14.
	,  14.
	,  16.
	,  18.
	,  20.
	,  20.
	,  24.
	,  24.
	,  28.
	,  28.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      -1.
	,  0.
	,  1.
	, -2.
	, -1.
	,  1.
	, -1.
	,  1.
	, -2.
	, -2.
	,  2.
	, -3.
	, -2.
	,  0.
	,  3.
	, -6.
	, -8.
	, -3.
	, -10.
	, -8.
	, -5.
	, -10.
	, -12.
	, -12.
	, -10.
	, -12.
	, -6.
	, -12.
	, -5.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.111371317395540e-3
	,  0.100342892423685e1
	,  0.530615581928979e1
	,  0.179058760078792e-5
	, -0.728541958464774e-3
	, -0.187576133371704e2
	,  0.199060874071849e-2
	,  0.243574755377290e2
	, -0.177040785499444e-3
	, -0.259680385227130e-2
	, -0.198704578406823e3
	,  0.738627790224287e-4
	, -0.236264692844138e-2
	, -0.161023121314333e1
	,  0.622322971786473e4
	, -0.960754116701669e-8
	, -0.510572269720488e-10
	,  0.767373781404211e-2
	,  0.663855469485254e-14
	, -0.717590735526745e-9
	,  0.146564542926508e-4
	,  0.309029474277013e-11
	, -0.464216300971708e-15
	, -0.390499637961161e-13
	, -0.236716126781431e-9
	,  0.454652854268717e-11
	, -0.422271787482497e-2
	,  0.283911742354706e-10
	,  0.270929002720228e1
}; //}}}
    constexpr auto v_s = 0.0054;
    constexpr auto p_s = 25e6  ;
    constexpr auto t_s = 670.  ;
    constexpr auto N   = 29.;
    constexpr auto a   = 0.875;
    constexpr auto b   = 0.964;
    constexpr auto c   = 0.5;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3j
    namespace
v_3k
{
    constexpr auto
I = array_t //{{{
{
      -2.
	, -2.
	, -1.
	, -1.
	,  0.
	,  0.
	,  0.
	,  0.
	,  0.
	,  0.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  1.
	,  1.
	,  1.
	,  2.
	,  2.
	,  2.
	,  2.
	,  2.
	,  2.
	,  5.
	,  5.
	,  5.
	,  6.
	,  6.
	,  6.
	,  6.
	,  8.
	,  10.
	,  12.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       10.
	,  12.
	, -5.
	,  6.
	, -12.
	, -6.
	, -2.
	, -1.
	,  0.
	,  1.
	,  2.
	,  3.
	,  14.
	, -3.
	, -2.
	,  0.
	,  1.
	,  2.
	, -8.
	, -6.
	, -3.
	, -2.
	,  0.
	,  4.
	, -12.
	, -6.
	, -3.
	, -12.
	, -10.
	, -8.
	, -5.
	, -12.
	, -12.
	, -10.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.401215699576099e9
	,  0.484501478318406e11
	,  0.394721471363678e-14
	,  0.372629967374147e5
	, -0.369794374168666e-29
	, -0.380436407012452e-14
	,  0.475361629970233e-6
	, -0.879148916140706e-3
	,  0.844317863844331
	,  0.122433162656600e2
	, -0.104529634830279e3
	,  0.589702771277429e3
	, -0.291026851164444e14
	,  0.170343072841850e-5
	, -0.277617606975748e-3
	, -0.344709605486686e1
	,  0.221333862447095e2
	, -0.194646110037079e3
	,  0.808354639772825e-15
	, -0.180845209145470e-10
	, -0.696664158132412e-5
	, -0.181057560300994e-2
	,  0.255830298579027e1
	,  0.328913873658481e4
	, -0.173270241249904e-18
	, -0.661876792558034e-6
	, -0.395688923421250e-2
	,  0.604203299819132e-17
	, -0.400879935920517e-13
	,  0.160751107464958e-8
	,  0.383719409025556e-4
	, -0.649565446702457e-14
	, -0.149095328506000e-11
	,  0.541449377329581e-8
}; //}}}
    constexpr auto v_s = 0.0077;
    constexpr auto p_s = 25e6  ;
    constexpr auto t_s = 680.  ;
    constexpr auto N   = 34.;
    constexpr auto a   = 0.802;
    constexpr auto b   = 0.935;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3k
    namespace
v_3l
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -8.
	, -8.
	, -8.
	, -8.
	, -6.
	, -5.
	, -5.
	, -4.
	, -4.
	, -3.
	, -3.
	, -3.
	, -3.
	, -2.
	, -2.
	, -2.
	, -1.
	, -1.
	, -1.
	,  0.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  2.
	,  4.
	,  5.
	,  5.
	,  6.
	,  10.
	,  10.
	,  14.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       14.
	,  16.
	,  18.
	,  20.
	,  22.
	,  14.
	,  24.
	,  6.
	,  10.
	,  12.
	,  14.
	,  18.
	,  24.
	,  36.
	,  8.
	,  4.
	,  5.
	,  7.
	,  16.
	,  1.
	,  3.
	,  18.
	,  20.
	,  2.
	,  3.
	,  10.
	,  0.
	,  1.
	,  3.
	,  0.
	,  1.
	,  2.
	,  12.
	,  0.
	,  16.
	,  1.
	,  0.
	,  0.
	,  1.
	,  14.
	,  4.
	,  12.
	,  10.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.260702058647537e10
	, -0.188277213604704e15
	,  0.554923870289667e19
	, -0.758966946387758e23
	,  0.413865186848908e27
	, -0.815038000738060e12
	, -0.381458260489955e33
	, -0.123239564600519e-1
	,  0.226095631437174e8
	, -0.495017809506720e12
	,  0.529482996422863e16
	, -0.444359478746295e23
	,  0.521635864527315e35
	, -0.487095672740742e55
	, -0.714430209937547e6
	,  0.127868634615495
	, -0.100752127917598e2
	,  0.777451437960990e7
	, -0.108105480796471e25
	, -0.357578581169659e-5
	, -0.212857169423484e1
	,  0.270706111085238e30
	, -0.695953622348829e33
	,  0.110609027472280
	,  0.721559163361354e2
	, -0.306367307532219e15
	,  0.265839618885530e-4
	,  0.253392392889754e-1
	, -0.214443041836579e3
	,  0.937846601489667
	,  0.223184043101700e1
	,  0.338401222509191e2
	,  0.494237237179718e21
	, -0.198068404154428
	, -0.141415349881140e31
	, -0.993862421613651e2
	,  0.125070534142731e3
	, -0.996473529004439e3
	,  0.473137909872765e5
	,  0.116662121219322e33
	, -0.315874976271533e16
	, -0.445703369196945e33
	,  0.642794932373694e33
}; //}}}
    constexpr auto v_s = 0.0026;
    constexpr auto p_s = 24e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 43.;
    constexpr auto a   = 0.908;
    constexpr auto b   = 0.989;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3l
    namespace
v_3m
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  3.
	,  8.
	,  20.
	,  1.
	,  3.
	,  4.
	,  5.
	,  1.
	,  6.
	,  2.
	,  4.
	,  14.
	,  2.
	,  5.
	,  3.
	,  0.
	,  1.
	,  1.
	,  1.
	,  28.
	,  2.
	,  16.
	,  0.
	,  5.
	,  0.
	,  3.
	,  4.
	,  12.
	,  16.
	,  1.
	,  8.
	,  14.
	,  0.
	,  2.
	,  3.
	,  4.
	,  8.
	,  14.
	,  24.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
	,  0.
	,  0.
	,  2.
	,  5.
	,  5.
	,  5.
	,  5.
	,  6.
	,  6.
	,  7.
	,  8.
	,  8.
	,  10.
	,  10.
	,  12.
	,  14.
	,  14.
	,  18.
	,  20.
	,  20.
	,  22.
	,  22.
	,  24.
	,  24.
	,  28.
	,  28.
	,  28.
	,  28.
	,  28.
	,  32.
	,  32.
	,  32.
	,  36.
	,  36.
	,  36.
	,  36.
	,  36.
	,  36.
	,  36.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.811384363481847
	, -0.568199310990094e4
	, -0.178657198172556e11
	,  0.795537657613427e32
	, -0.814568209346872e5
	, -0.659774567602874e8
	, -0.152861148659302e11
	, -0.560165667510446e12
	,  0.458384828593949e6
	, -0.385754000383848e14
	,  0.453735800004273e8
	,  0.939454935735563e12
	,  0.266572856432938e28
	, -0.547578313899097e10
	,  0.200725701112386e15
	,  0.185007245563239e13
	,  0.185135446828337e9
	, -0.170451090076385e12
	,  0.157890366037614e15
	, -0.202530509748774e16
	,  0.368193926183570e60
	,  0.170215539458936e18
	,  0.639234909918741e42
	, -0.821698160721956e15
	, -0.795260241872306e24
	,  0.233415869478510e18
	, -0.600079934586803e23
	,  0.594584382273384e25
	,  0.189461279349492e40
	, -0.810093428842645e46
	,  0.188813911076809e22
	,  0.111052244098768e36
	,  0.291133958602503e46
	, -0.329421923951460e22
	, -0.137570282536696e26
	,  0.181508996303902e28
	, -0.346865122768353e30
	, -0.211961148774260e38
	, -0.128617899887675e49
	,  0.479817895699239e65
}; //}}}
    constexpr auto v_s = 0.0028;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 40.;
    constexpr auto a   = 1.0;
    constexpr auto b   = 0.997;
    constexpr auto c   = 1.0;
    constexpr auto d   = 0.25;
    constexpr auto e   = 1.;
} // namespace v_3m
    namespace
v_3n
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  3.
	,  4.
	,  6.
	,  7.
	,  10.
	,  12.
	,  14.
	,  18.
	,  0.
	,  3.
	,  5.
	,  6.
	,  8.
	,  12.
	,  0.
	,  3.
	,  7.
	,  12.
	,  2.
	,  3.
	,  4.
	,  2.
	,  4.
	,  7.
	,  4.
	,  3.
	,  5.
	,  6.
	,  0.
	,  0.
	,  3.
	,  1.
	,  0.
	,  1.
	,  0.
	,  1.
	,  0.
	,  1.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -12.
	, -10.
	, -10.
	, -10.
	, -10.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -8.
	, -6.
	, -6.
	, -6.
	, -5.
	, -5.
	, -5.
	, -4.
	, -3.
	, -3.
	, -3.
	, -2.
	, -1.
	, -1.
	,  0.
	,  1.
	,  1.
	,  2.
	,  4.
	,  5.
	,  6.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.280967799943151e-38
	,  0.614869006573609e-30
	,  0.582238667048942e-27
	,  0.390628369238462e-22
	,  0.821445758255119e-20
	,  0.402137961842776e-14
	,  0.651718171878301e-12
	, -0.211773355803058e-7
	,  0.264953354380072e-2
	, -0.135031446451331e-31
	, -0.607246643970893e-23
	, -0.402352115234494e-18
	, -0.744938506925544e-16
	,  0.189917206526237e-12
	,  0.364975183508473e-5
	,  0.177274872361946e-25
	, -0.334952758812999e-18
	, -0.421537726098389e-8
	, -0.391048167929649e-1
	,  0.541276911564176e-13
	,  0.705412100773699e-11
	,  0.258585887897486e-8
	, -0.493111362030162e-10
	, -0.158649699894543e-5
	, -0.525037427886100
	,  0.220019901729615e-2
	, -0.643064132636925e-2
	,  0.629154149015048e2
	,  0.135147318617061e3
	,  0.240560808321713e-6
	, -0.890763306701305e-3
	, -0.440209599407714e4
	, -0.302807107747776e3
	,  0.159158748314599e4
	,  0.232534272709876e6
	, -0.792681207132600e6
	, -0.869871364662769e11
	,  0.354542769185671e12
	,  0.400849240129329e15
}; //}}}
    constexpr auto v_s = 0.0031;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 39.;
    constexpr auto a   = 0.976;
    constexpr auto b   = 0.997;
} // namespace v_3n
    namespace
v_3o
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  0.
	,  0.
	,  2.
	,  3.
	,  4.
	,  4.
	,  4.
	,  4.
	,  4.
	,  5.
	,  5.
	,  6.
	,  7.
	,  8.
	,  8.
	,  8.
	,  10.
	,  10.
	,  14.
	,  14.
	,  20.
	,  20.
	,  24.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      -12.
	, -4.
	, -1.
	, -1.
	, -10.
	, -12.
	, -8.
	, -5.
	, -4.
	, -1.
	, -4.
	, -3.
	, -8.
	, -12.
	, -10.
	, -8.
	, -4.
	, -12.
	, -8.
	, -12.
	, -8.
	, -12.
	, -10.
	, -12.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.128746023979718e-34
	, -0.735234770382342e-11
	,  0.289078692149150e-2
	,  0.244482731907223
	,  0.141733492030985e-23
	, -0.354533853059476e-28
	, -0.594539202901431e-17
	, -0.585188401782779e-8
	,  0.201377325411803e-5
	,  0.138647388209306e1
	, -0.173959365084772e-4
	,  0.137680878349369e-2
	,  0.814897605805513e-14
	,  0.425596631351839e-25
	, -0.387449113787755e-17
	,  0.139814747930240e-12
	, -0.171849638951521e-2
	,  0.641890529513296e-21
	,  0.118960578072018e-10
	, -0.155282762571611e-17
	,  0.233907907347507e-7
	, -0.174093247766213e-12
	,  0.377682649089149e-8
	, -0.516720236575302e-10
}; //}}}
    constexpr auto v_s = 0.0034;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 24.;
    constexpr auto a   = 0.974;
    constexpr auto b   = 0.996;
    constexpr auto c   = 0.5;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3o
    namespace
v_3p
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  0.
	,  0.
	,  0.
	,  1.
	,  2.
	,  3.
	,  3.
	,  4.
	,  6.
	,  7.
	,  7.
	,  8.
	,  10.
	,  12.
	,  12.
	,  12.
	,  14.
	,  14.
	,  14.
	,  16.
	,  18.
	,  20.
	,  22.
	,  24.
	,  24.
	,  36.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      -1.
	,  0.
	,  1.
	,  2.
	,  1.
	, -1.
	, -3.
	,  0.
	, -2.
	, -2.
	, -5.
	, -4.
	, -2.
	, -3.
	, -12.
	, -6.
	, -5.
	, -10.
	, -8.
	, -3.
	, -8.
	, -8.
	, -10.
	, -10.
	, -12.
	, -8.
	, -12.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.982825342010366e-4
	,  0.105145700850612e1
	,  0.116033094095084e3
	,  0.324664750281543e4
	, -0.123592348610137e4
	, -0.561403450013495e-1
	,  0.856677401640869e-7
	,  0.236313425393924e3
	,  0.972503292350109e-2
	, -0.103001994531927e1
	, -0.149653706199162e-8
	, -0.215743778861592e-4
	, -0.834452198291445e1
	,  0.586602660564988
	,  0.343480022104968e-25
	,  0.816256095947021e-5
	,  0.294985697916798e-2
	,  0.711730466276584e-16
	,  0.400954763806941e-9
	,  0.107766027032853e2
	, -0.409449599138182e-6
	, -0.729121307758902e-5
	,  0.677107970938909e-8
	,  0.602745973022975e-7
	, -0.382323011855257e-10
	,  0.179946628317437e-2
	, -0.345042834640005e-3
}; //}}}
    constexpr auto v_s = 0.0041;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 27.;
    constexpr auto a   = 0.972;
    constexpr auto b   = 0.997;
    constexpr auto c   = 0.5;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3p
    namespace
v_3q
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -10.
	, -10.
	, -10.
	, -10.
	, -8.
	, -6.
	, -5.
	, -5.
	, -4.
	, -4.
	, -3.
	, -2.
	, -2.
	, -2.
	, -2.
	, -1.
	, -1.
	, -1.
	,  0.
	,  1.
	,  1.
	,  1.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       10.
	,  12.
	,  6.
	,  7.
	,  8.
	,  10.
	,  8.
	,  6.
	,  2.
	,  5.
	,  3.
	,  4.
	,  3.
	,  0.
	,  1.
	,  2.
	,  4.
	,  0.
	,  1.
	,  2.
	,  0.
	,  0.
	,  1.
	,  3.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.820433843259950e5
	,  0.473271518461586e11
	, -0.805950021005413e-1
	,  0.328600025435980e2
	, -0.356617029982490e4
	, -0.172985781433335e10
	,  0.351769232729192e8
	, -0.775489259985144e6
	,  0.710346691966018e-4
	,  0.993499883820274e5
	, -0.642094171904570
	, -0.612842816820083e4
	,  0.232808472983776e3
	, -0.142808220416837e-4
	, -0.643596060678456e-2
	, -0.428577227475614e1
	,  0.225689939161918e4
	,  0.100355651721510e-2
	,  0.333491455143516
	,  0.109697576888873e1
	,  0.961917379376452
	, -0.838165632204598e-1
	,  0.247795908411492e1
	, -0.319114969006533e4
}; //}}}
    constexpr auto v_s = 0.0022;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 24.;
    constexpr auto a   = 0.848;
    constexpr auto b   = 0.983;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3q
    namespace
v_3r
{
    constexpr auto
I = array_t //{{{
{
      -8.
	, -8.
	, -3.
	, -3.
	, -3.
	, -3.
	, -3.
	,  0.
	,  0.
	,  0.
	,  0.
	,  3.
	,  3.
	,  8.
	,  8.
	,  8.
	,  8.
	,  10.
	,  10.
	,  10.
	,  10.
	,  10.
	,  10.
	,  10.
	,  10.
	,  12.
	,  14.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       6.
	,  14.
	, -3.
	,  3.
	,  4.
	,  5.
	,  8.
	, -1.
	,  0.
	,  1.
	,  5.
	, -6.
	, -2.
	, -12.
	, -10.
	, -8.
	, -5.
	, -12.
	, -10.
	, -8.
	, -6.
	, -5.
	, -4.
	, -3.
	, -2.
	, -12.
	, -12.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.144165955660863e-2
	, -0.701438599628258e13
	, -0.830946716459219e-16
	,  0.261975135368109
	,  0.393097214706245e3
	, -0.104334030654021e5
	,  0.490112654154211e9
	, -0.147104222772069e-3
	,  0.103602748043408e1
	,  0.305308890065089e1
	, -0.399745276971264e7
	,  0.569233719593750e-11
	, -0.464923504407778e-1
	, -0.535400396512906e-17
	,  0.399988795693162e-12
	, -0.536479560201811e-6
	,  0.159536722411202e-1
	,  0.270303248860217e-14
	,  0.244247453858506e-7
	, -0.983430636716454e-5
	,  0.663513144224454e-1
	, -0.993456957845006e1
	,  0.546491323528491e3
	, -0.143365406393758e5
	,  0.150764974125511e6
	, -0.337209709340105e-9
	,  0.377501980025469e-8
}; //}}}
    constexpr auto v_s = 0.0054;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 27.;
    constexpr auto a   = 0.874;
    constexpr auto b   = 0.982;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3r
    namespace
v_3s
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -10.
	, -8.
	, -6.
	, -5.
	, -5.
	, -4.
	, -4.
	, -3.
	, -3.
	, -2.
	, -1.
	, -1.
	, -1.
	,  0.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  3.
	,  3.
	,  3.
	,  4.
	,  4.
	,  4.
	,  5.
	,  14.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       20.
	,  24.
	,  22.
	,  14.
	,  36.
	,  8.
	,  16.
	,  6.
	,  32.
	,  3.
	,  8.
	,  4.
	,  1.
	,  2.
	,  3.
	,  0.
	,  1.
	,  4.
	,  28.
	,  0.
	,  32.
	,  0.
	,  1.
	,  2.
	,  3.
	,  18.
	,  24.
	,  4.
	,  24.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.532466612140254e23
	,  0.100415480000824e32
	, -0.191540001821367e30
	,  0.105618377808847e17
	,  0.202281884477061e59
	,  0.884585472596134e8
	,  0.166540181638363e23
	, -0.313563197669111e6
	, -0.185662327545324e54
	, -0.624942093918942e-1
	, -0.504160724132590e10
	,  0.187514491833092e5
	,  0.121399979993217e-2
	,  0.188317043049455e1
	, -0.167073503962060e4
	,  0.965961650599775
	,  0.294885696802488e1
	, -0.653915627346115e5
	,  0.604012200163444e50
	, -0.198339358557937
	, -0.175984090163501e58
	,  0.356314881403987e1
	, -0.575991255144384e3
	,  0.456213415338071e5
	, -0.109174044987829e8
	,  0.437796099975134e34
	, -0.616552611135792e46
	,  0.193568768917797e10
	,  0.950898170425042e54
}; //}}}
    constexpr auto v_s = 0.0022;
    constexpr auto p_s = 21e6  ;
    constexpr auto t_s = 640.  ;
    constexpr auto N   = 29.;
    constexpr auto a   = 0.886;
    constexpr auto b   = 0.990;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 4.0;
} // namespace v_3s
    namespace
v_3t
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  2.
	,  2.
	,  2.
	,  3.
	,  3.
	,  4.
	,  4.
	,  7.
	,  7.
	,  7.
	,  7.
	,  7.
	,  10.
	,  10.
	,  10.
	,  10.
	,  10.
	,  18.
	,  20.
	,  22.
	,  22.
	,  24.
	,  28.
	,  32.
	,  32.
	,  32.
	,  36.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       0.
	,  1.
	,  4.
	,  12.
	,  0.
	,  10.
	,  0.
	,  6.
	,  14.
	,  3.
	,  8.
	,  0.
	,  10.
	,  3.
	,  4.
	,  7.
	,  20.
	,  36.
	,  10.
	,  12.
	,  14.
	,  16.
	,  22.
	,  18.
	,  32.
	,  22.
	,  36.
	,  24.
	,  28.
	,  22.
	,  32.
	,  36.
	,  36.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.155287249586268e1
	,  0.664235115009031e1
	, -0.289366236727210e4
	, -0.385923202309848e13
	, -0.291002915783761e1
	, -0.829088246858083e12
	,  0.176814899675218e1
	, -0.534686695713469e9
	,  0.160464608687834e18
	,  0.196435366560186e6
	,  0.156637427541729e13
	, -0.178154560260006e1
	, -0.229746237623692e16
	,  0.385659001648006e8
	,  0.110554446790543e10
	, -0.677073830687349e14
	, -0.327910592086523e31
	, -0.341552040860644e51
	, -0.527251339709047e21
	,  0.245375640937055e24
	, -0.168776617209269e27
	,  0.358958955867578e29
	, -0.656475280339411e36
	,  0.355286045512301e39
	,  0.569021454413270e58
	, -0.700584546433113e48
	, -0.705772623326374e65
	,  0.166861176200148e53
	, -0.300475129680486e61
	, -0.668481295196808e51
	,  0.428432338620678e69
	, -0.444227367758304e72
	, -0.281396013562745e77
}; //}}}
    constexpr auto v_s = 0.0088;
    constexpr auto p_s = 20e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 33.;
    constexpr auto a   = 0.803;
    constexpr auto b   = 1.020;
    constexpr auto c   = 1.0;
    constexpr auto d   = 1.0;
    constexpr auto e   = 1.0;
} // namespace v_3t
    namespace
v_3u
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -10.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -6.
	, -6.
	, -5.
	, -5.
	, -5.
	, -3.
	, -1.
	, -1.
	, -1.
	, -1.
	,  0.
	,  0.
	,  1.
	,  2.
	,  2.
	,  3.
	,  5.
	,  5.
	,  5.
	,  6.
	,  6.
	,  8.
	,  8.
	,  10.
	,  12.
	,  12.
	,  12.
	,  14.
	,  14.
	,  14.
	,  14.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       14.
	,  10.
	,  12.
	,  14.
	,  10.
	,  12.
	,  14.
	,  8.
	,  12.
	,  4.
	,  8.
	,  12.
	,  2.
	, -1.
	,  1.
	,  12.
	,  14.
	, -3.
	,  1.
	, -2.
	,  5.
	,  10.
	, -5.
	, -4.
	,  2.
	,  3.
	, -5.
	,  2.
	, -8.
	,  8.
	, -4.
	, -12.
	, -4.
	,  4.
	, -12.
	, -10.
	, -6.
	,  6.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.122088349258355e18
	,  0.104216468608488e10
	, -0.882666931564652e16
	,  0.259929510849499e20
	,  0.222612779142211e15
	, -0.878473585050085e18
	, -0.314432577551552e22
	, -0.216934916996285e13
	,  0.159079648196849e21
	, -0.339567617303423e3
	,  0.884387651337836e13
	, -0.843405926846418e21
	,  0.114178193518022e2
	, -0.122708229235641e-3
	, -0.106201671767107e3
	,  0.903443213959313e25
	, -0.693996270370852e28
	,  0.648916718965575e-8
	,  0.718957567127851e4
	,  0.105581745346187e-2
	, -0.651903203602581e15
	, -0.160116813274676e25
	, -0.510254294237837e-8
	, -0.152355388953402
	,  0.677143292290144e12
	,  0.276378438378930e15
	,  0.116862983141686e-1
	, -0.301426947980171e14
	,  0.169719813884840e-7
	,  0.104674840020929e27
	, -0.108016904560140e5
	, -0.990623601934295e-12
	,  0.536116483602738e7
	,  0.226145963747881e22
	, -0.488731565776210e-9
	,  0.151001548880670e-4
	, -0.227700464643920e5
	, -0.781754507698846e28
}; //}}}
    constexpr auto v_s = 0.0026;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 38.;
    constexpr auto a   = 0.902;
    constexpr auto b   = 0.988;
    constexpr auto c   = 1.;
    constexpr auto d   = 1.;
    constexpr auto e   = 1.;
} // namespace v_3u
    namespace
v_3v
{
    constexpr auto
I = array_t //{{{
{
      -10.
	, -8.
	, -6.
	, -6.
	, -6.
	, -6.
	, -6.
	, -6.
	, -5.
	, -5.
	, -5.
	, -5.
	, -5.
	, -5.
	, -4.
	, -4.
	, -4.
	, -4.
	, -3.
	, -3.
	, -3.
	, -2.
	, -2.
	, -1.
	, -1.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  3.
	,  4.
	,  4.
	,  4.
	,  5.
	,  8.
	,  10.
	,  12.
	,  14.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      -8.
	, -12.
	, -12.
	, -3.
	,  5.
	,  6.
	,  8.
	,  10.
	,  1.
	,  2.
	,  6.
	,  8.
	,  10.
	,  14.
	, -12.
	, -10.
	, -6.
	,  10.
	, -3.
	,  10.
	,  12.
	,  2.
	,  4.
	, -2.
	,  0.
	, -2.
	,  6.
	,  10.
	, -12.
	, -10.
	,  3.
	, -6.
	,  3.
	,  10.
	,  2.
	, -12.
	, -2.
	, -3.
	,  1.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.415652812061591e-54
	,  0.177441742924043e-60
	, -0.357078668203377e-54
	,  0.359252213604114e-25
	, -0.259123736380269e2
	,  0.594619766193460e5
	, -0.624184007103158e11
	,  0.313080299915944e17
	,  0.105006446192036e-8
	, -0.192824336984852e-5
	,  0.654144373749937e6
	,  0.513117462865044e13
	, -0.697595750347391e19
	, -0.103977184454767e29
	,  0.119563135540666e-47
	, -0.436677034051655e-41
	,  0.926990036530639e-29
	,  0.587793105620748e21
	,  0.280375725094731e-17
	, -0.192359972440634e23
	,  0.742705723302738e27
	, -0.517429682450605e2
	,  0.820612048645469e7
	, -0.188214882341448e-8
	,  0.184587261114837e-1
	, -0.135830407782663e-5
	, -0.723681885626348e17
	, -0.223449194054124e27
	, -0.111526741826431e-34
	,  0.276032601145151e-28
	,  0.134856491567853e15
	,  0.652440293345860e-9
	,  0.510655119774360e17
	, -0.468138358908732e32
	, -0.760667491183279e16
	, -0.417247986986821e-18
	,  0.312545677756104e14
	, -0.100375333864186e15
	,  0.247761392329058e27
}; //}}}
    constexpr auto v_s = 0.0031;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 39.;
    constexpr auto a   = 0.960;
    constexpr auto b   = 0.995;
    constexpr auto c   = 1.;
    constexpr auto d   = 1.;
    constexpr auto e   = 1.;
} // namespace v_3v
    namespace
v_3w
{
    constexpr auto
I = array_t //{{{
{
      -12.
	, -12.
	, -10.
	, -10.
	, -8.
	, -8.
	, -8.
	, -6.
	, -6.
	, -6.
	, -6.
	, -5.
	, -4.
	, -4.
	, -3.
	, -3.
	, -2.
	, -2.
	, -1.
	, -1.
	, -1.
	,  0.
	,  0.
	,  1.
	,  2.
	,  2.
	,  3.
	,  3.
	,  5.
	,  5.
	,  5.
	,  8.
	,  8.
	,  10.
	,  10.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       8.
	,  14.
	, -1.
	,  8.
	,  6.
	,  8.
	,  14.
	, -4.
	, -3.
	,  2.
	,  8.
	, -10.
	, -1.
	,  3.
	, -10.
	,  3.
	,  1.
	,  2.
	, -8.
	, -4.
	,  1.
	, -12.
	,  1.
	, -1.
	, -1.
	,  2.
	, -12.
	, -5.
	, -10.
	, -8.
	, -6.
	, -12.
	, -10.
	, -12.
	, -8.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.586219133817016e-7
	, -0.894460355005526e11
	,  0.531168037519774e-30
	,  0.109892402329239
	, -0.575368389425212e-1
	,  0.228276853990249e5
	, -0.158548609655002e19
	,  0.329865748576503e-27
	, -0.634987981190669e-24
	,  0.615762068640611e-8
	, -0.961109240985747e8
	, -0.406274286652625e-44
	, -0.471103725498077e-12
	,  0.725937724828145
	,  0.187768525763682e-38
	, -0.103308436323771e4
	, -0.662552816342168e-1
	,  0.579514041765710e3
	,  0.237416732616644e-26
	,  0.271700235739893e-14
	, -0.907886213483600e2
	, -0.171242509570207e-36
	,  0.156792067854621e3
	,  0.923261357901470
	, -0.597865988422577e1
	,  0.321988767636389e7
	, -0.399441390042203e-29
	,  0.493429086046981e-7
	,  0.812036983370565e-19
	, -0.207610284654137e-11
	, -0.340821291419719e-6
	,  0.542000573372233e-17
	, -0.856711586510214e-12
	,  0.266170454405981e-13
	,  0.858133791857099e-5
}; //}}}
    constexpr auto v_s = 0.0039;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 35.;
    constexpr auto a   = 0.959;
    constexpr auto b   = 0.995;
    constexpr auto c   = 1.;
    constexpr auto d   = 1.;
    constexpr auto e   = 4.;
} // namespace v_3w
    namespace
v_3x
{
    constexpr auto
I = array_t //{{{
{
      -8.
	, -6.
	, -5.
	, -4.
	, -4.
	, -4.
	, -3.
	, -3.
	, -1.
	,  0.
	,  0.
	,  0.
	,  1.
	,  1.
	,  2.
	,  3.
	,  3.
	,  3.
	,  4.
	,  5.
	,  5.
	,  5.
	,  6.
	,  8.
	,  8.
	,  8.
	,  8.
	,  10.
	,  12.
	,  12.
	,  12.
	,  12.
	,  14.
	,  14.
	,  14.
	,  14.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       14.
	,  10.
	,  10.
	,  1.
	,  2.
	,  14.
	, -2.
	,  12.
	,  5.
	,  0.
	,  4.
	,  10.
	, -10.
	, -1.
	,  6.
	, -12.
	,  0.
	,  8.
	,  3.
	, -6.
	, -2.
	,  1.
	,  1.
	, -6.
	, -3.
	,  1.
	,  8.
	, -8.
	, -10.
	, -8.
	, -5.
	, -4.
	, -12.
	, -10.
	, -8.
	, -6.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.377373741298151e19
	, -0.507100883722913e13
	, -0.103363225598860e16
	,  0.184790814320773e-5
	, -0.924729378390945e-3
	, -0.425999562292738e24
	, -0.462307771873973e-12
	,  0.107319065855767e22
	,  0.648662492280682e11
	,  0.244200600688281e1
	, -0.851535733484258e10
	,  0.169894481433592e22
	,  0.215780222509020e-26
	, -0.320850551367334
	, -0.382642448458610e17
	, -0.275386077674421e-28
	, -0.563199253391666e6
	, -0.326068646279314e21
	,  0.397949001553184e14
	,  0.100824008584757e-6
	,  0.162234569738433e5
	, -0.432355225319745e11
	, -0.592874245598610e12
	,  0.133061647281106e1
	,  0.157338197797544e7
	,  0.258189614270853e14
	,  0.262413209706358e25
	, -0.920011937431142e-1
	,  0.220213765905426e-2
	, -0.110433759109547e2
	,  0.847004870612087e7
	, -0.592910695762536e9
	, -0.183027173269660e-4
	,  0.181339603516302
	, -0.119228759669889e4
	,  0.430867658061468e7
}; //}}}
    constexpr auto v_s = 0.0049;
    constexpr auto p_s = 23e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 36.;
    constexpr auto a   = 0.910;
    constexpr auto b   = 0.988;
    constexpr auto c   = 1.;
    constexpr auto d   = 1.;
    constexpr auto e   = 1.;
} // namespace v_3x
    namespace
v_3y
{
    constexpr auto
I = array_t //{{{
{
       0.
	,  0.
	,  0.
	,  0.
	,  1.
	,  2.
	,  2.
	,  2.
	,  2.
	,  3.
	,  3.
	,  3.
	,  4.
	,  4.
	,  5.
	,  5.
	,  8.
	,  8.
	,  10.
	,  12.
}; //}}}
    constexpr auto
J = array_t //{{{
{
      -3.
	,  1.
	,  5.
	,  8.
	,  8.
	, -4.
	, -1.
	,  4.
	,  5.
	, -8.
	,  4.
	,  8.
	, -6.
	,  6.
	, -2.
	,  1.
	, -8.
	, -2.
	, -5.
	, -8.
}; //}}}
    constexpr auto
n = array_t //{{{
{
      -0.525597995024633e-9
	,  0.583441305228407e4
	, -0.134778968457925e17
	,  0.118973500934212e26
	, -0.159096490904708e27
	, -0.315839902302021e-6
	,  0.496212197158239e3
	,  0.327777227273171e19
	, -0.527114657850696e22
	,  0.210017506281863e-16
	,  0.705106224399834e21
	, -0.266713136106469e31
	, -0.145370512554562e-7
	,  0.149333917053130e28
	, -0.149795620287641e8
	, -0.381881906271100e16
	,  0.724660165585797e-4
	, -0.937808169550193e14
	,  0.514411468376383e10
	, -0.828198594040141e5
}; //}}}
    constexpr auto v_s = 0.0031;
    constexpr auto p_s = 22e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 20.;
    constexpr auto a   = 0.996;
    constexpr auto b   = 0.994;
    constexpr auto c   = 1.;
    constexpr auto d   = 1.;
    constexpr auto e   = 4.;
} // namespace v_3y
    namespace
v_3z
{
    constexpr auto
I = array_t //{{{
{
      -8.
	, -6.
	, -5.
	, -5.
	, -4.
	, -4.
	, -4.
	, -3.
	, -3.
	, -3.
	, -2.
	, -1.
	,  0.
	,  1.
	,  2.
	,  3.
	,  3.
	,  6.
	,  6.
	,  6.
	,  6.
	,  8.
	,  8.
}; //}}}
    constexpr auto
J = array_t //{{{
{
       3.
	,  6.
	,  6.
	,  8.
	,  5.
	,  6.
	,  8.
	, -2.
	,  5.
	,  6.
	,  2.
	, -6.
	,  3.
	,  1.
	,  6.
	, -6.
	, -2.
	, -6.
	, -5.
	, -4.
	, -1.
	, -8.
	, -4.
}; //}}}
    constexpr auto
n = array_t //{{{
{
       0.244007892290650e-10
	, -0.463057430331242e7
	,  0.728803274777712e10
	,  0.327776302858856e16
	, -0.110598170118409e10
	, -0.323899915729957e13
	,  0.923814007023245e16
	,  0.842250080413712e-12
	,  0.663221436245506e12
	, -0.167170186672139e15
	,  0.253749358701391e4
	, -0.819731559610523e-20
	,  0.328380587890663e12
	, -0.625004791171543e8
	,  0.803197957462023e21
	, -0.204397011338353e-10
	, -0.378391047055938e4
	,  0.972876545938620e-2
	,  0.154355721681459e2
	, -0.373962862928643e4
	, -0.682859011374572e11
	, -0.248488015614543e-3
	,  0.394536049497068e7
}; //}}}
    constexpr auto v_s = 0.0038;
    constexpr auto p_s = 22e6  ;
    constexpr auto t_s = 650.  ;
    constexpr auto N   = 23.;
    constexpr auto a   = 0.993;
    constexpr auto b   = 0.994;
    constexpr auto c   = 1.;
    constexpr auto d   = 1.;
    constexpr auto e   = 4.;
} //namespace v_3z

    namespace
detail
{
    constexpr auto a = 1;
    constexpr auto b = 2;
    constexpr auto c = 3;
    constexpr auto d = 4;
    constexpr auto e = 5;
    constexpr auto f = 6;
    constexpr auto g = 7;
    constexpr auto h = 8;
    constexpr auto i = 9;
    constexpr auto j = 10;
    constexpr auto k = 11;
    constexpr auto l = 12;
    constexpr auto m = 13;
    constexpr auto n = 14;
    constexpr auto o = 15;
    constexpr auto p = 16;
    constexpr auto q = 17;
    constexpr auto r = 18;
    constexpr auto s = 19;
    constexpr auto t = 20;
    constexpr auto u = 21;
    constexpr auto v = 22;
    constexpr auto w = 23;
    constexpr auto x = 24;
    constexpr auto y = 25;
    constexpr auto z = 26;

        constexpr auto
    subregion_pt (auto const& pressure, auto const& temperature)
    {
        if (pressure > (40e6))
        {
            if (temperature <= t_3ab_p (pressure)) return a;
            return b;
        }
        if (pressure > (25e6))
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= t_3ab_p (pressure)) return d;
            if (temperature <= t_3ef_p (pressure)) return e;
            return f;
        }
        if (pressure > (23.5e6))
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= t_3gh_p (pressure)) return g;
            if (temperature <= t_3ef_p (pressure)) return h;
            if (temperature <= t_3ij_p (pressure)) return i;
            if (temperature <= t_3jk_p (pressure)) return j;
            return k;
        }
        if (pressure > (22.5e6))
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= t_3gh_p (pressure)) return l;
            if (temperature <= t_3mn_p (pressure)) return m;
            if (temperature <= t_3ef_p (pressure)) return n;
            if (temperature <= t_3op_p (pressure)) return o;
            if (temperature <= t_3ij_p (pressure)) return p;
            if (temperature <= t_3jk_p (pressure)) return j;
            return k;
        }
        // supercritical
        if (pressure > 22.11e6)
        {
            // sec. 5
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= t_3qu_p (pressure)) return q;
            // sec. 6
            if (temperature <= t_3uv_p (pressure)) return u;
            if (temperature <= t_3ef_p (pressure)) return v;
            if (temperature <= t_3wx_p (pressure)) return w;
            if (temperature <= t_3rx_p (pressure)) return x;
            // sec. 5
            if (temperature <= t_3jk_p (pressure)) return r;
            return k;
        }
        if (pressure > 22.064e6)
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= t_3qu_p (pressure)) return q;
            if (temperature <= t_3uv_p (pressure)) return u;
            if (temperature <= t_3ef_p (pressure)) return y;
            if (temperature <= t_3wx_p (pressure)) return z;
            if (temperature <= t_3rx_p (pressure)) return x;
            if (temperature <= t_3jk_p (pressure)) return r;
            return k;
        }
        // subcritical
        if (pressure > 21.93161551e6)
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= t_3qu_p (pressure)) return q;
            if (temperature <= t_3uv_p (pressure)) return u;
            if (temperature <= saturation_temperature_p (pressure)) return y;
            if (temperature <= t_3wx_p (pressure)) return z;
            if (temperature <= t_3rx_p (pressure)) return x;
            if (temperature <= t_3jk_p (pressure)) return r;
            return k;
        }
        if (pressure > 21.93161551e6)
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= t_3qu_p (pressure)) return q;
            if (temperature <= saturation_temperature_p (pressure)) return u;
            if (temperature <= t_3rx_p (pressure)) return x;
            if (temperature <= t_3jk_p (pressure)) return r;
            return k;
        }
        if (pressure > 21.04336732e6)
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= t_3qu_p (pressure)) return q;
            if (temperature <= saturation_temperature_p (pressure)) return u;
            if (temperature <= t_3rx_p (pressure)) return x;
            if (temperature <= t_3jk_p (pressure)) return r;
            return k;
        }
        if (pressure > 20.5e6)
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= saturation_temperature_p (pressure)) return s;
            if (temperature <= t_3jk_p (pressure)) return r;
            return k;
        }
        if (pressure > 1.900881189173929e7)
        {
            if (temperature <= t_3cd_p (pressure)) return c;
            if (temperature <= saturation_temperature_p (pressure)) return s;
            return t;
        }
        if (temperature <= saturation_temperature_p (pressure)) return c;
        return t;
    }
} // namespace detail
  //
#define ISTO_IAPWS_GEN_DERIVED                                                                        \
    constexpr auto                                                                                    \
density_pt (auto const& pressure, auto const& temperature)                                            \
{                                                                                                     \
    return 1 / massic_volume_pt (pressure, temperature);                                              \
}                                                                                                     \
    constexpr auto                                                                                    \
massic_internal_energy_pt (auto const& pressure, auto const& temperature)                             \
{                                                                                                     \
    return massic_internal_energy_dt (density_pt (pressure, temperature), temperature);               \
}                                                                                                     \
    constexpr auto                                                                                    \
massic_entropy_pt (auto const& pressure, auto const& temperature)                                     \
{                                                                                                     \
    return massic_entropy_dt (density_pt (pressure, temperature), temperature);                       \
}                                                                                                     \
    constexpr auto                                                                                    \
massic_enthalpy_pt (auto const& pressure, auto const& temperature)                                    \
{                                                                                                     \
    return massic_enthalpy_dt (density_pt (pressure, temperature), temperature);                      \
}                                                                                                     \
    constexpr auto                                                                                    \
massic_isochoric_heat_capacity_pt (auto const& pressure, auto const& temperature)                     \
{                                                                                                     \
    return massic_isochoric_heat_capacity_dt (density_pt (pressure, temperature), temperature);       \
}                                                                                                     \
    constexpr auto                                                                                    \
massic_isobaric_heat_capacity_pt (auto const& pressure, auto const& temperature)                      \
{                                                                                                     \
    return massic_isobaric_heat_capacity_dt (density_pt (pressure, temperature), temperature);        \
}                                                                                                     \
    constexpr auto                                                                                    \
speed_of_sound_pt (auto const& pressure, auto const& temperature)                                     \
{                                                                                                     \
    return speed_of_sound_dt (density_pt (pressure, temperature), temperature);                       \
}                                                                                                     \
    constexpr auto                                                                                    \
isobaric_cubic_expansion_coefficient_pt (auto const& pressure, auto const& temperature)               \
{                                                                                                     \
    return isobaric_cubic_expansion_coefficient_dt (density_pt (pressure, temperature), temperature); \
}                                                                                                     \
    constexpr auto                                                                                    \
isothermal_compressibility_pt (auto const& pressure, auto const& temperature)                         \
{                                                                                                     \
    return isothermal_compressibility_dt (density_pt (pressure, temperature), temperature);           \
}                                                                                                     \
    constexpr auto                                                                                    \
relative_pressure_coefficient_pt (auto const& pressure, auto const& temperature)                      \
{                                                                                                     \
    return relative_pressure_coefficient_dt (density_pt (pressure, temperature), temperature);        \
}                                                                                                     \
    constexpr auto                                                                                    \
isothermal_stress_coefficient_pt (auto const& pressure, auto const& temperature)                      \
{                                                                                                     \
    return isothermal_stress_coefficient_dt (density_pt (pressure, temperature), temperature);        \
}                                                                                                     \
    constexpr auto                                                                                    \
massic_helmholtz_free_energy_pt (auto const& pressure, auto const& temperature)                       \
{                                                                                                     \
    return massic_helmholtz_free_energy_dt (density_pt (pressure, temperature), temperature);         \
}                                                                                                     \
    constexpr auto                                                                                    \
massic_gibbs_free_energy_pt (auto const& pressure, auto const& temperature)                           \
{                                                                                                     \
    return massic_gibbs_free_energy_dt (density_pt (pressure, temperature), temperature);             \
}

#define ISTO_IAPWS_GEN_R3_FUNC(X)     \
    namespace                         \
r##X                                  \
{                                     \
    constexpr auto                    \
massic_volume_pt (                    \
      auto const& pressure            \
    , auto const& temperature         \
){                                    \
        using namespace v_3##X;       \
        using std::pow;               \
        auto const                    \
    pi = pressure / p_s;              \
        auto const                    \
    theta = temperature / t_s;        \
    return pow (sum (                 \
          n                           \
        * pow (pow (pi - a, c), I)    \
        * pow (pow (theta - b, d), J) \
    ), e) * v_s;                      \
}                                     \
ISTO_IAPWS_GEN_DERIVED                \
} // namespace X
ISTO_IAPWS_GEN_R3_FUNC(a)
ISTO_IAPWS_GEN_R3_FUNC(b)
ISTO_IAPWS_GEN_R3_FUNC(c)
ISTO_IAPWS_GEN_R3_FUNC(d)
ISTO_IAPWS_GEN_R3_FUNC(e)
ISTO_IAPWS_GEN_R3_FUNC(f)
ISTO_IAPWS_GEN_R3_FUNC(g)
ISTO_IAPWS_GEN_R3_FUNC(h)
ISTO_IAPWS_GEN_R3_FUNC(i)
ISTO_IAPWS_GEN_R3_FUNC(j)
ISTO_IAPWS_GEN_R3_FUNC(k)
ISTO_IAPWS_GEN_R3_FUNC(l)
ISTO_IAPWS_GEN_R3_FUNC(m)
ISTO_IAPWS_GEN_R3_FUNC(o)
ISTO_IAPWS_GEN_R3_FUNC(p)
ISTO_IAPWS_GEN_R3_FUNC(q)
ISTO_IAPWS_GEN_R3_FUNC(r)
ISTO_IAPWS_GEN_R3_FUNC(s)
ISTO_IAPWS_GEN_R3_FUNC(t)
ISTO_IAPWS_GEN_R3_FUNC(u)
ISTO_IAPWS_GEN_R3_FUNC(v)
ISTO_IAPWS_GEN_R3_FUNC(w)
ISTO_IAPWS_GEN_R3_FUNC(x)
ISTO_IAPWS_GEN_R3_FUNC(y)
ISTO_IAPWS_GEN_R3_FUNC(z)
#undef ISTO_IAPWS_GEN_R3_FUNC
    namespace
rn
{
    constexpr auto
massic_volume_pt (
      auto const& pressure
    , auto const& temperature
){
        using namespace v_3n;
        using std::exp, std::pow;
        auto const
    pi = pressure / p_s;
        auto const
    theta = temperature / t_s;
    return exp (sum (n * pow (pi - a, I) * pow (theta - b, J))) * v_s;
}
} // namespace n

    constexpr auto
massic_volume_pt (auto const& pressure, auto const& temperature)
{
        using std::pow;
    switch (detail::subregion_pt (pressure, temperature))
    {
#define ISTO_IAPWS_GEN_V_CASE(X)              \
        case detail::X:                       \
        {                                     \
                using namespace v_3##X;       \
                auto const                    \
            pi = pressure / p_s;              \
                auto const                    \
            theta = temperature / t_s;        \
            return pow (sum (                 \
                  n                           \
                * pow (pow (pi - a, c), I)    \
                * pow (pow (theta - b, d), J) \
            ), e) * v_s;                      \
        }
        ISTO_IAPWS_GEN_V_CASE(a)
        ISTO_IAPWS_GEN_V_CASE(b)
        ISTO_IAPWS_GEN_V_CASE(c)
        ISTO_IAPWS_GEN_V_CASE(d)
        ISTO_IAPWS_GEN_V_CASE(e)
        ISTO_IAPWS_GEN_V_CASE(f)
        ISTO_IAPWS_GEN_V_CASE(g)
        ISTO_IAPWS_GEN_V_CASE(h)
        ISTO_IAPWS_GEN_V_CASE(i)
        ISTO_IAPWS_GEN_V_CASE(j)
        ISTO_IAPWS_GEN_V_CASE(k)
        ISTO_IAPWS_GEN_V_CASE(l)
        ISTO_IAPWS_GEN_V_CASE(m)
        case detail::n:
        {
                using namespace v_3n;
                using std::exp;
                auto const
            pi = pressure / p_s;
                auto const
            theta = temperature / t_s;
            return exp (sum (n * pow (pi - a, I) * pow (theta - b, J))) * v_s;
        }
        ISTO_IAPWS_GEN_V_CASE(o)
        ISTO_IAPWS_GEN_V_CASE(p)
        ISTO_IAPWS_GEN_V_CASE(q)
        ISTO_IAPWS_GEN_V_CASE(r)
        ISTO_IAPWS_GEN_V_CASE(s)
        ISTO_IAPWS_GEN_V_CASE(t)
        ISTO_IAPWS_GEN_V_CASE(u)
        ISTO_IAPWS_GEN_V_CASE(v)
        ISTO_IAPWS_GEN_V_CASE(w)
        ISTO_IAPWS_GEN_V_CASE(x)
        ISTO_IAPWS_GEN_V_CASE(y)
        ISTO_IAPWS_GEN_V_CASE(z)
#undef ISTO_IAPWS_GEN_V_CASE
    }
    throw internal_error_e {};
}

ISTO_IAPWS_GEN_DERIVED
#undef ISTO_IAPWS_GEN_ALL
#undef ISTO_IAPWS_GEN_DERIVED
} // namespace r3
    namespace
r5
{
    namespace
detail
{
    constexpr auto
J_0 = array_t //{{{
{
       0.
    ,  1.
    , -3.
    , -2.
    , -1.
    ,  2.
}; //}}}
    constexpr auto
n_0 = array_t //{{{
{
      -0.13179983674201e2
    ,  0.68540841634434e1
    , -0.24805148933466e-1
    ,  0.36901534980333
    , -0.31161318213925e1
    , -0.32961626538917
}; //}}}
    constexpr auto
I_r = array_t //{{{
{
      1.
    , 1.
    , 1.
    , 2.
    , 2.
    , 3.
}; //}}}
    constexpr auto
J_r = array_t //{{{
{
      1.
    , 2.
    , 3.
    , 3.
    , 9.
    , 7.
}; //}}}
    constexpr auto
n_r = array_t //{{{
{
       0.15736404855259e-2
    ,  0.90153761673944e-3
    , -0.50270077677648e-2
    ,  0.22440037409485e-5
    , -0.41163275453471e-5
    ,  0.37919454822955e-7
}; //}}}

// low-level dimensionless {{{
    constexpr auto
gamma_0 (auto const& pi, auto const& tau)
{
        using std::log;
    return log (pi) + sum (n_0 * pow (tau, J_0));
}
    constexpr auto
gamma_0_p (auto const& pi, auto const& /* tau */)
{
    return 1 / pi;
}
    constexpr auto
gamma_0_pp (auto const& pi, auto const& /* tau */)
{
    return - 1 / pi / pi;
}
    constexpr auto
gamma_0_t (auto const& /* pi */, auto const& tau)
{
    return sum (n_0 * J_0 * pow (tau, J_0 - 1.));
}
    constexpr auto
gamma_0_tt (auto const& /* pi */, auto const& tau)
{
    return sum (n_0 * J_0 * (J_0 - 1.) * pow (tau, J_0 - 2.));
}
    constexpr auto
gamma_r (auto const& pi, auto const& tau)
{
    return sum (n_r * pow (pi, I_r) * pow (tau, J_r));
}
    constexpr auto
gamma_r_p (auto const& pi, auto const& tau)
{
    return sum (n_r * I_r * pow (pi, I_r - 1.) * pow (tau, J_r));
}
    constexpr auto
gamma_r_pp (auto const& pi, auto const& tau)
{
    return sum (n_r * I_r * (I_r - 1.) * pow (pi, I_r - 2.) * pow (tau, J_r));
}
    constexpr auto
gamma_r_t (auto const& pi, auto const& tau)
{
    return sum (n_r * pow (pi, I_r) * J_r * pow (tau, J_r - 1.));
}
    constexpr auto
gamma_r_tt (auto const& pi, auto const& tau)
{
    return sum (n_r * pow (pi, I_r) * J_r * (J_r - 1.) * pow (tau, J_r - 2.));
}
    constexpr auto
gamma_r_pt (auto const& pi, auto const& tau)
{
    return sum (n_r * I_r * pow (pi, I_r - 1.) * J_r * pow (tau, J_r - 1.));
}
    constexpr auto
gamma (auto const& pi, auto const& tau)
{
    return gamma_0 (pi, tau) + gamma_r (pi, tau);
}
    constexpr auto
gamma_p (auto const& pi, auto const& tau)
{
    return gamma_0_p (pi, tau) + gamma_r_p (pi, tau);
}
    constexpr auto
gamma_t (auto const& pi, auto const& tau)
{
    return gamma_0_t (pi, tau) + gamma_r_t (pi, tau);
}
    constexpr auto
gamma_pp (auto const& pi, auto const& tau)
{
    return gamma_0_pp (pi, tau) + gamma_r_pp (pi, tau);
}
    constexpr auto
gamma_tt (auto const& pi, auto const& tau)
{
    return gamma_0_tt (pi, tau) + gamma_r_tt (pi, tau);
}
    constexpr auto
gamma_pt (auto const& pi, auto const& tau)
{
    return /*gamma_0_pt (pi, tau) +*/ gamma_r_pt (pi, tau);
}
// }}}
} // namespace detail

    using namespace detail;

    constexpr auto
massic_volume_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 1000. * (1) / temperature;
    return pi * gamma_p (pi, tau) * massic_gas_constant * temperature / pressure;
}
    constexpr auto
massic_enthalpy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 1000. * (1) / temperature;
    return tau * gamma_t (pi, tau) * massic_gas_constant * temperature;
}
    constexpr auto
massic_internal_energy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 1000. * (1) / temperature;
    return (tau * gamma_t (pi, tau) - pi * gamma_p (pi, tau)) * massic_gas_constant * temperature;
}
    constexpr auto
massic_entropy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 1000. * (1) / temperature;
    return (tau * gamma_t (pi, tau) - gamma (pi, tau)) * massic_gas_constant;
}
    constexpr auto
massic_isobaric_heat_capacity_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 1000. * (1) / temperature;
    return - tau * tau * gamma_tt (pi, tau) * massic_gas_constant;
}
    constexpr auto
massic_isochoric_heat_capacity_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 1000. * (1) / temperature;
    return (- tau * tau * gamma_tt (pi, tau) + pow <2> (gamma_p (pi, tau) - tau * gamma_pt (pi, tau)) / gamma_pp (pi, tau)) * massic_gas_constant;
}
    constexpr auto
speed_of_sound_pt (auto const& pressure, auto const& temperature)
{
        using std::pow;
        using std::sqrt;
        auto const
    pi = pressure / 1e6;
        auto const
    tau = 1000. * (1) / temperature;
    return sqrt (
          pow (1. + pi * gamma_r_p (pi, tau), 2)
        / (
              (1. - pi * pi * gamma_r_pp (pi, tau))
            + pow (1. + pi * gamma_r_p (pi, tau) - tau * pi * gamma_r_pt (pi, tau), 2)
            / (tau * tau * gamma_tt (pi, tau))
          )
        * massic_gas_constant * temperature
    );
}
    constexpr auto
density_pt (auto const& pressure, auto const& temperature)
{
    return 1 / massic_volume_pt (pressure, temperature);
}
    constexpr auto
isobaric_cubic_expansion_coefficient_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / 1e6;
        auto const
    tau = (100.) / temperature;
    return ((1. + pi * gamma_r_p (pi, tau) - tau * pi * gamma_r_pt (pi, tau)) / (1. + pi * gamma_r_p (pi, tau))) / temperature;
}
    constexpr auto
isothermal_compressibility_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = 1386. * (1) / temperature;
    return ((1. - pi * pi * gamma_r_pp (pi, tau)) / (1. + pi * gamma_r_p (pi, tau))) / pressure;
}
    constexpr auto
relative_pressure_coefficient_pt (auto const& pressure, auto const& temperature)
{
    return isobaric_cubic_expansion_coefficient_pt (pressure, temperature) / pressure / isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
isothermal_stress_coefficient_pt (auto const& pressure, auto const& temperature)
{
    return 1. / pressure / massic_volume_pt (pressure, temperature) / isothermal_compressibility_pt (pressure, temperature);
}
    constexpr auto
massic_gibbs_free_energy_pt (auto const& pressure, auto const& temperature)
{
        auto const
    pi = pressure / (16.53 * 1e6);
        auto const
    tau = (1386. * (1)) / temperature;
    return gamma (pi, tau) * massic_gas_constant * temperature;
}
    constexpr auto
massic_helmholtz_free_energy_pt (auto const& pressure, auto const& temperature)
{
    return massic_internal_energy_pt (pressure, temperature) 
        - temperature * massic_entropy_pt(pressure, temperature);
}
} // namespace r5

// §4 Auxiliary Equation for the Boundary between Regions 2 and 3.
    constexpr auto
p_b23_t (auto const& temperature)
{
        constexpr double
    n1 =  0.34805185628969e3;
        constexpr double
    n2 = -0.11671859879975e1;
        constexpr double
    n3 =  0.10192970039326e-2;
        auto const
    t = temperature / (1);
    return (n1 + n2 * t + n3 * t * t) * 1e6;
}

    /*constexpr*/ auto // Someday, std::sqrt might be constexpr.
t_b23_p (auto const& pressure)
{
        using std::sqrt;
        constexpr double
    n3 = 0.10192970039326e-2;
        constexpr double
    n4 = 0.57254459862746e3;
        constexpr double
    n5 = 0.13918839778870e2;
        auto const
    p = pressure / 1e6;
    return (n4 + sqrt ((p - n5) / n3)) * (1);
}
    constexpr int
region_pt (auto const& pressure, auto const& temperature)
{
    if (temperature <= (623.15))
    {
            auto const
        ps = saturation_pressure_t (temperature);
        if (pressure > ps) return 1;
        if (pressure < ps) return 2;
        return 4;
    }
    if (temperature <= 1073.15)
    {
            auto const
        p23 = p_b23_t (temperature);
        if (pressure <= p23) return 2;
        return 3;
    }
    return 5;
}

// Region (h, s) (Supp-phs3-2014.pdf, sec. 4)
    namespace
detail
{
    namespace
h_p_1
{
    constexpr auto
I = array_t //{{{
{
        0.
    ,   0.
    ,   1.
    ,   1.
    ,   2.
    ,   2.
    ,   3.
    ,   3.
    ,   4.
    ,   4.
    ,   4.
    ,   5.
    ,   5.
    ,   7.
    ,   8.
    ,  12.
    ,  12.
    ,  14.
    ,  14.
    ,  16.
    ,  20.
    ,  20.
    ,  22.
    ,  24.
    ,  28.
    ,  32.
    ,  32.
}; //}}}

    constexpr auto
J = array_t //{{{
{
      14.
    , 36.
    ,  3.
    , 16.
    ,  0.
    ,  5.
    ,  4.
    , 36.
    ,  4.
    , 16.
    , 24.
    , 18.
    , 24.
    ,  1.
    ,  4.
    ,  2.
    ,  4.
    ,  1.
    , 22.
    , 10.
    , 12.
    , 28.
    ,  8.
    ,  3.
    ,  0.
    ,  6.
    ,  8.
}; //}}}

    constexpr auto
n = array_t //{{{
{
       0.332171191705237
    ,  0.611217706323496e-3
    , -0.882092478906822e1
    , -0.455628192543250
    , -0.263483840850452e-4
    , -0.223949661148062e2
    , -0.428398660164013e1
    , -0.616679338856916
    , -0.146823031104040e2
    ,  0.284523138727299e3
    , -0.113398503195444e3
    ,  0.115671380760859e4
    ,  0.395551267359325e3
    , -0.154891257229285e1
    ,  0.194486637751291e2
    , -0.357915139457043e1
    , -0.335369414148819e1
    , -0.664426796332460
    ,  0.323321885383934e5
    ,  0.331766744667084e4
    , -0.223501257931087e5
    ,  0.573953875852936e7
    ,  0.173226193407919e3
    , -0.363968822121321e-1
    ,  0.834596332878346e-6
    ,  0.503611916682674e1
    ,  0.655444787064505e2
}; //}}}
} // namespace h_p_1
    constexpr auto
h_p_1_s (auto const& massic_entropy)
{
        auto const
    sigma = massic_entropy / (3.8e3);
        using namespace h_p_1;
    return sum (n * pow (sigma - 1.09, I) * pow (sigma + 0.366e-4, J)) * (1700e3);
}
    namespace
h_p_3a
{
    constexpr auto
I = array_t //{{{
{
        0.
    ,   0.
    ,   0.
    ,   0.
    ,   2.
    ,   3.
    ,   4.
    ,   4.
    ,   5.
    ,   5.
    ,   6.
    ,   7.
    ,   7.
    ,   7.
    ,  10.
    ,  10.
    ,  10.
    ,  32.
    ,  32.
}; //}}}

    constexpr auto
J = array_t //{{{
{
       1.
    ,  4.
    , 10.
    , 16.
    ,  1.
    , 36.
    ,  3.
    , 16.
    , 20.
    , 36.
    ,  4.
    ,  2.
    , 28.
    , 32.
    , 14.
    , 32.
    , 36.
    ,  0.
    ,  6.
}; //}}}

    constexpr auto
n = array_t //{{{
{
       0.822673364673336
    ,  0.181977213534479
    , -0.112000260313624e-1
    , -0.746778287048033e-3
    , -0.179046263257381
    ,  0.424220110836657e-1
    , -0.341355823438768
    , -0.209881740853565e1
    , -0.822477343323596e1
    , -0.499684082076008e1
    ,  0.191413958471069
    ,  0.581062241093136e-1
    , -0.165505498701029e4
    ,  0.158870443421201e4
    , -0.850623535172818e2
    , -0.317714386511207e5
    , -0.945890406632871e5
    , -0.139273847088690e-5
    ,  0.631052532240980
}; //}}}
} // namespace h_p_3a
    constexpr auto
h_p_3a_s (auto const& massic_entropy)
{
        auto const
    sigma = massic_entropy / (3.8e3);
        using namespace h_p_3a;
    return sum (n * pow (sigma - 1.09, I) * pow (sigma + 0.366e-4, J)) * (1700e3);
}
    namespace
h_pp_2ab
{
    constexpr auto
I = array_t //{{{
{
        1.
    ,   1.
    ,   2.
    ,   2.
    ,   4.
    ,   4.
    ,   7.
    ,   8.
    ,   8.
    ,  10.
    ,  12.
    ,  12.
    ,  18.
    ,  20.
    ,  24.
    ,  28.
    ,  28.
    ,  28.
    ,  28.
    ,  28.
    ,  32.
    ,  32.
    ,  32.
    ,  32.
    ,  32.
    ,  36.
    ,  36.
    ,  36.
    ,  36.
    ,  36.
}; //}}}

    constexpr auto
J = array_t //{{{
{
       8.
    , 24.
    ,  4.
    , 32.
    ,  1.
    ,  2.
    ,  7.
    ,  5.
    , 12.
    ,  1.
    ,  0.
    ,  7.
    , 10.
    , 12.
    , 32.
    ,  8.
    , 12.
    , 20.
    , 22.
    , 24.
    ,  2.
    ,  7.
    , 12.
    , 14.
    , 24.
    , 10.
    , 12.
    , 20.
    , 22.
    , 28.
}; //}}}

    constexpr auto
n = array_t //{{{
{
      -0.524581170928788e3
    , -0.926947218142218e7
    , -0.237385107491666e3
    ,  0.210770155812776e11
    , -0.239494562010986e2
    ,  0.221802480294197e3
    , -0.510472533393438e7
    ,  0.124981396109147e7
    ,  0.200008436996201e10
    , -0.815158509791035e3
    , -0.157612685637523e3
    , -0.114200422332791e11
    ,  0.662364680776872e16
    , -0.227622818296144e19
    , -0.171048081348406e32
    ,  0.660788766938091e16
    ,  0.166320055886021e23
    , -0.218003784381501e30
    , -0.787276140295618e30
    ,  0.151062329700346e32
    ,  0.795732170300541e7
    ,  0.131957647355347e16
    , -0.325097068299140e24
    , -0.418600611419248e26
    ,  0.297478906557467e35
    , -0.953588761745473e20
    ,  0.166957699620939e25
    , -0.175407764869978e33
    ,  0.347581490626396e35
    , -0.710971318427851e39
}; //}}}
} // namespace h_pp_2ab
    constexpr auto
h_pp_2ab_s (auto const& massic_entropy)
{
        using namespace h_pp_2ab;
        auto const
    sigma1 = massic_entropy / (5.21e3);
        auto const
    sigma2 = massic_entropy / (9.2e3);
    return exp(sum (n * pow ((1 / sigma1) - 0.513, I) * pow (sigma2 - 0.524, J))) * (2800e3);
}
    namespace
h_pp_2c3b
{

    constexpr auto
I = array_t //{{{
{
        0.
    ,   0.
    ,   0.
    ,   1.
    ,   1.
    ,   5.
    ,   6.
    ,   7.
    ,   8.
    ,   8.
    ,  12.
    ,  16.
    ,  22.
    ,  22.
    ,  24.
    ,  36.
}; //}}}

    constexpr auto
J = array_t //{{{
{
       0.
    ,  3.
    ,  4.
    ,  0.
    , 12.
    , 36.
    , 12.
    , 16.
    ,  2.
    , 20.
    , 32.
    , 36.
    ,  2.
    , 32.
    ,  7.
    , 20.
}; //}}}

    constexpr auto
n = array_t //{{{
{
       0.104351280732769e1
    , -0.227807912708513e1
    ,  0.180535256723202e1
    ,  0.420440834792042
    , -0.105721244834660e6
    ,  0.436911607493884e25
    , -0.328032702839753e12
    , -0.678686760804270e16
    ,  0.743957464645363e4
    , -0.356896445355761e20
    ,  0.167590585186801e32
    , -0.355028625419105e38
    ,  0.396611982166538e12
    , -0.414716268484468e41
    ,  0.359080103867382e19
    , -0.116994334851995e41
}; //}}}
} // namespace h_pp_2c3b
    constexpr auto
h_pp_2c3b_s (auto const& massic_entropy)
{
        using namespace h_pp_2c3b;
        auto const
    sigma = massic_entropy / (5.9e3);
    return pow (sum (n * pow (sigma - 1.02, I) * pow (sigma - 0.726, J)), 4) * (2800e3);
}
    namespace
h_b13
{
    constexpr auto
I = array_t //{{{
{
        0.
    ,   1.
    ,   1.
    ,   3.
    ,   5.
    ,   6.
}; //}}}

    constexpr auto
J = array_t //{{{
{
        0.
    ,  -2.
    ,   2.
    , -12.
    ,  -4.
    ,  -3.
}; //}}}

    constexpr auto
n = array_t //{{{
{
       0.913965547600543
    , -0.430944856041991e-4
    ,  0.603235694765419e2
    ,  0.117518273082168e-17
    ,  0.220000904781292
    , -0.690815545851641e2
}; //}}}
} // namespace h_b13
    constexpr auto
h_b13_s (auto const& massic_entropy)
{
        using namespace h_b13;
        auto const
    sigma = massic_entropy / (3.8e3);
    return sum (n * pow (sigma - 0.884, I) * pow (sigma - 0.864, J)) * (1700e3);
}
    namespace
t_b23
{
    constexpr auto
I = array_t //{{{
{
      -12.
    , -10.
    ,  -8.
    ,  -4.
    ,  -3.
    ,  -2.
    ,  -2.
    ,  -2.
    ,  -2.
    ,   0.
    ,   1.
    ,   1.
    ,   1.
    ,   3.
    ,   3.
    ,   5.
    ,   6.
    ,   6.
    ,   8.
    ,   8.
    ,   8.
    ,  12.
    ,  12.
    ,  14.
    ,  14.
}; //}}}

    constexpr auto
J = array_t //{{{
{
       10.
    ,   8.
    ,   3.
    ,   4.
    ,   3.
    ,  -6.
    ,   2.
    ,   3.
    ,   4.
    ,   0.
    ,  -3.
    ,  -2.
    ,  10.
    ,  -2.
    ,  -1.
    ,  -5.
    ,  -6.
    ,  -3.
    ,  -8.
    ,  -2.
    ,  -1.
    , -12.
    ,  -1.
    , -12.
    ,   1.
}; //}}}

    constexpr auto
n = array_t //{{{
{
       0.629096260829810e-3
    , -0.823453502583165e-3
    ,  0.515446951519474e-7
    , -0.117565945784945e1
    ,  0.348519684726192e1
    , -0.507837382408313e-11
    , -0.284637670005479e1
    , -0.236092263939673e1
    ,  0.601492324973779e1
    ,  0.148039650824546e1
    ,  0.360075182221907e-3
    , -0.126700045009952e-1
    , -0.122184332521413e7
    ,  0.149276502463272
    ,  0.698733471798484
    , -0.252207040114321e-1
    ,  0.147151930985213e-1
    , -0.108618917681849e1
    , -0.936875039816322e-3
    ,  0.819877897570217e2
    , -0.182041861521835e3
    ,  0.261907376402688e-5
    , -0.291626417025961e5
    ,  0.140660774926165e-4
    ,  0.783237062349385e7
}; //}}}
} // namespace t_b23
    constexpr auto
t_b23_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
        using namespace t_b23;
        auto const
    eta = massic_enthalpy / (3000e3);
        auto const
    sigma = massic_entropy / (5.3e3);
    return sum (n * pow (eta - 0.727, I) * pow (sigma - 0.864, J)) * (900);
}
} // namespace detail
    namespace
r2
{
    constexpr auto
pressure_c_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    constexpr auto
I = array_t ///{{{
{
       0.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  0.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  1.
    ,  2.
    ,  2.
    ,  2.
    ,  2.
    ,  2.
    ,  3.
    ,  3.
    ,  3.
    ,  3.
    ,  3.
    ,  4.
    ,  5.
    ,  5.
    ,  5.
    ,  5.
    ,  6.
    ,  6.
    , 10.
    , 12.
    , 16.
}; //}}}

    constexpr auto
J = array_t ///{{{
{
       0.
    ,  1.
    ,  2.
    ,  3.
    ,  4.
    ,  8.
    ,  0.
    ,  2.
    ,  5.
    ,  8.
    , 14.
    ,  2.
    ,  3.
    ,  7.
    , 10.
    , 18.
    ,  0.
    ,  5.
    ,  8.
    , 16.
    , 18.
    , 18.
    ,  1.
    ,  4.
    ,  6.
    , 14.
    ,  8.
    , 18.
    ,  7.
    ,  7.
    , 10.
}; //}}}

    constexpr auto
n = array_t ///{{{
{
       0.112225607199012
    , -0.339005953606712e1
    , -0.320503911730094e2
    , -0.197597305104900e3
    , -0.407693861553446e3
    ,  0.132943775222331e5
    ,  0.170846839774007e1
    ,  0.373694198142245e2
    ,  0.358144365815434e4
    ,  0.423014446424664e6
    , -0.751071025760063e9
    ,  0.523446127607898e2
    , -0.228351290812417e3
    , -0.960652417056937e6
    , -0.807059292526074e8
    ,  0.162698017225669e13
    ,  0.772465073604171
    ,  0.463929973837746e5
    , -0.137317885134128e8
    ,  0.170470392630512e13
    , -0.251104628187308e14
    ,  0.317748830835520e14
    ,  0.538685623675312e2
    , -0.553089094625169e5
    , -0.102861522421405e7
    ,  0.204249418756234e13
    ,  0.273918446626977e9
    , -0.263963146312685e16
    , -0.107890854108088e10
    , -0.296492620980124e11
    , -0.111754907323424e16
}; //}}}
        auto const
    eta = massic_enthalpy / (3500e3);
        auto const
    sigma = massic_entropy / (5.9e3);
    return pow (sum (n * pow (eta - 0.7, I) * pow (sigma - 1.1, J)), 4) * (100e6);
}
} // namespace r2;

    constexpr int
region_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    if (massic_entropy <= (3.778281340e3))
    {
        if (massic_enthalpy <= detail::h_p_1_s (massic_entropy)) return 4;
        if (massic_enthalpy <= detail::h_b13_s (massic_entropy)) return 1;
        return 3;
    }
    if (massic_entropy <= (4.412021482234e3))
    {
        if (massic_enthalpy <= detail::h_p_3a_s (massic_entropy)) return 4;
        return 3;
    }
    if (massic_entropy <= (5.85e3))
    {
        if (massic_enthalpy <= detail::h_pp_2c3b_s (massic_entropy)) return 4;
        if (
               massic_entropy  <= (5.048096828e3)  
            || massic_enthalpy <= (2.563592004e6)
        ){
            return 3;
        }
        if (
               massic_entropy  >= (5.260578707e3)
            || massic_enthalpy >= (2.812942061e6)
        ){
            return 2;
        }
            auto const
        t = detail::t_b23_hs (massic_enthalpy, massic_entropy);
            auto const
        p = r2::pressure_c_hs (massic_enthalpy, massic_entropy);
            auto const
        p23 = p_b23_t (t);
        if (p > p23) return 3;
        return 2;
        
    }
    if (massic_enthalpy <= detail::h_pp_2ab_s (massic_entropy)) return 4;
    return 2;
}

// Region (p, h) (Supp-phs3-2014.pdf, sec. 4)   
    namespace
detail
{
        const auto
    p_s_623_15 = saturation_pressure_t (623.15);
} // namespace detail
    constexpr int
region_ph (auto const& pressure, auto const& massic_enthalpy)
{
    if (pressure <= detail::p_s_623_15)
    {
        if (massic_enthalpy <= r1::massic_enthalpy_pt (pressure, saturation_temperature_p (pressure))) return 1;
        if (massic_enthalpy <= r2::massic_enthalpy_pt (pressure, saturation_temperature_p (pressure))) return 4;
        return 2;
    }
    if (massic_enthalpy <= r1::massic_enthalpy_pt (pressure, 623.15)) return 1;
    if (massic_enthalpy > r2::massic_enthalpy_pt (pressure, t_b23_p (pressure))) return 2;
    if (pressure <= saturation_pressure_h (massic_enthalpy)) return 4;
    return 3;
}
    constexpr int
region_ps (auto const& pressure, auto const& massic_entropy)
{
    if (pressure <= detail::p_s_623_15)
    {
        if (massic_entropy <= r1::massic_entropy_pt (pressure, saturation_temperature_p (pressure))) return 1;
        if (massic_entropy <= r2::massic_entropy_pt (pressure, saturation_temperature_p (pressure))) return 4;
        return 2;
    }
    if (massic_entropy <= r1::massic_entropy_pt (pressure, 623.15)) return 1;
    if (massic_entropy > r2::massic_entropy_pt (pressure, t_b23_p (pressure))) return 2;
    if (pressure <= saturation_pressure_s (massic_entropy)) return 4;
    return 3;
}

// Main API
// P-T
#define ISTO_IAPWS_R7_GEN_PT(NAME)                                                                      \
    constexpr auto                                                                                      \
NAME##_pt (auto const& pressure, auto const& temperature)                                               \
{                                                                                                       \
    switch (region_pt (pressure, temperature))                                                          \
    {                                                                                                   \
        case 1: return r1::NAME##_pt (pressure, temperature);                                           \
        case 2: return r2::NAME##_pt (pressure, temperature);                                           \
        case 3: return r3::NAME##_pt (pressure, temperature);                                           \
        case 4: throw not_implemented_e { /*"On the saturation line, this is probably unreachable"*/ }; \
        case 5: return r5::NAME##_pt (pressure, temperature);                                           \
        default: throw internal_error_e {};                                                             \
    }                                                                                                   \
}                                                                                                       \
    constexpr auto                                                                                      \
NAME##_tp (auto const& temperature, auto const& pressure)                                               \
{                                                                                                       \
    return NAME##_pt (pressure, temperature);                                                           \
}
ISTO_IAPWS_R7_GEN_PT(massic_volume)
ISTO_IAPWS_R7_GEN_PT(density)
ISTO_IAPWS_R7_GEN_PT(massic_enthalpy)
ISTO_IAPWS_R7_GEN_PT(massic_internal_energy)
ISTO_IAPWS_R7_GEN_PT(massic_entropy)
ISTO_IAPWS_R7_GEN_PT(massic_isobaric_heat_capacity)
ISTO_IAPWS_R7_GEN_PT(massic_isochoric_heat_capacity)
ISTO_IAPWS_R7_GEN_PT(speed_of_sound)
ISTO_IAPWS_R7_GEN_PT(isobaric_cubic_expansion_coefficient)
ISTO_IAPWS_R7_GEN_PT(isothermal_compressibility)
ISTO_IAPWS_R7_GEN_PT(relative_pressure_coefficient)
ISTO_IAPWS_R7_GEN_PT(isothermal_stress_coefficient)
ISTO_IAPWS_R7_GEN_PT(massic_gibbs_free_energy)
ISTO_IAPWS_R7_GEN_PT(massic_helmholtz_free_energy)
#undef ISTO_IAPWS_R7_GEN_PT

// H-S
    constexpr auto
pressure_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    switch (region_hs (massic_enthalpy, massic_entropy))
    {
        case 1: return r1::pressure_hs (massic_enthalpy, massic_entropy);
        case 2: return r2::pressure_hs (massic_enthalpy, massic_entropy);
        case 3: return r3::pressure_hs (massic_enthalpy, massic_entropy);
        case 4: throw not_implemented_e { /*TODO*/ };
        case 5: throw not_implemented_e {};
        default: throw internal_error_e {};
    }
}
    constexpr auto
pressure_sh (auto const& massic_entropy, auto const& massic_enthalpy)
{
    return pressure_hs (massic_enthalpy, massic_entropy);
}
    constexpr auto
temperature_hs (auto const& massic_enthalpy, auto const& massic_entropy)
{
    switch (region_hs (massic_enthalpy, massic_entropy))
    {
        case 1: return r1::temperature_hs (massic_enthalpy, massic_entropy);
        case 2: return r2::temperature_hs (massic_enthalpy, massic_entropy);
        case 3: return r3::temperature_hs (massic_enthalpy, massic_entropy);
        case 4: throw not_implemented_e { /*TODO*/ };
        case 5: throw not_implemented_e {};
        default: throw internal_error_e {};
    }
}
    constexpr auto
temperature_sh (auto const& massic_entropy, auto const& massic_enthalpy)
{
    return temperature_hs (massic_enthalpy, massic_entropy);
}
#define ISTO_IAPWS_R7_GEN_HS(NAME)                                  \
    constexpr auto                                                  \
NAME##_hs (auto const& massic_enthalpy, auto const& massic_entropy) \
{                                                                   \
    return NAME##_pt (                                              \
          pressure_hs    (massic_enthalpy, massic_entropy)          \
        , temperature_hs (massic_enthalpy, massic_entropy)          \
    );                                                              \
}                                                                   \
    constexpr auto                                                  \
NAME##_sh (auto const& massic_entropy, auto const& massic_enthalpy) \
{                                                                   \
    return NAME##_hs (massic_enthalpy, massic_entropy);             \
}
ISTO_IAPWS_R7_GEN_HS(massic_volume)
ISTO_IAPWS_R7_GEN_HS(density)
ISTO_IAPWS_R7_GEN_HS(massic_internal_energy)
ISTO_IAPWS_R7_GEN_HS(massic_isobaric_heat_capacity)
ISTO_IAPWS_R7_GEN_HS(speed_of_sound)
#undef ISTO_IAPWS_R7_GEN_HS

// P-H
    constexpr auto
temperature_ph (auto const& pressure, auto const& massic_enthalpy)
{
    switch (region_ph (pressure, massic_enthalpy))
    {
        case 1: return r1::temperature_ph (pressure, massic_enthalpy);
        case 2: return r2::temperature_ph (pressure, massic_enthalpy);
        case 3: return r3::temperature_ph (pressure, massic_enthalpy);
        case 4: throw not_implemented_e { /*TODO*/ };
        case 5: throw not_implemented_e {};
        default: throw internal_error_e {};
    }
}
    constexpr auto
temperature_hp (auto const& massic_enthalpy, auto const& pressure)
{
    return temperature_ph (pressure, massic_enthalpy);
}
#define ISTO_IAPWS_R7_GEN_PH(NAME)                            \
    constexpr auto                                            \
NAME##_ph (auto const& pressure, auto const& massic_enthalpy) \
{                                                             \
    return NAME##_pt (                                        \
          pressure                                            \
        , temperature_ph (pressure, massic_enthalpy)          \
    );                                                        \
}                                                             \
    constexpr auto                                            \
NAME##_hp (auto const& massic_enthalpy, auto const& pressure) \
{                                                             \
    return NAME##_ph (pressure, massic_enthalpy);             \
}
ISTO_IAPWS_R7_GEN_PH(density)
ISTO_IAPWS_R7_GEN_PH(massic_internal_energy)
ISTO_IAPWS_R7_GEN_PH(massic_entropy)
ISTO_IAPWS_R7_GEN_PH(massic_isobaric_heat_capacity)
ISTO_IAPWS_R7_GEN_PH(speed_of_sound)
#undef ISTO_IAPWS_R7_GEN_PH
    constexpr auto
massic_volume_ph (auto const& pressure, auto const& massic_enthalpy)
{
    if (region_ph (pressure, massic_enthalpy) == 3) 
        return r3::massic_volume_ph (pressure, massic_enthalpy);
    return massic_volume_pt (
          pressure
        , temperature_ph (pressure, massic_enthalpy)
    );
}
    constexpr auto
massic_volume_hp (auto const& massic_enthalpy, auto const& pressure)
{
    return massic_volume_ph (pressure, massic_enthalpy);
}

// P-S
    constexpr auto
temperature_ps (auto const& pressure, auto const& massic_entropy)
{
    switch (region_ps (pressure, massic_entropy))
    {
        case 1: return r1::temperature_ps (pressure, massic_entropy);
        case 2: return r2::temperature_ps (pressure, massic_entropy);
        case 3: return r3::temperature_ps (pressure, massic_entropy);
        case 4: throw not_implemented_e { /*TODO*/ };
        case 5: throw not_implemented_e {};
        default: throw internal_error_e {};
    }
}
    constexpr auto
temperature_sp (auto const& massic_entropy, auto const& pressure)
{
    return temperature_ps (pressure, massic_entropy);
}
#define ISTO_IAPWS_R7_GEN_PS(NAME)                           \
    constexpr auto                                           \
NAME##_ps (auto const& pressure, auto const& massic_entropy) \
{                                                            \
    return NAME##_pt (                                       \
          pressure                                           \
        , temperature_ps (pressure, massic_entropy)          \
    );                                                       \
}                                                            \
    constexpr auto                                           \
NAME##_sp (auto const& massic_entropy, auto const& pressure) \
{                                                            \
    return NAME##_ps (pressure, massic_entropy);             \
}
ISTO_IAPWS_R7_GEN_PS(density)
ISTO_IAPWS_R7_GEN_PS(massic_internal_energy)
ISTO_IAPWS_R7_GEN_PS(massic_enthalpy)
ISTO_IAPWS_R7_GEN_PS(massic_isobaric_heat_capacity)
ISTO_IAPWS_R7_GEN_PS(speed_of_sound)
#undef ISTO_IAPWS_R7_GEN_PS
    constexpr auto
massic_volume_ps (auto const& pressure, auto const& massic_entropy)
{
    if (region_ps (pressure, massic_entropy) == 3) return r3::massic_volume_ps (pressure, massic_entropy);
    return massic_volume_pt (
          pressure
        , temperature_ps (pressure, massic_entropy)
    );
}
    constexpr auto
massic_volume_sp (auto const& massic_entropy, auto const& pressure)
{
    return massic_volume_ps (pressure, massic_entropy);
}

} // inline namespace r7_97_2012
} // namespace calculisto::iapws::r7
// vim: foldmethod=marker
