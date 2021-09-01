#pragma once
#include "detail/common.hpp"

// TODO template <class T, class U, class V = std::common_type_t <T, U>>
// and long double coeffs
// ?

    namespace 
isto::iapws::r7
{
    inline namespace 
r7_97_2012
{

    struct internal_error_e {};

#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
#define ISTO_IAPWS_U_GC * unit::joule <> / unit::kelvin <> / unit::kilogram <>
#define ISTO_IAPWS_U_T  * unit::kelvin <>
#define ISTO_IAPWS_U_P  * unit::pascal <>
#define ISTO_IAPWS_U_D  * unit::kilogram <> / pow <3> (unit::metre <>)
#define ISTO_IAPWS_U_H  * unit::joule <> / unit::kilogram <>
#define ISTO_IAPWS_U_S  * unit::joule <> / unit::kilogram <> / unit::kelvin <>
#else
#define ISTO_IAPWS_U_GC
#define ISTO_IAPWS_U_T
#define ISTO_IAPWS_U_P
#define ISTO_IAPWS_U_D
#define ISTO_IAPWS_U_H
#define ISTO_IAPWS_U_S
#endif

    constexpr auto
massic_gas_constant = 0.461526e3 ISTO_IAPWS_U_GC;

    constexpr auto
critical_temperature = 647.096 ISTO_IAPWS_U_T;

    constexpr auto
critical_pressure = 22.064e6 ISTO_IAPWS_U_P;

    constexpr auto
critical_density = 322.0 ISTO_IAPWS_U_D;

#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
#define ISTO_IAPWS_P pressure_t <T>
#define ISTO_IAPWS_T temperature_t <T>
#define ISTO_IAPWS_D density_t <T>
#define ISTO_IAPWS_H massic_enthalpy_t <T>
#define ISTO_IAPWS_S massic_entropy_t <T>
#else
#define ISTO_IAPWS_P T
#define ISTO_IAPWS_T T
#define ISTO_IAPWS_D T
#define ISTO_IAPWS_H T
#define ISTO_IAPWS_S T
#endif

    namespace
r4
{
    namespace
detail
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
} // namespace r4

    template <class T>
    constexpr auto
saturation_pressure_t (ISTO_IAPWS_T const& temperature)
{
        using std::pow;
        using std::sqrt;
        using namespace r4::detail;
        auto const
    t = temperature / (1 ISTO_IAPWS_U_T);
        auto const
    theta = t + n[8] / (t - n[9]);
        auto const
    A = theta * theta + n[0] * theta + n[1];
        auto const
    B = n[2] * theta * theta + n[3] * theta + n[4];
        auto const
    C = n[5] * theta * theta + n[6] * theta + n[7];
    return pow (2. * C / (-B + sqrt (B * B - 4. * A * C)), 4) * 1e6 ISTO_IAPWS_U_P;
}

    template <class T>
    constexpr auto
saturation_temperature_p (ISTO_IAPWS_P const& pressure)
{
        using std::pow;
        using std::sqrt;
        using namespace r4::detail;
        auto const
    beta = pow (pressure / (1e6 ISTO_IAPWS_U_P), 0.25);
        auto const
    E = beta * beta + n[2] * beta + n[5];
        auto const
    F = n[0] * beta * beta + n[3] * beta + n[6];
        auto const
    G = n[1] * beta * beta + n[4] * beta + n[7];
        auto const
    D = 2. * G / (-F - sqrt (F * F - 4. * E * G));
    return (n[9] + D - sqrt (pow (n[9] + D, 2) - 4. * (n[8] + n[9] * D))) / 2. ISTO_IAPWS_U_T;
}
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T>
    constexpr auto
saturation_pressure (temperature_t <T> const& temperature)
{
    return saturation_pressure_t (temperature);
}
    template <class T>
    constexpr auto
saturation_temperature (pressure_t <T> const& pressure)
{
    return saturation_temperature_p (pressure);
}
#endif
// §4 Auxiliary Equation for the Boundary between Regions 2 and 3.
    template <class T>  
    constexpr auto
b23_t (ISTO_IAPWS_T const& temperature)
{
        constexpr double
    n1 =  0.34805185628969e3;
        constexpr double 
    n2 = -0.11671859879975e1;
        constexpr double 
    n3 =  0.10192970039326e-2;
        auto const
    t = temperature / (1 ISTO_IAPWS_U_T);
    return (n1 + n2 * t + n3 * t * t) * (1e6 ISTO_IAPWS_U_P);
}

    template <class T>
    /*constexpr*/ auto // Someday, std::sqrt might be constexpr.
b23_p (ISTO_IAPWS_P const& pressure)
{
        using std::sqrt;
        constexpr double 
    n3 = 0.10192970039326e-2;
        constexpr double 
    n4 = 0.57254459862746e3;
        constexpr double 
    n5 = 0.13918839778870e2;
        auto const
    p = pressure / (1e6 ISTO_IAPWS_U_P);
    return (n4 + sqrt ((p - n5) / n3)) * (1 ISTO_IAPWS_U_T);
}
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T>
    constexpr auto
b23 (temperature_t <T> const& temperature)
{
    return b23_t (temperature);
}
    template <class T>
    constexpr auto
b23 (pressure_t <T> const& pressure)
{
    return b23_p (pressure);
}
#endif

    template <class T>
    constexpr int
region_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    if (temperature <= (623.15 ISTO_IAPWS_U_T))
    {
            auto const
        ps = saturation_pressure_t (temperature);
        if (pressure > ps) return 1;
        if (pressure < ps) return 2;
        return 4;
    }
    if (temperature < 1073.15 ISTO_IAPWS_U_T)
    {
            auto const
        p23 = b23_t (temperature);
        if (pressure < p23) return 2;
        return 3;
    }
    return 5;
}
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T>
    constexpr int
region (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return region_pt (pressure, temperature);
}
#endif
// §5  Equations for Region 1.
    namespace
r1 
{
    constexpr array_t <double, 34>
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

    constexpr array_t <double, 34>
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

    constexpr array_t <double, 34>
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
    template <class T>
    T
gamma (const T& pi, const T& tau)
{
        using std::pow;
    return sum (
          n 
        * pow (7.1 - pi, I) 
        * pow (tau - 1.222, J)
    );
}
    template <class T>
    T
gamma_p (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
        - n 
        * I * pow (7.1 - pi, I - 1.) 
        * pow (tau - 1.222, J)
    );
}
    template <class T>
    T
gamma_t (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
          n 
        * pow (7.1 - pi, I) 
        * J * pow (tau - 1.222, J - 1.)
    );
}

    template <class T>
    T
gamma_pp (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
          n * I 
        * (I - 1.) * pow (7.1 - pi, I - 2.) 
        * pow (tau - 1.222, J)
    );
}
    template <class T>
    T
gamma_tt (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
          n 
        * pow (7.1 - pi, I) 
        * J * (J - 1.) * pow (tau - 1.222, J - 2.)
    );
}
    template <class T>
    T
gamma_pt (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
        - n * I
        * pow (7.1 - pi, I - 1.) 
        * J * pow (tau - 1.222, J - 1.)
    );
}
// }}}

    template <class T>
    auto
massic_volume_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (16.53 * (1e6 ISTO_IAPWS_U_P));
        auto const
    tau = (1386. * (1 ISTO_IAPWS_U_T)) / temperature;
    return pi * gamma_p (pi, tau) * massic_gas_constant * temperature / pressure;
}
    template <class T>
    auto
massic_enthalpy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (16.53 * (1e6 ISTO_IAPWS_U_P));
        auto const
    tau = 1386. * (1 ISTO_IAPWS_U_T) / temperature;
    return tau * gamma_t (pi, tau) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_internal_energy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (16.53 * (1e6 ISTO_IAPWS_U_P));
        auto const
    tau = 1386. * (1 ISTO_IAPWS_U_T) / temperature;
    return (tau * gamma_t (pi, tau) - pi * gamma_p (pi, tau)) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_entropy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (16.53 * (1e6 ISTO_IAPWS_U_P));
        auto const
    tau = 1386. * (1 ISTO_IAPWS_U_T) / temperature;
    return (tau * gamma_t (pi, tau) - gamma (pi, tau)) * massic_gas_constant;
}
    template <class T>
    auto
massic_isobaric_heat_capacity_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (16.53 * (1e6 ISTO_IAPWS_U_P));
        auto const
    tau = 1386. * (1 ISTO_IAPWS_U_T) / temperature;
    return - tau * tau * gamma_tt (pi, tau) * massic_gas_constant;
}
    template <class T>
    auto
massic_isochoric_heat_capacity_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (16.53 * (1e6 ISTO_IAPWS_U_P));
        auto const
    tau = 1386. * (1 ISTO_IAPWS_U_T) / temperature;
    return (- tau * tau * gamma_tt (pi, tau) + pow <2> (gamma_p (pi, tau) - tau * gamma_pp (pi, tau)) / gamma_pp (pi, tau)) * massic_gas_constant;
}
    template <class T>
    auto
speed_of_sound_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        using std::pow;
        using std::sqrt;
        auto const
    pi = pressure / (16.53 * (1e6 ISTO_IAPWS_U_P));
        auto const
    tau = 1386. * (1 ISTO_IAPWS_U_T) / temperature;
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
    template <class T>
    auto
temperature_ph (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_H const& enthalpy)
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
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    eta = enthalpy / (2500e3 ISTO_IAPWS_U_H);
    return sum (n * pow (pi, I) * pow (eta + 1., J)) * (1 ISTO_IAPWS_U_T);
}
    template <class T>
    auto
temperature_ps (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_S const& entropy)
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
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    sigma = entropy / (1e3 ISTO_IAPWS_U_S);
    return sum (n * pow (pi, I) * pow (sigma + 2., J)) * (1 ISTO_IAPWS_U_T);
}
    template <class T>
    auto
pressure_hs (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
        constexpr auto
    I = array_t <double, 19> //{{{
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
    J = array_t <double, 19> //{{{
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
    n = array_t <double, 19> //{{{
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
    eta = massic_enthalpy / (3400e3 ISTO_IAPWS_U_H);
        auto const
    sigma = massic_entropy / (7.6e3 ISTO_IAPWS_U_S);
    return sum (n * pow (eta + 0.05, I) * pow (sigma + 0.05, J)) * (100. ISTO_IAPWS_U_P);
}
    template <class T>
    constexpr auto
temperature_hs (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
    return temperature_ph (pressure_hs (massic_enthalpy, massic_entropy), massic_enthalpy);
}
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T>
    auto
massic_volume (pressure_t <T> const& pressure, temperature_t <T> const& temperature)
{
    return massic_volume_pt (pressure, temperature);
}
    template <class T>
    auto
massic_enthalpy (pressure_t <T> const& pressure, temperature_t <T> const& temperature)
{
    return massic_enthalpy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_internal_energy (pressure_t <T> const& pressure, temperature_t <T> const& temperature)
{
    return massic_internal_energy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_entropy (pressure_t <T> const& pressure, temperature_t <T> const& temperature)
{
    return massic_entropy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_isobaric_heat_capacity (pressure_t <T> const& pressure, temperature_t <T> const& temperature)
{
    return massic_isobaric_heat_capacity_pt (pressure, temperature);
}
    template <class T>
    auto
massic_isochoric_heat_capacity (pressure_t <T> const& pressure, temperature_t <T> const& temperature)
{
    return massic_isochoric_heat_capacity_pt (pressure, temperature);
}
    template <class T>
    auto
speed_of_sound (pressure_t <T> const& pressure, temperature_t <T> const& temperature)
{
    return speed_of_sound_pt (pressure, temperature);
}
    template <class T>
    constexpr auto
temperature (pressure_t <T> const& pressure, ISTO_IAPWS_H const& massic_enthalpy)
{
    return temperature_ph (pressure, massic_enthalpy);
}
    template <class T>
    constexpr auto
temperature (pressure_t <T> const& pressure, ISTO_IAPWS_S const& massic_entropy)
{
    return temperature_ps (pressure, massic_entropy);
}
    template <class T>
    auto
pressure (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
    return pressure_hs (massic_enthalpy, massic_entropy);
}
    template <class T>
    constexpr auto
temperature (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
    return temperature_hs (massic_enthalpy, massic_entropy);
}
#endif
} // namespace r1
// §5  Equations for Region 2.
    namespace
r2
{
    constexpr array_t <double, 9>
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
    constexpr array_t <double, 9>
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
    constexpr array_t <double, 43>
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
    constexpr array_t <double, 43>
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
    constexpr array_t <double, 43>
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
    template <class T>
    T
gamma_0 (const T& pi, const T& tau)
{
        using std::pow;
        using std::log;
    return log (pi) + sum (
          n_0
        * pow (tau, J_0)
    );
}
    template <class T>
    T
gamma_0_p (T const& pi, T const& /*tau*/)
{
    return static_cast <T> (1) / pi;
}
    template <class T>
    T
gamma_0_pp (T const& pi, T const& /*tau*/)
{
    return - static_cast <T> (1) / pi / pi;
}
    template <class T>
    T
gamma_0_t (T const& /*pi*/, T const& tau)
{
        using std::pow;
    return sum (
          n_0 
        * J_0 
        * pow (tau, J_0 - 1.)
    );
}
    template <class T>
    T
gamma_0_tt (T const& /*pi*/, T const& tau)
{
        using std::pow;
    return sum (
          n_0 
        * J_0 
        * (J_0 - 1.)
        * pow (tau, J_0 - 2.)
    );
}
    template <class T>
    T
gamma_0_pt (T const& pi, T const& tau)
{
    return static_cast <T> (0);
}
    template <class T>
    T
gamma_r (const T& pi, const T& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * pow (tau - 0.5, J_r)
    );
}
    template <class T>
    T
gamma_r_p (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
          n_r
        * I_r 
        * pow (pi, I_r - 1)
        * pow (tau - 0.5, J_r)
    );
}
    template <class T>
    T
gamma_r_pp (T const& pi, T const& tau)
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
    template <class T>
    T
gamma_r_t (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * J_r 
        * pow (tau - 0.5, J_r - 1.)
    );
}
    template <class T>
    T
gamma_r_tt (T const& pi, T const& tau)
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
    template <class T>
    T
gamma_r_pt (T const& pi, T const& tau)
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

    template <class T>
    T
gamma (const T& pi, const T& tau)
{
    return gamma_0 (pi, tau) + gamma_r (pi, tau);
}
    template <class T>
    T
gamma_p (T const& pi, T const& tau)
{
    return gamma_0_p (pi, tau) + gamma_r_p (pi, tau);
}
    template <class T>
    T
gamma_t (T const& pi, T const& tau)
{
    return gamma_0_t (pi, tau) + gamma_r_t (pi, tau);
}
    template <class T>
    T
gamma_pp (T const& pi, T const& tau)
{
    return gamma_0_pp (pi, tau) + gamma_r_pp (pi, tau);
}
    template <class T>
    T
gamma_tt (T const& pi, T const& tau)
{
    return gamma_0_tt (pi, tau) + gamma_r_tt (pi, tau);
}
    template <class T>
    T
gamma_pt (T const& pi, T const& tau)
{
    return gamma_0_pt (pi, tau) + gamma_r_pt (pi, tau);
}
// }}}
    template <class T>
    auto
massic_volume_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return pi * gamma_p (pi, tau) * massic_gas_constant * temperature / pressure;
}
    template <class T>
    auto
massic_enthalpy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return tau * gamma_t (pi, tau) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_internal_energy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return (tau * gamma_t (pi, tau) - pi * gamma_p (pi, tau)) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_entropy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return (tau * gamma_t (pi, tau) - gamma (pi, tau)) * massic_gas_constant;
}
    template <class T>
    auto
massic_isobaric_heat_capacity_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return - tau * tau * gamma_tt (pi, tau) * massic_gas_constant;
}
    template <class T>
    auto
massic_isochoric_heat_capacity_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return (- tau * tau * gamma_tt (pi, tau) + pow <2> (gamma_p (pi, tau) - tau * gamma_pp (pi, tau)) / gamma_pp (pi, tau)) * massic_gas_constant;
}
    template <class T>
    auto
speed_of_sound_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        using std::pow;
        using std::sqrt;
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
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
    template <class T>
    constexpr auto
b2bc_h (ISTO_IAPWS_H const& massic_enthalpy)
{
        auto const
    eta = massic_enthalpy / (1e3 ISTO_IAPWS_U_H);
    return (0.90584278514723e3 - 0.67955786399241 * eta + 0.12809002730136e-3 * eta * eta) * (1e6 ISTO_IAPWS_U_P);
}
    template <class T>
    constexpr auto
b2bc_p (ISTO_IAPWS_P const& pressure)
{
        using std::sqrt;
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
    return (0.26526571908428e4 + sqrt ((pi - 0.45257578905948e1 ) / 0.12809002730136e-3)) * (1e3 ISTO_IAPWS_U_H);
}
    template <class T>
    constexpr auto
region_ph (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_H const& massic_enthalpy)
{
    if (pressure < 4. * (1e6 ISTO_IAPWS_U_P)) return 1;
    if (pressure < b2bc_h (massic_enthalpy)) return 2;
    return 3;
}
    template <class T>
    constexpr auto
region_ps (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_S const& massic_entropy)
{
    if (pressure < 4. * (1e6 ISTO_IAPWS_U_P)) return 1;
    if (massic_entropy >= 5.85 * (1e3 ISTO_IAPWS_U_S)) return 2;
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

    template <class T>
    constexpr auto
temperature_ph (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_H const& massic_enthalpy)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    eta = massic_enthalpy / (2000e3 ISTO_IAPWS_U_H);
        using namespace t_ph;
    switch (region_ph (pressure, massic_enthalpy))
    {
        case 1: return sum (a::n * pow (pi,       a::I) * pow (eta - 2.1, a::J)) ISTO_IAPWS_U_T;
        case 2: return sum (b::n * pow (pi - 2.,  b::I) * pow (eta - 2.6, b::J)) ISTO_IAPWS_U_T;
        case 3: return sum (c::n * pow (pi + 25., c::I) * pow (eta - 1.8, c::J)) ISTO_IAPWS_U_T;
    }
    throw internal_error_e {};
}
    template <class T>
    constexpr auto
temperature_ps (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_S const& massic_entropy)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        using namespace t_ps;
    switch (region_ps (pressure, massic_entropy))
    {
        case 1: 
        {
                auto const
            sigma = massic_entropy / (2e3 ISTO_IAPWS_U_S);
            return sum (a::n * pow (pi, a::I) * pow (sigma - 2., a::J))  ISTO_IAPWS_U_T;
        }
        case 2: 
        {
                auto const
            sigma = massic_entropy / (0.7853e3 ISTO_IAPWS_U_S);
            return sum (b::n * pow (pi, b::I) * pow (10. - sigma, b::J)) ISTO_IAPWS_U_T;
        }
        case 3: 
        {
                auto const
            sigma = massic_entropy / (2.9251e3 ISTO_IAPWS_U_S);
            return sum (c::n * pow (pi, c::I) * pow (2. - sigma, c::J))  ISTO_IAPWS_U_T;
        }
    }
    throw internal_error_e {};
}
    template <class T>
    constexpr auto
b2ab_s (ISTO_IAPWS_S const& massic_entropy)
{
        auto const
    sigma = massic_entropy / (1e3 ISTO_IAPWS_U_S);
    return (
          -0.349898083432139e4
        +  0.257560716905876e4 * sigma
        + -0.421073558227969e3 * sigma * sigma
        +  0.276349063799944e2 * sigma * sigma * sigma
    ) * (1e3 ISTO_IAPWS_U_H);
}
    template <class T>
    constexpr auto
region_hs (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
    if (massic_enthalpy < b2ab_s (massic_entropy)) return 1;
    if (massic_entropy >= 5.85 * (1e3 ISTO_IAPWS_U_S)) return 2;
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
    template <class T>
    constexpr auto
pressure_hs (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
        using namespace p_hs;
    switch (region_hs (massic_enthalpy, massic_entropy))
    {
        case 1: 
        {
                auto const
            eta = massic_enthalpy / (4200e3 ISTO_IAPWS_U_H);
                auto const
            sigma = massic_entropy / (12e3 ISTO_IAPWS_U_S);
            return pow (sum (a::n * pow (eta - 0.5, a::I) * pow (sigma - 1.2, a::J)), 4) * (4e6 ISTO_IAPWS_U_P);
        }
        case 2: 
        {
                auto const
            eta = massic_enthalpy / (4100e3 ISTO_IAPWS_U_H);
                auto const
            sigma = massic_entropy / (7.9e3 ISTO_IAPWS_U_S);
            return pow (sum (b::n * pow (eta - 0.6, b::I) * pow (sigma - 1.01, b::J)), 4) * (100e6 ISTO_IAPWS_U_P);
        }
        case 3: 
        {
                auto const
            eta = massic_enthalpy / (3500e3 ISTO_IAPWS_U_H);
                auto const
            sigma = massic_entropy / (5.9e3 ISTO_IAPWS_U_S);
            return pow (sum (c::n * pow (eta - 0.7, c::I) * pow (sigma - 1.1, c::J)), 4) * (100e6 ISTO_IAPWS_U_P);
        }
    }
    throw internal_error_e {};
}
    template <class T>
    constexpr auto
temperature_hs (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
    temperature_ph (pressure_hs (massic_enthalpy, massic_entropy), massic_enthalpy);
}
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T>
    auto
massic_volume (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_volume_pt (pressure, temperature);
}
    template <class T>
    auto
massic_enthalpy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_enthalpy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_internal_energy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_internal_energy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_entropy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_entropy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_isobaric_heat_capacity (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_isobaric_heat_capacity_pt (pressure, temperature);
}
    template <class T>
    auto
massic_isochoric_heat_capacity (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_isochoric_heat_capacity_pt (pressure, temperature);
}
    template <class T>
    auto
speed_of_sound (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return speed_of_sound_pt (pressure, temperature);
}
    template <class T>
    constexpr auto
temperature (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_H const& massic_enthalpy)
{
    return temperature_ph (pressure, massic_enthalpy);
}
    template <class T>
    constexpr auto
temperature (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_S const& massic_entropy)
{
    return temperature_ps (pressure, massic_entropy);
}
    template <class T>
    constexpr auto
pressure (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
    return pressure_hs (massic_enthalpy, massic_entropy);
}
    template <class T>
    constexpr auto
temperature (ISTO_IAPWS_H const& massic_enthalpy, ISTO_IAPWS_S const& massic_entropy)
{
    temperature_hs (massic_enthalpy, massic_entropy);
}
#endif
} //namespace r2

    namespace
r2_metastable_vapor
{
    constexpr array_t <double, 9>
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
    constexpr array_t <double, 13>
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
    constexpr array_t <double, 13>
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
    constexpr array_t <double, 13>
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
    template <class T>
    T
gamma_0 (const T& pi, const T& tau)
{
        using std::pow;
        using std::log;
    return log (pi) + sum (
          n_0
        * pow (tau, J_0)
    );
}
    template <class T>
    T
gamma_0_p (T const& pi, T const& /*tau*/)
{
    return static_cast <T> (1) / pi;
}
    template <class T>
    T
gamma_0_pp (T const& pi, T const& /*tau*/)
{
    return - static_cast <T> (1) / pi / pi;
}
    template <class T>
    T
gamma_0_t (T const& /*pi*/, T const& tau)
{
        using std::pow;
    return sum (
          n_0 
        * J_0 
        * pow (tau, J_0 - 1.)
    );
}
    template <class T>
    T
gamma_0_tt (T const& /*pi*/, T const& tau)
{
        using std::pow;
    return sum (
          n_0 
        * J_0 
        * (J_0 - 1.)
        * pow (tau, J_0 - 2.)
    );
}
    template <class T>
    T
gamma_0_pt (T const& pi, T const& tau)
{
    return static_cast <T> (0);
}
    template <class T>
    T
gamma_r (const T& pi, const T& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * pow (tau - 0.5, J_r)
    );
}
    template <class T>
    T
gamma_r_p (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
          n_r
        * I_r 
        * pow (pi, I_r - 1)
        * pow (tau - 0.5, J_r)
    );
}
    template <class T>
    T
gamma_r_pp (T const& pi, T const& tau)
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
    template <class T>
    T
gamma_r_t (T const& pi, T const& tau)
{
        using std::pow;
    return sum (
          n_r
        * pow (pi, I_r)
        * J_r 
        * pow (tau - 0.5, J_r - 1.)
    );
}
    template <class T>
    T
gamma_r_tt (T const& pi, T const& tau)
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
    template <class T>
    T
gamma_r_pt (T const& pi, T const& tau)
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

    template <class T>
    T
gamma (const T& pi, const T& tau)
{
    return gamma_0 (pi, tau) + gamma_r (pi, tau);
}
    template <class T>
    T
gamma_p (T const& pi, T const& tau)
{
    return gamma_0_p (pi, tau) + gamma_r_p (pi, tau);
}
    template <class T>
    T
gamma_t (T const& pi, T const& tau)
{
    return gamma_0_t (pi, tau) + gamma_r_t (pi, tau);
}
    template <class T>
    T
gamma_pp (T const& pi, T const& tau)
{
    return gamma_0_pp (pi, tau) + gamma_r_pp (pi, tau);
}
    template <class T>
    T
gamma_tt (T const& pi, T const& tau)
{
    return gamma_0_tt (pi, tau) + gamma_r_tt (pi, tau);
}
    template <class T>
    T
gamma_pt (T const& pi, T const& tau)
{
    return gamma_0_pt (pi, tau) + gamma_r_pt (pi, tau);
}
// }}}
    template <class T>
    auto
massic_volume_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return pi * gamma_p (pi, tau) * massic_gas_constant * temperature / pressure;
}
    template <class T>
    auto
massic_enthalpy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return tau * gamma_t (pi, tau) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_internal_energy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return (tau * gamma_t (pi, tau) - pi * gamma_p (pi, tau)) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_entropy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return (tau * gamma_t (pi, tau) - gamma (pi, tau)) * massic_gas_constant;
}
    template <class T>
    auto
massic_isobaric_heat_capacity_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return - tau * tau * gamma_tt (pi, tau) * massic_gas_constant;
}
    template <class T>
    auto
massic_isochoric_heat_capacity_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
    return (- tau * tau * gamma_tt (pi, tau) + pow <2> (gamma_p (pi, tau) - tau * gamma_pp (pi, tau)) / gamma_pp (pi, tau)) * massic_gas_constant;
}
    template <class T>
    auto
speed_of_sound_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        using std::pow;
        using std::sqrt;
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 540 * (1 ISTO_IAPWS_U_T) / temperature;
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
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T>
    auto
massic_volume (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_volume_pt (pressure, temperature);
}
    template <class T>
    auto
massic_enthalpy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_enthalpy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_internal_energy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_internal_energy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_entropy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_entropy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_isobaric_heat_capacity (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_isobaric_heat_capacity_pt (pressure, temperature);
}
    template <class T>
    auto
massic_isochoric_heat_capacity (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_isochoric_heat_capacity_pt (pressure, temperature);
}
    template <class T>
    auto
speed_of_sound (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return speed_of_sound_pt (pressure, temperature);
}
#endif
} // namespace r2_metastable_vapor


// §5  Equations for Region 3.
    namespace
r3
{
    constexpr array_t <double, 39>
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
    constexpr array_t <double, 39>
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
    constexpr array_t <double, 39>
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

// low-level dimensionless {{{
    template <class T>
    auto
phi (T const& delta, T const& tau)
{
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
    template <class T>
    auto
phi_d (T const& delta, T const& tau)
{
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
    template <class T>
    auto
phi_dd (T const& delta, T const& tau)
{
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
    template <class T>
    auto
phi_t (T const& delta, T const& tau)
{
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
    template <class T>
    auto
phi_tt (T const& delta, T const& tau)
{
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
    template <class T>
    auto
phi_dt (T const& delta, T const& tau)
{
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
    template <class T>
    auto
pressure_dt (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (delta * phi_d (delta, tau)) * density * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_internal_energy_dt (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (tau * phi_t (delta, tau)) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_entropy_dt (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (tau * phi_t (delta, tau) - phi (delta, tau)) * massic_gas_constant;
}
    template <class T>
    auto
massic_enthalpy_dt (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (tau * phi_t (delta, tau) + delta * phi_d (delta, tau)) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_isochoric_heat_capacity_dt (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (- tau * tau * phi_tt (delta, tau)) * massic_gas_constant;
}
    template <class T>
    auto
massic_isobaric_heat_capacity_dt (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
        auto const
    delta = density / critical_density;
        auto const
    tau = critical_temperature / temperature;
    return (- tau * tau * phi_tt (delta, tau) + pow (delta * phi_d (delta, tau) - delta * tau * phi_dt (delta, tau), 2) / (2. * delta * phi_d (delta, tau) + delta * delta * phi_dd (delta, tau))) * massic_gas_constant;
}
    template <class T>
    auto
speed_of_sound_dt (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
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
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T>
    auto
pressure (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
    return pressure_dt (density, temperature);
}
    template <class T>
    auto
massic_internal_energy (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
    return massic_internal_energy_dt (density, temperature);
}
    template <class T>
    auto
massic_entropy (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
    return massic_entropy_dt (density, temperature);
}
    template <class T>
    auto
massic_enthalpy (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
    return massic_enthalpy_dt (density, temperature);
}
    template <class T>
    auto
massic_isochoric_heat_capacity (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
    return massic_isobaric_heat_capacity_dt (density, temperature);
}
    template <class T>
    auto
massic_isobaric_heat_capacity (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
    return massic_isobaric_heat_capacity_dt (density, temperature);
}
    template <class T>
    auto
speed_of_sound (ISTO_IAPWS_D const& density, ISTO_IAPWS_T const& temperature)
{
    return speed_of_sound_dt (density, temperature);
}
#endif
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
    template <class T>
    constexpr auto
gamma_0 (T const& pi, T const& tau)
{
        using std::log;
    return log (pi) + sum (n_0 * pow (tau, J_0));
}
    template <class T>
    constexpr auto
gamma_0_p (T const& pi, T const& /* tau */)
{
    return static_cast <T> (1) / pi;
}
    template <class T>
    constexpr auto
gamma_0_pp (T const& pi, T const& /* tau */)
{
    return - static_cast <T> (1) / pi / pi;
}
    template <class T>
    constexpr auto
gamma_0_t (T const& /* pi */, T const& tau)
{
    return sum (n_0 * J_0 * pow (tau, J_0 - 1.));
}
    template <class T>
    constexpr auto
gamma_0_tt (T const& /* pi */, T const& tau)
{
    return sum (n_0 * J_0 * (J_0 - 1.) * pow (tau, J_0 - 2.));
}
    template <class T>
    constexpr auto
gamma_0_pt (T const& /* pi */, T const& /* tau */)
{
    return static_cast <T> (0);
}
    template <class T>
    constexpr auto
gamma_r (T const& pi, T const& tau)
{
    return sum (n_r * pow (pi, I_r) * pow (tau, J_r));
}
    template <class T>
    constexpr auto
gamma_r_p (T const& pi, T const& tau)
{
    return sum (n_r * I_r * pow (pi, I_r - 1.) * pow (tau, J_r));
}
    template <class T>
    constexpr auto
gamma_r_pp (T const& pi, T const& tau)
{
    return sum (n_r * I_r * (I_r - 1.) * pow (pi, I_r - 2.) * pow (tau, J_r));
}
    template <class T>
    constexpr auto
gamma_r_t (T const& pi, T const& tau)
{
    return sum (n_r * pow (pi, I_r) * J_r * pow (tau, J_r - 1.));
}
    template <class T>
    constexpr auto
gamma_r_tt (T const& pi, T const& tau)
{
    return sum (n_r * pow (pi, I_r) * J_r * (J_r - 1.) * pow (tau, J_r - 2.));
}
    template <class T>
    constexpr auto
gamma_r_pt (T const& pi, T const& tau)
{
    return sum (n_r * I_r * pow (pi, I_r - 1.) * J_r * pow (tau, J_r - 1.));
}
    template <class T>
    T
gamma (const T& pi, const T& tau)
{
    return gamma_0 (pi, tau) + gamma_r (pi, tau);
}
    template <class T>
    T
gamma_p (T const& pi, T const& tau)
{
    return gamma_0_p (pi, tau) + gamma_r_p (pi, tau);
}
    template <class T>
    T
gamma_t (T const& pi, T const& tau)
{
    return gamma_0_t (pi, tau) + gamma_r_t (pi, tau);
}
    template <class T>
    T
gamma_pp (T const& pi, T const& tau)
{
    return gamma_0_pp (pi, tau) + gamma_r_pp (pi, tau);
}
    template <class T>
    T
gamma_tt (T const& pi, T const& tau)
{
    return gamma_0_tt (pi, tau) + gamma_r_tt (pi, tau);
}
    template <class T>
    T
gamma_pt (T const& pi, T const& tau)
{
    return gamma_0_pt (pi, tau) + gamma_r_pt (pi, tau);
}
// }}}
} // namespace detail

    using namespace detail;

    template <class T>
    auto
massic_volume_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 1000. * (1 ISTO_IAPWS_U_T) / temperature;
    return pi * gamma_p (pi, tau) * massic_gas_constant * temperature / pressure;
}
    template <class T>
    auto
massic_enthalpy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 1000. * (1 ISTO_IAPWS_U_T) / temperature;
    return tau * gamma_t (pi, tau) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_internal_energy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 1000. * (1 ISTO_IAPWS_U_T) / temperature;
    return (tau * gamma_t (pi, tau) - pi * gamma_p (pi, tau)) * massic_gas_constant * temperature;
}
    template <class T>
    auto
massic_entropy_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 1000. * (1 ISTO_IAPWS_U_T) / temperature;
    return (tau * gamma_t (pi, tau) - gamma (pi, tau)) * massic_gas_constant;
}
    template <class T>
    auto
massic_isobaric_heat_capacity_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 1000. * (1 ISTO_IAPWS_U_T) / temperature;
    return - tau * tau * gamma_tt (pi, tau) * massic_gas_constant;
}
    template <class T>
    auto
massic_isochoric_heat_capacity_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 1000. * (1 ISTO_IAPWS_U_T) / temperature;
    return (- tau * tau * gamma_tt (pi, tau) + pow <2> (gamma_p (pi, tau) - tau * gamma_pp (pi, tau)) / gamma_pp (pi, tau)) * massic_gas_constant;
}
    template <class T>
    auto
speed_of_sound_pt (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
        using std::pow;
        using std::sqrt;
        auto const
    pi = pressure / (1e6 ISTO_IAPWS_U_P);
        auto const
    tau = 1000. * (1 ISTO_IAPWS_U_T) / temperature;
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
#ifdef ISTO_IAPWS_FLAVOR_CONSTRAINED
    template <class T>
    auto
massic_volume (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_volume_pt (pressure, temperature);
}
    template <class T>
    auto
massic_enthalpy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_enthalpy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_internal_energy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_internal_energy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_entropy (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_entropy_pt (pressure, temperature);
}
    template <class T>
    auto
massic_isobaric_heat_capacity (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_isobaric_heat_capacity_pt (pressure, temperature);
}
    template <class T>
    auto
massic_isochoric_heat_capacity (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return massic_isochoric_heat_capacity_pt (pressure, temperature);
}
    template <class T>
    auto
speed_of_sound (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return speed_of_sound_pt (pressure, temperature);
}
#endif
} // namespace r5
#if 0

#define ISTO_THERMODYNAMICS_IAPWSR7_GEN(NAME)                                                           \
    template <class T>                                                                                  \
    auto                                                                                                \
NAME (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)                   \
{                                                                                                       \
        using namespace detail;                                                                         \
    switch (region (pressure, temperature))                                                             \
    {                                                                                                   \
        case 1: return r1::NAME (pressure, temperature);                                                \
        case 2: return r2::NAME (pressure, temperature);                                                \
        case 3: throw not_implemented_e { /*"You are in region 3, inversion needed"*/ };                \
        case 4: throw not_implemented_e { /*"On the saturation line, this is probably unreachable"*/ }; \
        case 5: return r5::NAME (pressure, temperature);                                                \
        default: throw not_yet_implemented_e {};                                                        \
    }                                                                                                   \
}                                                                                                       \
    template <class T>                                                                                  \
    auto                                                                                                \
NAME (ISTO_IAPWS_T const& temperature, ISTO_IAPWS_P const& pressure)                   \
{                                                                                                       \
    return NAME (pressure, temperature);                                                                \
}                                                                                                       \

ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_volume)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_enthalpy)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_internal_energy)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_entropy)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_isobaric_heat_capacity)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(massic_isochoric_heat_capacity)
ISTO_THERMODYNAMICS_IAPWSR7_GEN(speed_of_sound)

#undef ISTO_THERMODYNAMICS_IAPWSR7_GEN

    template <class T>
    auto
density (ISTO_IAPWS_P const& pressure, ISTO_IAPWS_T const& temperature)
{
    return pow <-1> (massic_volume (pressure, temperature));
}
    template <class T>
    auto
density (ISTO_IAPWS_T const& temperature, ISTO_IAPWS_P const& pressure)
{
    return density (pressure, temperature);
}

#endif

#undef ISTO_IAPWS_P
#undef ISTO_IAPWS_T
#undef ISTO_IAPWS_D
#undef ISTO_IAPWS_H
#undef ISTO_IAPWS_S

#undef ISTO_IAPWS_U_GC
#undef ISTO_IAPWS_U_T
#undef ISTO_IAPWS_U_P
#undef ISTO_IAPWS_U_D
#undef ISTO_IAPWS_U_H
#undef ISTO_IAPWS_U_S

} // inline namespace r7_97_2012
} // namespace isto::iapws::r7
// vim: foldmethod=marker
