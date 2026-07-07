#include <calculisto/array/array.hpp>
#include <type_traits>
    using calculisto::array::array_t;
#include <calculisto/auto_diff/dual.hpp>
    using calculisto::auto_diff::dual_t;


    namespace
calculisto::auto_diff
{

    template <
          std::size_t N
        , std::size_t M
        , class T
        , class U
        , class V = std::common_type_t <T, U>
    >
    auto
pow (dual_t <N, T> const& d, array_t <U, M> const& a)
{
        auto
    r = array_t <dual_t <N, V>, M> {};
    std::transform (
            a.begin (),
            a.end   (),
            r.begin (),
            [=](auto e){ return pow (d, e); }
    );
    return r;
}
    template <
          std::size_t N
        , std::size_t M
        , class T
        , class U
        , class V = std::common_type_t <T, U>
    >
    auto
operator - (dual_t <N, T> const& d, array_t <U, M> const& a)
{
        auto
    r = array_t <dual_t <N, V>, M> {};
    std::transform (
            a.begin (),
            a.end   (),
            r.begin (),
            [=](auto e){ return d - e; }
    );
    return r;
}
    template <
          std::size_t N
        , std::size_t M
        , class T
        , class U
        , class V = std::common_type_t <T, U>
    >
    auto
operator / (array_t <U, M> const& a, dual_t <N, T> const& d)
{
        auto
    r = array_t <dual_t <N, V>, M> {};
    std::transform (
            a.begin (),
            a.end   (),
            r.begin (),
            [=](auto e){ return e / d; }
    );
    return r;
}
    template <
          std::size_t N
        , std::size_t M
        , class T
        , class U
        , class V = std::common_type_t <T, U>
    >
    auto
operator * (array_t <U, M> const& a, dual_t <N, T> const& d)
{
        auto
    r = array_t <dual_t <N, V>, M> {};
    std::transform (
            a.begin (),
            a.end   (),
            r.begin (),
            [=](auto e){ return e * d; }
    );
    return r;
}
    template <
          std::size_t N
        , std::size_t M
        , class T
        , class U
        , class V = std::common_type_t <T, U>
    >
    auto
operator * (dual_t <N, T> const& d, array_t <dual_t <N, U>, M> const& a)
{
        auto
    r = array_t <dual_t <N, V>, M> {};
    std::transform (
            a.begin (),
            a.end   (),
            r.begin (),
            [=](auto e){ return d * e; }
    );
    return r;
}
    template <
          std::size_t N
        , std::size_t M
        , class T
        , class U
        , class V = std::common_type_t <T, U>
    >
    auto
operator / (array_t <dual_t <N, U>, M> const& a, dual_t <N, T> const& d)
{
        auto
    r = array_t <dual_t <N, V>, M> {};
    std::transform (
            a.begin (),
            a.end   (),
            r.begin (),
            [=](auto e){ return e / d; }
    );
    return r;
}
    template <
          std::size_t N
        , std::size_t M
        , class T
        , class U
        , class V = std::common_type_t <T, U>
    >
    auto
operator * (array_t <dual_t <N, U>, M> const& a, dual_t <N, T> const& d)
{
        auto
    r = array_t <dual_t <N, V>, M> {};
    std::transform (
            a.begin (),
            a.end   (),
            r.begin (),
            [=](auto e){ return e * d; }
    );
    return r;
}
} // namespace calculisto::auto_diff
