#pragma once
#define CHECK_U(a, b, eps)\
    static_assert (decltype (a)::dimension == decltype (b)::dimension); \
    CHECK(a.magnitude == Approx { b.magnitude }.epsilon (eps))

#define CHECK_US(a, b, sca, eps)\
    static_assert (decltype (a)::dimension == decltype (b)::dimension); \
    CHECK(a.magnitude == Approx { b.magnitude }.scale(sca).epsilon (eps))
