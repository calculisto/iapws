#pragma once
#define CHECK_U(a, b, eps)\
    static_assert (decltype (a)::dimension == decltype (b)::dimension); \
    CHECK(a.magnitude == Approx { b.magnitude }.epsilon (eps))
