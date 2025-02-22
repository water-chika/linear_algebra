#pragma once

#include <numeric>
#include <cassert>
#include <optional>
#include <iostream>
#include <vector>
#include<source_location>

#include "prime_number.hpp"

namespace modular_arithmetic {
    using namespace prime_number;

    constexpr auto debug_modular_arithmetic = false;

    constexpr auto mod(auto v, auto m) {
        auto r = v%m;
        if (r < 0) {
            r += m;
        }
        return r;
    }

    auto is_invertible(auto v, auto m) {
        using std::gcd;
        return gcd(v, m) == 1;
    }

    template<class L, class R, class M>
    struct divides_return_type {
        L l;
        R r;
        M m;
        using type = decltype((l*m+l)/r);
    };
    template<
        typename Lhs,
        typename Rhs,
        typename Mod,
        typename Res = divides_return_type<Lhs, Rhs, Mod>::type
            >
    Res&& divides(Lhs lhs, Rhs rhs, Mod m, Res&& res = Res{}) {
        if (rhs == 1 || m == 1) {
            res = lhs;
        }
        else if (rhs == 0) {
            throw std::runtime_error{"0 does not have inverse"};
        }
        else {
            auto n = divides(mod(-lhs, rhs), mod(m, rhs), rhs);
            res = (n*m+lhs) / rhs;
        }
        return std::forward<Res>(res);
    }

    template<
        typename Value,
        typename Mod,
        typename Res = divides_return_type<Value, Value, Mod>::type
        >
    Res&& inverse(Value value, Mod m, Res&& res = Res{}) {
        return std::forward<Res>(res = divides(1, value, m));
    }
}
