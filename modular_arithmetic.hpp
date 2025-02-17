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

    template<class T, class M>
    struct inverse_return_type {
        T t;
        M m;
        using type = decltype((t*m+1)/t);
    };

    template<
        typename Value,
        typename Mod,
        typename Res = inverse_return_type<Value, Mod>::type
        >
    Res&& inverse(Value value, Mod m, Res&& res = Res{}) {
        if (value == 1 || m == 1) {
            res = 1;
        }
        else if (value == 0) {
            throw std::runtime_error{"0 does not have inverse"};
        }
        else {
            auto n = mod((value-1) * inverse(mod(m, value), value), value) ;
            res = (n*m+1) / value;
        }
        return std::forward<Res>(res);
    }
    template<
        typename Lhs,
        typename Rhs,
        typename Res = Lhs,
        typename Mod = Lhs>
    auto&& divides(Lhs lhs, Rhs rhs, Mod m, Res&& res = Res{}) {
        auto primes = factor_primes(rhs);
        res = lhs;
        for (auto prime : primes) {
            res *= inverse(prime, m);
        }
        return std::forward<Res>(res);
    }
}
