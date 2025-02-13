#pragma once

#include <numeric>
#include <cassert>
#include <optional>
#include <iostream>
#include <vector>
#include<source_location>

namespace modular_arithmetic {
    constexpr auto debug_modular_arithmetic = false;

    constexpr bool is_prime(auto n) {
        using std::gcd;
        for (decltype(n) i = 2; i*i < n; i++) {
            if (gcd(i, n) > 1) {
                return false;
            }
        }
        return n>=2;
    }
    template<class T, T V>
    constexpr bool is_prime(std::integral_constant<T, V> n) {
        using std::gcd;
        for (T i = 2; i*i < n; i++) {
            if (gcd(i, static_cast<T>(n)) > 1) {
                return false;
            }
        }
        return n>=2;
    }
    template<class T>
    class prime_number {
    public:
        prime_number(T v) : m_value{v} {
            if (debug_modular_arithmetic) assert(is_prime(v));
        }
        operator T() const {
            return m_value;
        }
        auto value() const {
            return m_value;
        }
    private:
        T m_value;
    };

    template<class T>
    std::ostream& operator<<(std::ostream& out, prime_number<T> p) {
        return out << static_cast<T>(p);
    }
    template<class T>
    auto operator+(const prime_number<T>& lhs, const auto& rhs) {
        return static_cast<T>(lhs) + rhs;
    }
    template<class T>
    auto operator+(const auto& lhs, const prime_number<T>& rhs) {
        return lhs + static_cast<T>(rhs);
    }
    template<class T>
    auto& operator+=(auto& lhs, const prime_number<T>& rhs) {
        return lhs += static_cast<T>(rhs);
    }
    template<class T>
    auto operator-(const prime_number<T>& lhs, const auto& rhs) {
        return static_cast<T>(lhs) - rhs;
    }
    template<class T>
    auto operator-(const auto& lhs, const prime_number<T>& rhs) {
        return lhs - static_cast<T>(rhs);
    }
    template<class T>
    auto operator*(const prime_number<T>& lhs, const auto& rhs) {
        return static_cast<T>(lhs) * rhs;
    }
    template<class T>
    auto operator*(const auto& lhs, const prime_number<T>& rhs) {
        return lhs * static_cast<T>(rhs);
    }
    template<class T,class U>
    auto operator*(const prime_number<T>& lhs, const prime_number<U>& rhs) {
        return static_cast<T>(lhs) * static_cast<U>(rhs);
    }
    template<class T>
    auto operator/(const auto& lhs, const prime_number<T>& rhs) {
        return lhs / static_cast<T>(rhs);
    }
    template<class T>
    auto operator/(const prime_number<T>& lhs, const auto& rhs) {
        return lhs / static_cast<T>(rhs);
    }
    template<class T,class U>
    auto operator/(const prime_number<T>& lhs, const prime_number<U>& rhs) {
        return static_cast<T>(lhs) / static_cast<U>(rhs);
    }
    template<class T,class U>
    auto operator%(const prime_number<T>& lhs, const prime_number<U>& rhs) {
        return static_cast<T>(lhs) % static_cast<U>(rhs);
    }
    template<class T>
    auto operator%(const auto& lhs, const prime_number<T>& rhs) {
        return lhs % static_cast<T>(rhs);
    }
    template<class T>
    auto operator%(const prime_number<T>& lhs, const auto& rhs) {
        return lhs.value() % rhs;
    }
    template<class T>
    auto operator<=>(const auto& lhs, const prime_number<T>& rhs) {
        return lhs <=> static_cast<T>(rhs);
    }
    template<class T>
    auto operator<(const prime_number<T>& lhs, const auto& rhs) {
        return lhs <=> static_cast<T>(rhs);
    }
    template<class T,class U>
    auto operator<=>(const prime_number<T>& lhs, const prime_number<U>& rhs) {
        return static_cast<T>(lhs) <=> static_cast<U>(rhs);
    }


    template<class T>
    constexpr bool is_prime(prime_number<T> p) {
        return true;
    }

    template<class T>
    std::vector<prime_number<T>> factor_primes(T n) {
        std::vector<prime_number<T>> primes{};

        T prime = 2;
        while (n > 1) {
            if (n % prime == 0) {
                primes.push_back(prime);
                n /= prime;
            }
            else {
                prime++;
            }
        }
        return primes;
    }

    constexpr auto mod(auto v, auto m) {
        auto r = v%m;
        if (r < 0) {
            r += m;
        }
        return r;
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
            auto n = inverse(mod(m, value), value);
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
