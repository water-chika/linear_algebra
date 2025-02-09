#pragma once

#include <numeric>
#include <cassert>
#include <optional>
#include <iostream>
#include <vector>
#include<source_location>

namespace linear_algebra{
    constexpr bool debug_modular_arithmetic =
        true;
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
        bool operator==(const prime_number& rhs) const = default;
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
    template<class T>
    auto operator%(const auto& lhs, const prime_number<T>& rhs) {
        return lhs % static_cast<T>(rhs);
    }
    template<class T>
    auto operator%(const prime_number<T>& lhs, const auto& rhs) {
        return lhs % static_cast<T>(rhs);
    }
    template<class T,class U>
    auto operator%(const prime_number<T>& lhs, const prime_number<U>& rhs) {
        return static_cast<T>(lhs) % static_cast<U>(rhs);
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
    class modular_arithmetic {
    public:
        constexpr modular_arithmetic(T v, M m) : m_value{v}, m_mod{m}{
            if (debug_modular_arithmetic) assert(m_value >= 0 && m_value < m_mod);
        }
        explicit constexpr operator T() const {
            return m_value;
        }
        constexpr auto value() const {
            return m_value;
        }
        constexpr auto modular_value() const {
            return m_mod;
        }
        constexpr bool operator==(const modular_arithmetic& rhs) const = default;
    private:
        T m_value;
        M m_mod;
    };

    template<class T, class M>
    auto make_modular_arithmetic(T v, M m) {
        return modular_arithmetic{mod(v,m), m};
    }

    template<class T, class M>
    std::ostream& operator<<(std::ostream& out, modular_arithmetic<T, M> v) {
        return out << static_cast<T>(v) << "(mod " << v.modular_value() << ")";
    }

    template<class T1, class M1, class T2, class M2>
    auto modular_value(modular_arithmetic<T1, M1> lhs, modular_arithmetic<T2, M2> rhs) {
        if (debug_modular_arithmetic) assert(lhs.modular_value() == rhs.modular_value());
        return lhs.modular_value();
    }
    template<class T1, class M1, class T2, class M2>
    auto binary_op(const modular_arithmetic<T1, M1>& lhs, const modular_arithmetic<T2, M2>& rhs,
            auto&& op) {
        auto m = modular_value(lhs, rhs);
        return make_modular_arithmetic(
                mod(op(lhs.value(), rhs.value()), m),
                m
        );
    }
    template<class T, class M>
    auto operator+(const modular_arithmetic<T,M> lhs, const modular_arithmetic<T,M> rhs) {
        return binary_op(lhs, rhs,
                [](auto lhs, auto rhs) {
                    return lhs + rhs;
                }
                );
    }
    template<class T, class M>
    auto& operator+=(modular_arithmetic<T,M>& lhs, const modular_arithmetic<T,M>& rhs) {
        lhs = lhs + rhs;
        return lhs;
    }
    template<class T, class M>
    auto operator-(modular_arithmetic<T,M> lhs, modular_arithmetic<T,M> rhs) {
        return binary_op(lhs, rhs,
                [m=rhs.modular_value()](auto lhs, auto rhs) {
                    return lhs + m - rhs;
                }
                );
    }
    template<class T1, class M1, class T2, class M2>
    auto operator*(const modular_arithmetic<T1,M1>& lhs, const modular_arithmetic<T2,M2>& rhs) {
        return binary_op(lhs, rhs,
                [](auto lhs, auto rhs) {
                    return lhs * rhs;
                }
                );
    }
    template<class T, class M>
    auto& operator*=(modular_arithmetic<T,M>& lhs, modular_arithmetic<T,M> rhs) {
        lhs = lhs * rhs;
        return lhs;
    }

    template<class T, class M>
    struct inverse_return_type {
        T t;
        M m;
        using type = decltype((t*m+1)/t);
    };

    constexpr auto debug_modular_arithmetic_inverse = true;
    template<class T, class M>
    auto inverse(modular_arithmetic<prime_number<T>, prime_number<M>> v)
        ->
        modular_arithmetic<typename inverse_return_type<T,M>::type, prime_number<M>> {
        if (debug_modular_arithmetic_inverse) {
            std::cout << "inverse(" << v << ")" << std::endl;
            std::cout << std::source_location::current().function_name() << std::endl;
        }
        auto mod_v = v.modular_value();
        auto p = static_cast<prime_number<T>>(v);
        auto m_p = make_modular_arithmetic(mod_v, p);
        auto n = inverse(m_p).value();
        auto res = modular_arithmetic{(n*mod_v+1)/p, mod_v};
        if (debug_modular_arithmetic_inverse) {
            assert((res * v).value() == 1);
        }
        return res;
    }
    template<class T, class M>
    auto inverse(modular_arithmetic<T, prime_number<M>> v) {
        auto primes = factor_primes(static_cast<T>(v));
        if (debug_modular_arithmetic_inverse) {
            std::cout << "inverse(" << v << ")" << std::endl;
            std::cout << "primes: ";
            for (auto prime : primes) {
                std::cout << prime << ", ";
            }
            std::cout << std::endl;
        }
        auto mod_v = v.modular_value();
        auto res = modular_arithmetic<T, prime_number<M>>{1, mod_v};
        for (auto prime : primes) {
            res *= inverse(modular_arithmetic<prime_number<T>, prime_number<M>>{prime, mod_v});
        }
        return res;
    }
    template<class T, class M>
    auto operator/(modular_arithmetic<T,prime_number<M>> lhs, modular_arithmetic<T,prime_number<M>> rhs) {
        auto res = lhs*inverse(rhs);
        if (debug_modular_arithmetic) {
            assert(res * rhs == lhs);
        }
        return res;
    }
}
