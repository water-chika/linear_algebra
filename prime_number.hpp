#pragma once

namespace prime_number {
    constexpr auto debug_prime_number = false;
    template<
        typename Int,
        typename IntIterate = Int>
    constexpr bool is_prime(Int n) {
        using std::gcd;
        for (IntIterate i = 2; i*i < n; i++) {
            if (gcd(i, n) > 1) {
                return false;
            }
        }
        return n>=2;
    }
    template<class T>
    class prime_number {
    public:
        prime_number(T v) : m_value{v} {
            if (debug_prime_number) assert(is_prime(v));
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
}
