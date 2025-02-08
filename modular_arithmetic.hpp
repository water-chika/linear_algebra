#pragma once

namespace linear_algebra{
    template<class T, T M>
    class modular_arithmetic {
    public:
        modular_arithmetic(T v) : m_value{mod(v,M)} {}
        operator T() {
            return m_value;
        }
    private:
        static T mod(T v, T m) {
            auto r = v%m;
            if (r < 0) {
                r += m;
            }
            return r;
        }
        T m_value;
    };
    template<class T, T M>
    auto operator+(modular_arithmetic<T,M> lhs, modular_arithmetic<T,M> rhs) {
        return modular_arithmetic<T,M>{static_cast<T>(lhs) + static_cast<T>(rhs)};
    }
    template<class T, T M>
    auto operator-(modular_arithmetic<T,M> lhs, modular_arithmetic<T,M> rhs) {
        return modular_arithmetic<T,M>{static_cast<T>(lhs) + (M - static_cast<T>(rhs))};
    }
    template<class T, T M>
    auto operator*(modular_arithmetic<T,M> lhs, modular_arithmetic<T,M> rhs) {
        return modular_arithmetic<T,M>{static_cast<T>(lhs) * static_cast<T>(rhs)};
    }
}
