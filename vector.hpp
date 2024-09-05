#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <cassert>

#include <fixsized_vector.hpp>

namespace linear_algebra {
    template<class T>
    constexpr bool is_vector_type = false;
    template<class T, size_t Size>
    constexpr bool is_vector_type<fixsized_vector<T, Size>> = true;
    template<class T>
    concept vector = is_vector_type<std::remove_cvref_t<T>> && requires (T t1, T t2) {
        t1.size();
        t1[0];
    };
    template<class Vector, class ElementType>
    concept vector_element_type = vector<Vector> && requires(Vector v) {
        { v[0] } -> std::convertible_to<ElementType>;
    };

    // Support ranged for loop.
    template<vector Vector>
    struct vector_iterator {
        Vector& v;
        size_t index;
        vector_iterator& operator++() {
            ++index;
            return *this;
        }
        vector_iterator operator++(int) {
            auto res = *this;
            index++;
            return res;
        }
        auto& operator*() {
            return v[index];
        }
        bool operator==(const vector_iterator& lhs) const {
            return (&v) == (&lhs.v) && index == lhs.index;
        }
    };
    auto begin(vector auto&& v) {
        return vector_iterator{v, 0};
    }
    auto end(vector auto&& v) {
        return vector_iterator{v, v.size()};
    }

    template<vector Vector, class T>
        requires vector_element_type<Vector, T>
    Vector operator*(const Vector& lhs, const T& rhs) {
        Vector res{ lhs };
        for (auto& e : res) {
            e *= rhs;
        }
        return res;
    }
    template<vector Vector, class T>
        requires vector_element_type<Vector, T>
    Vector operator*(const T& lhs, const Vector& rhs) {
        return rhs * lhs;
    }
    template<vector Vector, class T>
        requires vector_element_type<Vector, T>
    Vector& operator*=(Vector& lhs, const T& rhs) {
        lhs = lhs * rhs;
        return lhs;
    }
    bool operator==(const vector auto& lhs, const vector auto& rhs) {
        bool equal{lhs.size() == rhs.size()};
        for (size_t i = 0; equal && i < lhs.size(); i++) {
            equal = lhs[i] == rhs[i];
        }
        return equal;
    }
    vector auto operator+(const vector auto& lhs, const vector auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        auto res{lhs};
        for (size_t i = 0; i < lhs.size(); i++) {
            res[i] += rhs[i];
        }
        return res;
    }
    vector auto operator-(const vector auto& lhs, const vector auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        auto res{ lhs };
        for (size_t i = 0; i < res.size(); i++) {
            res[i] -= rhs[i];
        }
        return res;
    }
    vector auto operator-(const vector auto& lhs) {
        auto res{ lhs };
        for (size_t i = 0; i < res.size(); i++) {
            res[i] = -lhs[i];
        }
        return res;
    }
    auto dot_product(const vector auto& lhs, const vector auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        auto res{ lhs[0]*rhs[0] };
        for (size_t i = 1; i < lhs.size(); i++) {
            res += lhs[i] * rhs[i];
        }
        return res;
    }
    auto element_multi(const vector auto& lhs, const vector auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        auto res{ lhs };
        for (size_t i = 0; i < lhs.size(); i++) {
            res[i] *= rhs[i];
        }
        return res;
    }
    template<vector Vector, class T>
        requires vector_element_type<Vector, T>
    Vector operator/(const Vector& lhs, const T& rhs) {
        Vector res{ lhs };
        for (size_t i = 0; i < res.size(); i++) {
            res[i] /= rhs;
        }
        return res;
    }
    std::ostream& operator<<(std::ostream& out, const vector auto& v) {
        out << "[";
        auto size = v.size();
        for (size_t i = 0; i+1 < size; i++) {
            out << v[i] << ", ";
        }
        out << v[size - 1];
        out << "]";
        return out;
    }
}
