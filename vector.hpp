#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <cassert>
#include <cmath>

namespace linear_algebra {
    template<class T>
    constexpr bool is_vector_type = false;
    template<class T>
    constexpr bool is_vector_reference_type = false;
    template<class T>
    using referenced_type = typename T::referenced_type;

    template<class T>
    concept vector = is_vector_type<std::remove_cvref_t<T>> && requires (T t1, T t2) {
        t1.size();
        t1[0];
    };
    template<class Vector, class ElementType>
    concept vector_element_type = vector<Vector> && requires(Vector v) {
        { v[0] } -> std::convertible_to<ElementType>;
    };
    template<class T>
    concept vector_reference = is_vector_reference_type<std::remove_cvref_t<T>> && vector<referenced_type<std::remove_cvref_t<T>>>;
    template<class Vector, class VectorReference>
    concept vector_reference_type = vector<Vector> && vector_reference<VectorReference> && std::same_as<Vector, referenced_type<Vector>>;
    template<class T>
    concept vector_or_vector_reference = vector<T> || vector_reference<T>;

    // Support ranged for loop.
    template<vector_or_vector_reference Vector>
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
    auto begin(vector_or_vector_reference auto& v) {
        return vector_iterator{v, 0};
    }
    auto end(vector_or_vector_reference auto& v) {
        return vector_iterator{v, v.size()};
    }

    void foreach_index(vector_or_vector_reference auto&& v, auto&& f) {
        std::remove_cvref_t<decltype(v.size())> i = 0;
        while (i < v.size()) {
            f(i);
            ++i;
        }
    }
    template<vector Vector, class T>
        requires vector_element_type<Vector, T>
    Vector operator*(const Vector& lhs, const T& rhs) {
        auto res{ lhs };
        for (auto& e : res) {
            e *= rhs;
        }
        assert(res[0] == lhs[0] * rhs);
        return res;
    }
    template<vector_reference VectorRef, class T>
        requires vector_element_type<referenced_type<VectorRef>, T>
    auto operator*(const VectorRef& lhs, const T& rhs) {
        referenced_type<VectorRef> res{};
        foreach_index(res,
            [&res, &lhs, &rhs](auto i) {
                res[i] = lhs[i] * rhs;
            }
        );
        assert(res[0] == lhs[0] * rhs);
        return res;
    }
    template<vector Vector, class T>
        requires vector_element_type<Vector, T>
    Vector operator*(const T& lhs, const Vector& rhs) {
        return rhs * lhs;
    }
    auto&& operator*=(vector_or_vector_reference auto&& lhs, const auto& rhs) {
        for (auto& e : lhs) {
            e *= rhs;
        }
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
    auto&& operator+=(vector_or_vector_reference auto&& lhs, const vector_or_vector_reference auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        for (size_t i = 0; i < lhs.size(); i++) {
            lhs[i] += rhs[i];
        }
        return std::forward<decltype(lhs)>(lhs);
    }
    auto&& operator-=(vector_or_vector_reference auto&& lhs, const vector_or_vector_reference auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        for (size_t i = 0; i < lhs.size(); i++) {
            lhs[i] -= rhs[i];
        }
        return std::forward<decltype(lhs)>(lhs);
    }
    auto dot_product(const vector_or_vector_reference auto& lhs, const vector_or_vector_reference auto& rhs) {
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
    std::ostream& operator<<(std::ostream& out, const vector_or_vector_reference auto& v) {
        out << "[";
        auto size = v.size();
        for (size_t i = 0; i+1 < size; i++) {
            out << v[i] << ", ";
        }
        out << v[size - 1];
        out << "]";
        return out;
    }
    using std::sqrt;
    template<vector Vector>
    auto length(const Vector& v) {
        std::remove_cvref_t<decltype(v[0])> len{0};
        for (auto e : v) {
            len += e*e;
        }
        return sqrt(len);
    }
    template<vector Vector>
    auto normalize(Vector v) {
        return v / length(v);
    }
    template<vector Vector>
    auto construct_orthonormal_vectors(std::vector<Vector> vectors) {
        std::vector<Vector> orthonormal_vectors{};
        for (auto v : vectors) {
            for (auto q : orthonormal_vectors) {
                v -= q * dot_product(v, q);
            }
            Vector zero{};
            if (v != zero) {
                auto q = normalize(v);
                orthonormal_vectors.emplace_back(q);
            }
        }
        return orthonormal_vectors;
    }
    template<vector Vector>
    auto gram_schmidt(std::vector<Vector> vectors) {
        std::vector<Vector> orthonormal_vectors{};
        for (auto v : vectors) {
            for (auto q : orthonormal_vectors) {
                v -= q * dot_product(v, q);
            }
            auto q = normalize(v);
            orthonormal_vectors.emplace_back(q);
        }
        return orthonormal_vectors;
    }

    void swap(vector_or_vector_reference auto&& lhs, vector_or_vector_reference auto&& rhs) {
        foreach_index(lhs,
            [&lhs, &rhs](auto i) {
                using std::swap;
                swap(lhs[i], rhs[i]);
            }
        );
    }
}
