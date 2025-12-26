#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <cassert>
#include <cmath>

#include <complex_number.hpp>

namespace linear_algebra {
    template<class T>
    constexpr bool is_vector_type = false;
    template<class T>
    constexpr bool is_vector_reference_type = false;
    template<class T>
    using referenced_type = typename T::referenced_type;

    template<class T>
    concept vectorlike = is_vector_type<std::remove_cvref_t<T>> && requires (T t1, T t2) {
        t1.size();
        t1[0];
    };
    template<class T>
    concept vector_reference = is_vector_reference_type<std::remove_cvref_t<T>> && vectorlike<referenced_type<std::remove_cvref_t<T>>>;
    template<class Vector, class VectorReference>
    concept vector_reference_type = vectorlike<Vector> && vector_reference<VectorReference> && std::same_as<Vector, referenced_type<Vector>>;
    template<class T>
    concept vector_or_vector_reference = vectorlike<T> || vector_reference<T>;

    template<class T>
    concept with_element_type = requires(T t) {
        typename T::element_type;
    };

    template<with_element_type T>
    struct element_type_struct {
    public:
        using type = typename T::element_type;
    };
    template<class T>
    using element_type = element_type_struct<std::remove_cvref_t<T>>::type;

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

    template<vector_or_vector_reference Vec>
    struct target_type_struct {
        using type = Vec;
    };
    template<vector_reference VecRef>
    struct target_type_struct<VecRef> {
        using type = referenced_type<VecRef>;
    };

    template<vector_or_vector_reference Vec>
    using target_type = typename target_type_struct<Vec>::type;

    template<
        vector_or_vector_reference Lhs,
        typename Res = target_type<Lhs>,
        typename Element_mul = std::multiplies<void>
            >
    auto multiplies(const Lhs& lhs, const element_type<Lhs>& rhs,
            Res&& res = Res{}, Element_mul&& element_mul = std::multiplies<void>{}) {
        res = lhs;
        foreach_index(res,
                [&res, &lhs, &rhs, &element_mul](auto i) {
                    res[i] = element_mul(lhs[i] , rhs);
                }
                );
        return res;
    }
    template<
        vector_or_vector_reference Lhs,
        vector_or_vector_reference Rhs,
        typename Res = target_type<Lhs>,
        typename Element_sub = std::minus<void>
            >
    auto subtract(const Lhs& lhs, const Rhs& rhs,
            Res&& res = Res{}, Element_sub&& element_sub = std::minus<void>{}) {
        res = lhs;
        foreach_index(res,
                [&res, &lhs, &rhs, &element_sub](auto i) {
                    res[i] = element_sub(lhs[i] , rhs[i]);
                }
                );
        return res;
    }
    template<vectorlike Vector, class T>
        requires std::same_as<element_type<Vector>, T>
    Vector operator*(const Vector& lhs, const T& rhs) {
        return multiplies(lhs, rhs);
    }
    template<vector_reference VectorRef, class T>
        requires std::same_as<element_type<VectorRef>, T>
    auto operator*(const VectorRef& lhs, const T& rhs) {
        return multiplies(lhs, rhs, referenced_type<VectorRef>{});
    }
    template<vectorlike Vector, class T>
        requires std::same_as<element_type<Vector>, T>
    Vector operator*(const T& lhs, const Vector& rhs) {
        return rhs * lhs;
    }
    template<vector_reference VectorRef, class T>
        requires std::same_as<element_type<VectorRef>, T>
    auto operator*(const T& lhs, const VectorRef& rhs) {
        return rhs * lhs;
    }
    auto&& operator*=(vector_or_vector_reference auto&& lhs, const auto& rhs) {
        for (auto& e : lhs) {
            e *= rhs;
        }
        return lhs;
    }
    bool operator==(const vectorlike auto& lhs, const vectorlike auto& rhs) {
        bool equal{lhs.size() == rhs.size()};
        for (size_t i = 0; equal && i < lhs.size(); i++) {
            equal = lhs[i] == rhs[i];
        }
        return equal;
    }
    vectorlike auto operator+(const vectorlike auto& lhs, const vectorlike auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        auto res{lhs};
        for (size_t i = 0; i < lhs.size(); i++) {
            res[i] += rhs[i];
        }
        return res;
    }
    vectorlike auto operator-(const vectorlike auto& lhs, const vectorlike auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        auto res{ lhs };
        for (size_t i = 0; i < res.size(); i++) {
            res[i] -= rhs[i];
        }
        return res;
    }
    vectorlike auto operator-(const vectorlike auto& lhs) {
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
    template<
        typename Element_add = std::plus<void>,
        typename Element_multiplies = std::multiplies<void>
        >
    auto dot_product(
            const vector_or_vector_reference auto& lhs,
            const vector_or_vector_reference auto& rhs,
            Element_add&& element_add = std::plus<void>{},
            Element_multiplies&& element_multiplies = std::multiplies<void>{}
            ) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        auto res{ element_multiplies(lhs[0], rhs[0]) };
        for (size_t i = 1; i < lhs.size(); i++) {
            res = element_add(res, element_multiplies(lhs[i], rhs[i]));
        }
        return res;
    }
    auto element_multi(const vectorlike auto& lhs, const vectorlike auto& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error{ "size not equal" };
        }
        auto res{ lhs };
        for (size_t i = 0; i < lhs.size(); i++) {
            res[i] *= rhs[i];
        }
        return res;
    }

    template<vector_or_vector_reference Vector, class T, vector_or_vector_reference VectorRes = Vector>
    auto divides(const Vector& lhs, const T& rhs, VectorRes&& res = {}) {
        res = lhs;
        for (size_t i = 0; i < res.size(); i++) {
            res[i] /= rhs;
        }
        return std::forward<VectorRes>(res);
    }
    template<vector_or_vector_reference Vector, class T>
        requires std::convertible_to<T, element_type<Vector>>
    Vector operator/(const Vector& lhs, const T& rhs) {
        Vector res{};
        return divides(lhs, rhs, res);
    }
    template<vectorlike Vector, class T>
        requires std::convertible_to<T, element_type<Vector>>
    Vector& operator/=(Vector& lhs, const T& rhs) {
        for (size_t i = 0; i < lhs.size(); i++) {
            lhs[i] /= rhs;
        }
        return lhs;
    }
    template<vector_reference Vector, class T>
        requires std::convertible_to<T, element_type<Vector>>
    Vector operator/=(Vector lhs, const T& rhs) {
        for (size_t i = 0; i < lhs.size(); i++) {
            lhs[i] /= rhs;
        }
        return lhs;
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


    using complex_number::length_square;
    template<vector_or_vector_reference Vector>
    auto length_square(const Vector& v) {
        element_type<Vector> t{0};
        auto len_2 = length_square(t);
        for (auto e : v) {
            len_2 += length_square(e);
        }
        return len_2;
    }

    using std::sqrt;
    template<vector_or_vector_reference Vector>
    auto length(const Vector& v) {
        return sqrt(length_square(v));
    }
    template<vector_or_vector_reference Vector>
    auto normalize(Vector v) {
        return v / length(v);
    }
    template<vector_or_vector_reference Vector>
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
    template<vector_or_vector_reference Vector>
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
