#pragma once

#include <vector>
#include "vector.hpp"

namespace linear_algebra {
    template<class T>
    class resizeable_vector {
    public:
        using element_type = T;
        resizeable_vector() : m_elements{} {}
        resizeable_vector(std::initializer_list<T> data) : m_elements{ data } {}
        resizeable_vector(vector_or_vector_reference auto& v) : m_elements(v.size())
        {
            for (size_t i = 0; i < v.size(); i++) {
                m_elements[i] = v[i];
            }
        }
        size_t size() const {
            return m_elements.size();
        }
        auto&& operator[](this auto&& child, size_t i) {
            auto& self = parent_cast<resizeable_vector&>(child);
            return self.m_elements[i];
        }
    private:
        std::vector<T> m_elements;
    };
    template<class T>
    constexpr bool is_vector_type<resizeable_vector<T>> = true;

    template<class T>
    class reference_vector {
    public:
        using element_type = T::element_type;
        using referenced_type = resizeable_vector<element_type>;
        reference_vector(T ref) : m_ref{ref} {}
        size_t size() const {
            return m_ref.size();
        }
        auto& operator[](size_t i) {
            return m_ref[i];
        }
        const auto& operator[](size_t i) const {
            return m_ref[i];
        }
    private:
        T m_ref;
    };
    template<class T>
    constexpr bool is_vector_reference_type<reference_vector<T>> = true;
}
