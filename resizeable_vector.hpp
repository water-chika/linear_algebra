#pragma once

#include <vector>
#include "vector.hpp"

namespace linear_algebra {
    template<class T>
    class resizeable_vector {
    public:
        resizeable_vector() : m_elements{} {}
        resizeable_vector(std::initializer_list<T> data) : m_elements{ data } {}
        size_t size() const {
            return m_elements.size();
        }
        auto& operator[](size_t i) {
            return m_elements[i];
        }
        const auto& operator[](size_t i) const {
            return m_elements[i];
        }
    private:
        std::vector<T> m_elements;
    };
    template<class T>
    constexpr bool is_vector_type<resizeable_vector<T>> = true;
}