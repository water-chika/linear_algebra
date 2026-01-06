#pragma once

#include <vector>
#include "vector.hpp"

namespace linear_algebra {
    template<class T>
    class resizeable_vector {
    public:
        using element_type = T;
        __device__ __host__
        resizeable_vector() : m_elements{} {}
        __device__ __host__
        resizeable_vector(std::initializer_list<T> data) : m_elements{ data } {}
        __device__ __host__
        resizeable_vector(vector_or_vector_reference auto& v) : m_elements(v.size())
        {
            for (size_t i = 0; i < v.size(); i++) {
                m_elements[i] = v[i];
            }
        }
        __device__ __host__
        size_t size() const {
            return m_elements.size();
        }
        __device__ __host__
        auto&& operator[](size_t i) {
            return m_elements[i];
        }
        __device__ __host__
        auto&& operator[](size_t i) const {
            return m_elements[i];
        }
        /*auto&& operator[](this auto&& child, size_t i) {
            auto& self = parent_cast<resizeable_vector&>(child);
            return self.m_elements[i];
        }*/
    private:
        std::vector<T> m_elements;
    };
    template<class T>
    constexpr bool is_vector_type<resizeable_vector<T>> = true;
}
