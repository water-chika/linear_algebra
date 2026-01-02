#pragma once

#include "vector.hpp"

#include "type_cast.hpp"

#include <cpp_helper.hpp>

namespace linear_algebra {
    template<class T, size_t Size>
    class fixsized_vector {
    public:
        using element_type = T;
        template<typename T_, size_t N>
        using array = cpp_helper::array<T_, N>;

        constexpr fixsized_vector() : m_data{} {}
        constexpr fixsized_vector(std::initializer_list<T> data) : m_data{ array_from_initializer_list(data) } {}
        template<vector_or_vector_reference Vector>
        constexpr fixsized_vector(Vector v) : m_data{}
        {
            assert(v.size() == Size);
            for (size_t i = 0; i < Size; i++) {
                m_data[i] = v[i];
            }
        }
        constexpr size_t size() const {
            return Size;
        }
        constexpr auto&& operator[](size_t i) {
            return m_data[i];
        }
        constexpr auto&& operator[](size_t i) const {
            return m_data[i];
        }
        /*constexpr auto&& operator[](this auto&& self, size_t i) {
            return std::forward_like<decltype(self)>(
                        parent_cast<fixsized_vector&>(self)
                        .m_data[i]
                    );
        }*/
    private:
        __device__ __host__
        array<T, Size> array_from_initializer_list(std::initializer_list<T> data) {
            array<T, Size> res{};
            auto ite = data.begin();
            for (uint32_t i = 0; i < data.size(); i++) {
                res[i] = *ite;
                ++ite;
            }
            return res;
        }
        array<T, Size> m_data;
    };
    template<class T, size_t Size>
    class fixsized_pointer_vector {
    public:
        using element_type = T;
        using referenced_type = fixsized_vector<T, Size>;
        fixsized_pointer_vector() : m_data{} {}
        fixsized_pointer_vector(std::initializer_list<T*> data) : m_data{ array_from_initializer_list(data) } {}

        template<vectorlike Vector>
        auto& operator=(const Vector& v) {
            for (size_t i = 0; i < size(); i++) {
                this->operator[](i) = v[i];
            }
            return *this;
        }
        size_t size() const {
            return Size;
        }
        void set(size_t i, T* p) {
            m_data[i] = p;
        }
        constexpr auto&& operator[](size_t i) {
            return *m_data[i];
        }
        constexpr auto&& operator[](size_t i) const {
            return *m_data[i];
        }
        /*auto&& operator[](this auto&& child, size_t i) {
            auto& self = parent_cast<fixsized_pointer_vector&>(child);
            return *self.m_data[i];
        }*/
    private:
        std::array<T*, Size> array_from_initializer_list(std::initializer_list<T*> data) {
            std::array<T*, Size> res{};
            std::ranges::copy(data, res.begin());
            return res;
        }
        std::array<T*, Size> m_data;
    };

    template<class T, size_t Size>
    constexpr bool is_vector_type<fixsized_vector<T, Size>> = true;
    template<class T, size_t Size>
    constexpr bool is_vector_reference_type<fixsized_pointer_vector<T, Size>> = true;

}
