#include "vector.hpp"
#include "resizeable_vector.hpp"

namespace linear_algebra {
    template<class T>
    class reference_vector {
    public:
        using element_type = T::element_type;
        using referenced_type = resizeable_vector<element_type>;

        __device__ __host__
        reference_vector(T ref) : m_ref{ref} {}

        __device__ __host__
        auto& operator=(const vector_or_vector_reference auto& v) {
            return copy_vector(*this, v);
        }
        __device__ __host__
        auto& operator=(const reference_vector& v) {
            return copy_vector(*this, v);
        }
        __device__ __host__
        size_t size() const {
            return m_ref.size();
        }
        __device__ __host__
        auto& operator[](size_t i) {
            return m_ref[i];
        }
        __device__ __host__
        const auto& operator[](size_t i) const {
            return m_ref[i];
        }
    private:
        T m_ref;
    };
    template<class T>
    constexpr bool is_vector_reference_type<reference_vector<T>> = true;
}