#pragma once
#include <vector.hpp>

namespace linear_algebra {
    template<vector V1, vector V2>
    class combined_reference_vector {
    public:
        combined_reference_vector(
                V1& left, V2& right)
            : left_half{left}, right_half{right}
        {}
        size_t size() const {
            return left_half.size() + right_half.size();
        }
        auto& operator[](size_t i) {
            auto l_size = left_half.size();
            if (i < l_size) {
                return left_half[i];
            }
            else {
                return right_half[i-l_size];
            }
        }
    private:
        V1& left_half;
        V2& right_half;
    };
}
