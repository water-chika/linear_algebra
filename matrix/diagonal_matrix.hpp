#pragma once

#include <matrix/matrix.hpp>

namespace linear_algebra{
    template<class T, size_t SIZE>
    class diagonal_matrix {
    public:
        using element_type = T;
        diagonal_matrix() = default;
        auto size() {
            return std::pair<size_t, size_t>{SIZE, SIZE};
        }
        auto operator[](std::pair<size_t, size_t> index) {
            if (index.first == index.second) {
                return m_diagonals[index.first];
            }
            else {
                return static_cast<T>(0);
            }
        }
        static diagonal_matrix from_diagonals(std::array<T, SIZE> diagonals) {
            diagonal_matrix D;
            D.m_diagonals = diagonals;
            return D;
        }
    private:
        std::array<T, SIZE> m_diagonals;
    };
    template<class T, size_t SIZE>
    diagonal_matrix<T,SIZE> make_diagonal_matrix(std::array<T, SIZE> diagonals) {
        return diagonal_matrix<T,SIZE>::from_diagonals(diagonals);
    }
    template<linear_algebra::concept_helper::matrix Matrix, class T, size_t SIZE>
    auto operator*(
            Matrix&& A,
            diagonal_matrix<T, SIZE> D) {
        auto B = A;
        for (size_t i = 0; i < SIZE; i++) {
            auto d = D[{i,i}];
            B.column(i) *= d;
        }
        return B;
    }
}
