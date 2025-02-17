#pragma once
#include <matrix/matrix.hpp>
#include <vector/combined_reference_vector.hpp>

namespace linear_algebra{
    template<class T, class Vector>
    struct row_type {
        T& A;
        size_t row;
        using element_type = typename T::element_type;
        using referenced_type = Vector;
        size_t size() const {
            return A.size().get_column();
        }
        auto& operator[](size_t i) {
            return A[A.size().set_column(i).set_row(row)];
        }
        auto operator[](size_t i) const {
            return A[A.size().set_column(i).set_row(row)];
        }

        auto& operator=(const vector_or_vector_reference auto& rhs) {
            foreach_index(*this,
                    [this, &rhs](auto i) {
                        this->operator[](i) = rhs[i];
                    }
                    );
            return *this;
        }
    };
    template<class T, class U>
    constexpr bool is_vector_reference_type<row_type<T,U>> = true;

    template<concept_helper::matrix M1, concept_helper::matrix M2>
    class combined_reference_matrix {
    public:
        using element_type = typename M1::element_type;
        using index_type = typename M1::index_type;
        using left_matrix_type = M1;
        using right_matrix_type = M2;
        combined_reference_matrix(M1& left, M2& right)
            : left_half{left}, right_half{right}
        {
            assert(left_half.size().get_row()
                    == right_half.size().get_row());
        }
        auto size() {
            auto l_size = left_half.size();
            auto r_size = right_half.size();
            auto res = l_size;
            return res.set_column(l_size.get_column()+r_size.get_column());
        }
        auto& column(size_t i) {
            const auto l_size = left_half.size().get_column();
            if (i < l_size) {
                return left_half.column(i);
            }
            else {
                return right_half.column(i-l_size);
            }
        }
        auto row(size_t i) const {
            return
                row_type<
                    combined_reference_matrix,
                    fixsized_vector<
                        std::remove_cvref_t<decltype(left_half[{0, 0}])>,
                        fixsized_matrix_column_count<M1>() + fixsized_matrix_column_count<M2>()
                    >
                >
                {*this, i};
        }
        auto row(size_t i) {
            return
                row_type<
                    combined_reference_matrix,
                    fixsized_vector<
                        std::remove_cvref_t<decltype(left_half[{0, 0}])>,
                        fixsized_matrix_column_count<M1> + fixsized_matrix_column_count<M2>
                    >
                >
                {*this, i};
        }
        auto& operator[](index_type i) {
            auto l_size = left_half.size().get_column();
            if (i.get_column() < l_size) {
                return left_half[i];
            }
            else {
                return right_half[i.set_column(i.get_column()-l_size)];
            }
        }
    private:
        M1& left_half;
        M2& right_half;
    };

    template<
        typename Element_subtract = std::minus<void>,
        typename Element_multiplies = std::multiplies<void>,
        typename Element_divides = std::divides<void>,
        typename Element_inverse = inverse,
        typename Element_is_invertible = is_invertible>
    auto inverse(concept_helper::matrix auto A,
            Element_subtract element_subtract = std::minus<void>{},
            Element_multiplies element_multiplies = std::multiplies<void>{},
            Element_divides element_divides = std::divides<void>{},
            Element_inverse element_inverse = inverse{},
            Element_is_invertible element_is_invertible = is_invertible{}
            ) {
        std::remove_cvref_t<decltype(A)> res = I;
        auto A_I = combined_reference_matrix{A, res};
        do_back_substitution(
                do_eliminate(A_I,
                    element_subtract,
                    element_multiplies,
                    element_divides,
                    element_inverse,
                    element_is_invertible
                    ),
                element_subtract,
                element_multiplies,
                element_divides,
                element_inverse,
                element_is_invertible
                );
        return res;
    }
    template<class Matrix>
    concept is_combined_reference_matrix_type = concept_helper::matrix<Matrix> && requires (Matrix m) {
        typename Matrix::left_matrix_type;
        typename Matrix::right_matrix_type;
    };
}

