#pragma once
#include <matrix.hpp>
#include <combined_reference_vector.hpp>

namespace linear_algebra{
    template<class T, class Vector>
    struct row_type {
        T& A;
        size_t row;
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
    };
    template<class T, class U>
    constexpr bool is_vector_reference_type<row_type<T,U>> = true;

    template<class T>
    constexpr size_t fixsized_matrix_column_count = 0;
    template<class T, size_t COLUMN, size_t ROW>
    constexpr size_t fixsized_matrix_column_count<fixsized_matrix<T, COLUMN, ROW>> = COLUMN;
    template<class T, size_t COLUMN, size_t ROW>
    constexpr size_t fixsized_matrix_column_count<matrix<T, COLUMN, ROW>> = COLUMN;

    template<concept_helper::matrix M1, concept_helper::matrix M2>
    class combined_reference_matrix {
    public:
        using index_type = typename M1::index_type;
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

    auto inverse(concept_helper::matrix auto A) {
        std::remove_cvref_t<decltype(A)> res = I;
        auto A_I = combined_reference_matrix{A, res};
        do_back_substitution(do_eliminate(A_I));
        return res;
    }
}

