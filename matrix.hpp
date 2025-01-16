#pragma once
#include "vector.hpp"
#include "fixsized_vector.hpp"
#include <algorithm>
#include <iomanip>

namespace linear_algebra {
    namespace concept_helper {
        template <class Matrix_Size>
        concept matrix_size = requires (Matrix_Size matrix_size) {
            matrix_size.get_column();
            matrix_size.get_row();
        };
        template<class Matrix_Element>
        concept matrix_element = requires (Matrix_Element matrix_element) {
            matrix_element* matrix_element;
            matrix_element + matrix_element;
            matrix_element - matrix_element;
            matrix_element / matrix_element;
        };
        template<class Matrix_Index>
        concept matrix_index = requires (Matrix_Index matrix_index) {
            matrix_index.get_row();
            matrix_index.get_column();
        };
        template <class Matrix>
        concept matrix = requires (Matrix m, int i) {
            m.size();
            m.column(i);
            m.row(i);
        };
    }
    template<class T>
    struct matrix_index {
        constexpr auto set_row(T row) {
            m_row = row;
            return *this;
        }
        constexpr auto set_column(T column) {
            m_column = column;
            return *this;
        }
        constexpr auto get_row() const {
            return m_row;
        }
        constexpr auto get_column() const {
            return m_column;
        }
        T m_row;
        T m_column;
    };

    template<class T, size_t size>
    class column_vector : public fixsized_vector<T, size> {
    public:
        column_vector() : fixsized_vector<T, size>{} {}
        column_vector(fixsized_vector<T, size> v) : fixsized_vector<T, size>{ v } {}
    };

    struct identity_matrix_type {};
    identity_matrix_type I;

    template<class T>
    constexpr size_t fixsized_matrix_column_count = 0;

    template<class T>
    void iterate_from_0_to(T end, auto&& f) {
        for (T i = 0; i < end; i++) {
            f(i);
        }
    }
    template<concept_helper::matrix Matrix, class F>
        requires std::invocable<F, typename Matrix::index_type>
    void foreach_index(Matrix& A, F fn) {
        auto [row, column] = A.size();
        iterate_from_0_to(row,
            [&A, &fn, column](auto i) {
                iterate_from_0_to(column,
                    [&A, &fn, i](auto j) {
                        using index_t = typename Matrix::index_type;
                        fn(index_t{i, j});
                    });
            }
        );
    }
    template<class M>
    auto get_element() {
        M A;
        return A[typename M::index_type{}];
    }
    template<concept_helper::matrix Matrix>
    struct element_type_struct<Matrix> {
    public:
        using type = std::remove_cvref_t<typeof(get_element<Matrix>())>;
    };

    template<concept_helper::matrix Matrix, class F>
        requires std::invocable<F, element_type<Matrix>&>
    void foreach_element(Matrix& A, F fn) {
        auto [row, column] = A.size();
        iterate_from_0_to(row,
            [&A, &fn, column](auto i) {
                iterate_from_0_to(column,
                    [&A, &fn, i](auto j) {
                        using index_t = typename Matrix::index_type;
                        fn(A[index_t{i, j}]);
                    });
            }
        );
    }
    bool operator==(const concept_helper::matrix_index auto& i, const concept_helper::matrix_index auto& j) {
        return i.get_column() == j.get_column() && i.get_row() == j.get_row();
    }
    bool operator==(const concept_helper::matrix auto& A, const concept_helper::matrix auto& B) {
        if (A.size() != B.size()) {
            return false;
        }
        auto [row_count, column_count] = A.size();
        for (decltype(row_count) i = 0; i < row_count; i++) {
            for (decltype(column_count) j = 0; j < column_count; j++) {
                if (A[{i,j}] != B[{i,j}]) {
                    return false;
                }
            }
        }
        return true;
    }

    auto operator/(
            linear_algebra::concept_helper::matrix auto&& A,
            linear_algebra::concept_helper::matrix_element auto t) {
        auto B = A;
        foreach_element(B, [t](element_type<typeof(B)>& e) {
                e /= t;
                });
        return B;
    }

    auto& operator<<(std::ostream& out, linear_algebra::concept_helper::matrix auto&& A) {
        out << "{" << std::endl;
        auto [row_count, column_count] = A.size();
        iterate_from_0_to(row_count,
            [&out, &A, column_count](auto row) {
                out << "{";
                iterate_from_0_to(column_count,
                    [&out, &A, &row](auto column) {
                        out << std::setw(8) << A[{row, column}] << ",";
                    });
                out << "}," << std::endl;
            });
        out << "}";
        return out;
    }


    template<std::integral Int>
    void for_index_range(std::convertible_to<Int> auto&& begin, Int end, std::invocable<Int> auto&& f) {
        for (Int i = begin; i < end; i++) {
            f(i);
        }
    }

    auto& do_eliminate(concept_helper::matrix auto& A, concept_helper::matrix_index auto&& index) {
        auto column = index.get_column();
        auto row = index.get_row();
        auto pivot = A[{column, column}];
        if (pivot == static_cast<typeof(pivot)>(0)) {
            throw std::runtime_error{ "matrix is not invertible" };
        }
        auto multier = A[index] / pivot;
        A.row(row) -= A.row(column)*multier;
        return A;
    }

    template<concept_helper::matrix Matrix>
    auto& do_eliminate(Matrix& A) {
        auto [row_count, column_count] = A.size();
        for_index_range(0, std::min(row_count, column_count),
            [&A](auto column) {
                auto pivot = A[{column, column}];
                auto zero = static_cast<typeof(pivot)>(0);
                if (pivot == zero) {
                    for (decltype(column) row = column + 1; row < A.size().get_row(); row++) {
                        if (A[{row, column}] != zero) {
                            swap(A.row(column), A.row(row));
                            break;
                        }
                    }
                }
                if (A[{column, column}] != zero) {
                    for_index_range(column + 1, A.size().get_row(),
                        [&A, column](auto row) {
                            do_eliminate(A, typename Matrix::index_type{ row, column });
                        });
                }
            }
        );
        return A;
    }
    template<concept_helper::matrix Matrix>
    Matrix eliminate(Matrix A) {
        return do_eliminate(A);
    }
    
    template<concept_helper::matrix Matrix>
    Matrix& do_back_substitution(Matrix& A) {
        for (auto row = A.size().get_row() - 1; true; row--) {
            for (auto column = A.size().get_row() - 1; true; column--) {
                if (column != row) {
                    A.row(row) -= A.row(column) *
                        A[{row, column}];
                }
                else {
                    auto pivot = A[std::remove_cvref_t<decltype(A.size())>{}.set_row(row).set_column(row)];
                    auto zero = typeof(pivot){0};
                    if (pivot == zero) {
                        throw std::runtime_error{ "matrix is not invertible" };
                    }
                    auto one = typeof(pivot){1};
                    A.row(row) *= one / pivot;
                }
                if (column == row) {
                    break;
                }
            }
            if (row == 0) {
                break;
            }
        }
        return A;
    }
    template<concept_helper::matrix Matrix>
    Matrix back_substitution(Matrix A) {
        return do_back_substitution(A);
    }

    template<concept_helper::matrix Matrix>
    auto get_columns(Matrix&& A) {
        auto size = A.size().get_column();
        std::vector<std::remove_cvref_t<decltype(A[0])>> columns(size);
        iterate_from_0_to(size,
            [&columns, &A](auto i) {
                columns[i] = A.column(i);
            });
        return columns;
    }

    template<concept_helper::matrix Matrix>
    auto gram_schmidt(Matrix&& A) {
        auto res = A;
        std::remove_cvref_t<Matrix> R = I;
        iterate_from_0_to(A.size().get_column(),
            [&res, &R](auto i) {
                iterate_from_0_to(i,
                    [&res, &R, i](auto j) {
                        auto base = res.column(j);
                        auto b = res.column(i);
                        auto scale = dot_product(base, b)/dot_product(base,base);
                        res.column(i) -= scale*base;
                        R[{j,i}] = scale;
                    });
            });
        return std::pair{res, R};
    }

    template<concept_helper::matrix Matrix>
    auto determinant(Matrix A) {
        auto U = eliminate(A);
        std::remove_cvref_t<decltype(U[{0, 0}]) > res = 1;
        iterate_from_0_to(U.size().get_row(),
            [&res, &U](auto i) {
                res *= U[{i, i}];
            });
        return res;
    }

}
