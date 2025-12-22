#pragma once
#include "vector/vector.hpp"
#include "vector/fixsized_vector.hpp"
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
        concept matrix_element = requires (Matrix_Element element) {
            element * element;
            element + element;
            element - element;
            element / element;
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
    template <class Matrix>
    concept matrix = requires (Matrix m, int i) {
        m.size();
        m.column(i);
        m.row(i);
    };
    template<class T, class U>
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

        matrix_index() = default;

        matrix_index(T r, U c) 
        : m_row{r}, m_column{c}
        {}

        template<class T_, class U_>
        matrix_index(matrix_index<T_,U_> i) 
        : m_row{i.m_row}, m_column{i.m_column}
        {}

        T m_row;
        U m_column;
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
            [&fn, column](auto i) {
                iterate_from_0_to(column,
                    [&fn, i](auto j) {
                        using index_t = typename Matrix::index_type;
                        fn(index_t{i, j});
                    });
            }
        );
    }
    template<concept_helper::matrix_index MatrixIndex, class F>
        requires std::invocable<F, MatrixIndex>
    void foreach_index(MatrixIndex size, F fn) {
        auto [row, column] = size;
        iterate_from_0_to(row,
            [&fn, column](auto i) {
                iterate_from_0_to(column,
                    [&fn, i](auto j) {
                        fn(MatrixIndex{i, j});
                    });
            }
        );
    }

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
    template<concept_helper::matrix Matrix1, concept_helper::matrix Matrix2, class F>
        requires std::invocable<F, element_type<Matrix1>&, element_type<Matrix2>&>
    void foreach_element(Matrix1& A, Matrix2& B, F fn) {
        auto [row, column] = A.size();
        assert(A.size() == B.size());
        iterate_from_0_to(row,
            [&A, &B, &fn, column](auto i) {
                iterate_from_0_to(column,
                    [&A, &B, &fn, i](auto j) {
                        using index_t = typename Matrix1::index_type;
                        auto index = index_t{i,j};
                        fn(A[index], B[index]);
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

    constexpr auto debug_matrix = false;

    template<concept_helper::matrix MatrixLhs,
        concept_helper::matrix MatrixRhs,
        concept_helper::matrix MatrixRes,
        class ElementOp>
    auto element_op(const MatrixLhs& lhs, const MatrixRhs& rhs, ElementOp&& element_op) {
        if (debug_matrix) {
            assert(lhs.size() == rhs.size());
        }
        MatrixRes res{};
        foreach_index(res,
                [&res, &lhs, &rhs, &element_op](auto i) {
                    res[i] = element_op(lhs[i],rhs[i]);
                }
                );
        return res;
    }

    template<concept_helper::matrix MatrixLhs,
        concept_helper::matrix MatrixRhs,
        concept_helper::matrix MatrixRes = MatrixLhs,
        class ElementAdd = std::plus<void>>
    auto add(const MatrixLhs& lhs, const MatrixRhs& rhs, ElementAdd&& element_add) {
        return element_op(lhs, rhs, element_add);
    }
    template<concept_helper::matrix MatrixLhs,
        concept_helper::matrix MatrixRhs,
        concept_helper::matrix MatrixRes = MatrixLhs,
        class ElementSubtract = std::minus<void>>
    auto subtract(const MatrixLhs& lhs, const MatrixRhs& rhs, ElementSubtract&& element_sub) {
        return element_op(lhs, rhs, element_sub);
    }
    template<concept_helper::matrix MatrixLhs,
        concept_helper::matrix MatrixRhs,
        concept_helper::matrix MatrixRes = MatrixLhs,
        class ElementMultiplies = std::multiplies<void>>
    auto multiplies(const MatrixLhs& lhs, const MatrixRhs& rhs, ElementMultiplies&& element_mul) {
        return element_op(lhs, rhs, element_mul);
    }
    template<concept_helper::matrix MatrixLhs,
        concept_helper::matrix MatrixRhs,
        concept_helper::matrix MatrixRes = MatrixLhs,
        class ElementDivides = std::divides<void>>
    auto divides(const MatrixLhs& lhs, const MatrixRhs& rhs, ElementDivides&& element_div) {
        return element_op(lhs, rhs, element_div);
    }

    auto operator+(const concept_helper::matrix auto& A, const concept_helper::matrix auto& B) {
        auto res = A;
        foreach_element(res, B,
                [](auto& r, auto& b) {
                    r += b;
                });
        return res;
    }
    auto operator-(const concept_helper::matrix auto& A, const concept_helper::matrix auto& B) {
        auto res = A;
        foreach_element(res, B,
                [](auto& r, auto& b) {
                    r -= b;
                });
        return res;
    }

    template<linear_algebra::concept_helper::matrix Matrix>
    auto operator*(
            Matrix&& A,
            element_type<Matrix> t) {
        auto B = A;
        foreach_element(B, [t](element_type<decltype(B)>& e) {
                e *= t;
                });
        return B;
    }
    template<linear_algebra::concept_helper::matrix Matrix>
    auto operator*(
            element_type<Matrix> t,
            Matrix&& A) {
        return A*t;
    }

    template<linear_algebra::concept_helper::matrix Matrix, typename T>
        requires std::invocable<std::divides<void>, linear_algebra::element_type<Matrix>, T>
    auto matrix_element_divides(
            Matrix&& A,
            T t) {
        auto B = A;
        foreach_element(B, [t](linear_algebra::element_type<decltype(B)>& e) {
                e /= t;
                });
        return B;
    }
    template<linear_algebra::concept_helper::matrix Matrix, typename T>
        requires std::invocable<std::divides<void>, linear_algebra::element_type<Matrix>, T>
    auto operator/(
            Matrix&& A,
            T t) {
        return matrix_element_divides(A, t);
    }
    auto& operator/=(
            linear_algebra::concept_helper::matrix auto&& A,
            auto t) {
        A = A / t;
        return A;
    }

    auto matrix_divides(
            linear_algebra::concept_helper::matrix auto&& A,
            linear_algebra::concept_helper::matrix auto&& B) {
        return A * inverse(B);
    }

    auto operator/(
            linear_algebra::concept_helper::matrix auto&& A,
            linear_algebra::concept_helper::matrix auto&& B) {
        return matrix_divides(A,B);
    }


    auto& operator<<(std::ostream& out, linear_algebra::concept_helper::matrix auto&& A) {
        auto [row_count, column_count] = A.size();
        bool row_is_1 = row_count  == 1;
        out << "{";
        if (!row_is_1) out<< std::endl;
        iterate_from_0_to(row_count,
            [&out, &A, column_count, row_is_1](auto row) {
                out << "{";
                iterate_from_0_to(column_count,
                    [&out, &A, &row](auto column) {
                        out << std::setw(8) << A[{row, column}] << ",";
                    });
                out << "},";
                if (!row_is_1) out<< std::endl;
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

    class inverse{
    public:
        auto operator()(auto v) {
            return static_cast<decltype(v)>(1)/v;
        }
    };
    template<typename T>
        requires std::integral<T> || std::floating_point<T> || complex_type<T>
    auto determinant(T x) {
        return x;
    }
    class is_invertible{
    public:
        auto operator()(auto v) {
            return determinant(v) != decltype(determinant(v)){0};
        }
    };

    constexpr auto debug_matrix_eliminate = false;

    template<concept_helper::matrix Matrix,
        typename Element_subtract = std::minus<void>,
        typename Element_multiplies = std::multiplies<void>,
        typename Element_divides = std::divides<void>,
        typename Element_inverse = inverse,
        typename Element_is_invertible = is_invertible>
    auto& do_eliminate(Matrix& A,
            Element_subtract element_subtract = std::minus<void>{},
            Element_multiplies element_multiplies = std::multiplies<void>{},
            Element_divides element_divides = std::divides<void>{},
            Element_inverse element_inverse = inverse{},
            Element_is_invertible element_is_invertible = is_invertible{}
            ) {
        if (debug_matrix_eliminate) {
            std::cout << A << std::endl;
        }
        auto [row_count, column_count] = A.size();
        for_index_range(0, std::min(row_count, column_count),
            [&A,
             &element_subtract,
             &element_multiplies,
             &element_divides,
             &element_inverse,
             &element_is_invertible](auto column) {
                if (!element_is_invertible(A[{column, column}])) {
                    for (decltype(column) row = column + 1; row < A.size().get_row(); row++) {
                        if (element_is_invertible(A[{row, column}])) {
                            swap(A.row(column), A.row(row));
                            if (debug_matrix_eliminate) {
                                std::cout << "swap rows: " << column << ", " << row << std::endl;
                                std::cout << A << std::endl;
                            }
                            break;
                        }
                    }
                }
                if (element_is_invertible(A[{column, column}])) {
                    for_index_range(column + 1, A.size().get_row(),
                        [&A, &element_subtract, &element_multiplies, &element_divides, column](auto row) {
                            auto pivot = A[{column, column}];
                            auto multier = element_divides(A[{row, column}] , pivot);
                            A.row(row) = 
                                subtract(
                                    A.row(row),
                                    multiplies(A.row(column), multier, {}, element_multiplies),
                                    {},
                                    element_subtract
                                );
                            if (debug_matrix_eliminate) {
                                std::cout << "row " << row << " subtract row " << column << " * " << multier << std::endl;
                                std::cout << A << std::endl;
                            }
                        });
                }
                if (debug_matrix_eliminate) {
                    std::cout << A << std::endl;
                }
            }
        );
        return A;
    }
    template<concept_helper::matrix Matrix>
    Matrix eliminate(Matrix A) {
        return do_eliminate(A);
    }

    constexpr auto debug_matrix_back_substitution = false;
    
    template<concept_helper::matrix Matrix,
        typename Element_subtract = std::minus<void>,
        typename Element_multiplies = std::multiplies<void>,
        typename Element_divides = std::divides<void>,
        typename Element_inverse = inverse,
        typename Element_is_invertible = is_invertible>
    auto& do_back_substitution(Matrix& A,
            Element_subtract element_subtract = std::minus<void>{},
            Element_multiplies element_multiplies = std::multiplies<void>{},
            Element_divides element_divides = std::divides<void>{},
            Element_inverse element_inverse = inverse{},
            Element_is_invertible element_is_invertible = is_invertible{}
            ) {
        for (auto row = A.size().get_row() - 1; true; row--) {
            for (auto column = A.size().get_row() - 1; true; column--) {
                if (column != row) {
                    A.row(row) = subtract(
                            A.row(row),
                            multiplies(
                                A.row(column),
                                A[{row, column}],
                                {},
                                element_multiplies
                                ),
                            {},
                            element_subtract
                            );
                }
                else {
                    auto pivot = A[std::remove_cvref_t<decltype(A.size())>{}.set_row(row).set_column(row)];
                    if (!element_is_invertible(pivot)) {
                        throw std::runtime_error{ std::format("matrix is not invertible") };
                    }
                    auto inv_pivot = element_inverse(pivot);
                    A.row(row) = multiplies(
                            A.row(row),
                            inv_pivot,
                            {},
                            element_multiplies
                            );
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

    template<class T>
    constexpr bool is_approximate_number = std::floating_point<T>;
    template<class T>
    constexpr bool is_approximate_number<std::complex<T>> = std::floating_point<T>;
    template<class T>
    constexpr bool is_exact_number = std::integral<T>;
    template<class T>
    constexpr T accuracy = 1;
    template<>
    constexpr double accuracy<double> = 0.00000001;
    template<>
    constexpr float accuracy<float> = 0.000001f;
    template<class T>
    constexpr T accuracy<std::complex<T>> = accuracy<T>;

    template<concept_helper::matrix Matrix>
    auto gram_schmidt(Matrix&& A) {
        auto res = A;
        auto R = A;
        using element_t = element_type<Matrix>;
        R = I;
        iterate_from_0_to(A.size().get_column(),
            [&res, &R](auto i) {
                iterate_from_0_to(i,
                    [&res, &R, i](auto j) {
                        auto base = res.column(j);
                        auto b = res.column(i);
                        auto scale = dot_product(base, b);
                        res.column(i) -= scale*base;
                        R[{j,i}] = scale;
                    });
                auto base = res.column(i);
                auto len = length(base);
                if (is_approximate_number<element_t>) {
                    if (len <= i*100*accuracy<element_t>) {
                        len = 0;
                        res.column(i) *= 0;
                    }
                    else {
                        res.column(i) /= len;
                    }
                }
                else if (is_exact_number<element_t>) {
                    if (len != 0) res.column(i) /= len;
                }
                else {
                    assert(0);
                }
                R[{i,i}] = len;
            });
        return std::pair{res, R};
    }

    template<typename T>
    struct unit_value_struct {
        static constexpr T value = 1;
    };
    template<typename T>
    constexpr T unit_value = unit_value_struct<T>::value;

    template<concept_helper::matrix Matrix>
    struct unit_value_struct<Matrix> {
        static constexpr Matrix value = I;
    };

    template<concept_helper::matrix Matrix>
    auto determinant(Matrix A) {
        auto U = eliminate(A);
        auto res = unit_value<std::remove_cvref_t<decltype(U[{0, 0}])>>;
        iterate_from_0_to(U.size().get_row(),
            [&res, &U](auto i) {
                res *= U[{i, i}];
            });
        return determinant(res);
    }

    template<concept_helper::matrix Matrix>
    auto eigenvalues(const Matrix& A) {
        if (A.size().get_row() != A.size().get_column()) {
            throw std::runtime_error{"A is not squre matrix, eigenvalues only apply to squre matrix"};
        }
        auto A_i = A;
        for (size_t i = 0; i < 100; i++) {
            auto A_prev = A_i;
            auto [Q, R] = gram_schmidt(A_i);
            A_i = R*Q;

        }
        std::vector<element_type<Matrix>> eigenvalues(A.size().get_row());
        for (size_t j = 0; j < eigenvalues.size(); j++) {
            eigenvalues[j] = A_i[{j,j}];
        }
        return eigenvalues;
    }

    template<concept_helper::matrix Matrix>
    auto eigenvector_matrix_and_eigenvalue_matrix(const Matrix& A) {
        if (A.size().get_row() != A.size().get_column()) {
            throw std::runtime_error{"A is not squre matrix, eigenvalues only apply to squre matrix"};
        }
        auto A_i = A;
        decltype(A_i) X = I;
        for (size_t i = 0; i < 100; i++) {
            auto A_prev = A_i;
            auto [Q, R] = gram_schmidt(A_i);
            A_i = R*Q;
            X = X * Q;
        }
        return std::pair{X, A_i};
    }

    template<concept_helper::matrix Matrix>
    auto singular_value_decompose(const Matrix& A) {
        using std::sqrt;
        auto AT_A = transpose(A) * A;
        auto A_AT = A * transpose(A);
        auto [X, L] = eigenvector_matrix_and_eigenvalue_matrix(AT_A);
        auto [U_X, L_] = eigenvector_matrix_and_eigenvalue_matrix(A_AT);
        decltype(X) V{};
        decltype(U_X) U{};
        std::vector<std::pair<element_type<decltype(L)>, size_t>> singular_values(L.size().get_column());
        for (size_t i = 0; i < singular_values.size(); i++) {
            auto s = sqrt(L[{i,i}]);
            singular_values[i] = {s, i};
        }
        std::ranges::sort(singular_values, [](auto left, auto right) { return length_square(left.first) > length_square(right.first); });
        Matrix S{};
        for (size_t i = 0; i < singular_values.size(); i++) {
            auto [s, index] = singular_values[i];
            if (s == static_cast<decltype(s)>(0)) break;
            S[{i, i}] = s;
            V.column(i) = X.column(index);
            U.column(i) = A*V.column(i)/s;
        }
        return std::tuple{U, S, transpose(V)};
    }
    auto svd(const concept_helper::matrix auto& A) {
        return singular_value_decompose(A);
    }
}

template<linear_algebra::concept_helper::matrix Matrix>
struct std::formatter<Matrix, char> {
};
