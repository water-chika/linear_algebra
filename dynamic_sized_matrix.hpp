#pragma once
#include <matrix.hpp>

namespace linear_algebra {
    using dynamic_sized_matrix_index = matrix_index<size_t,size_t>;
    template<class T>
    class dynamic_sized_matrix {
    public:
        using index_type = dynamic_sized_matrix_index;
        dynamic_sized_matrix() = default;
        dynamic_sized_matrix(dynamic_sized_matrix_index s)
            : m_size{s}, m_elements(m_size.get_column() * m_size.get_row())
        {}
        class column_ref {
        public:
            column_ref(dynamic_sized_matrix* m, size_t i)
                : m_matrix{m}, m_column_index{i}
            {}
            auto size() {
                auto [m, n] = m_matrix->size();
                return m;
            }
            auto& operator[](size_t i) {
                return m_matrix->operator[]({i,m_column_index});
            }
            const auto& operator[](size_t i) const {
                return m_matrix->operator[]({i,m_column_index});
            }
        private:
            dynamic_sized_matrix* m_matrix;
            size_t m_column_index;
        };
        auto column(size_t i) {
            return reference_vector{column_ref{this, i}};
        }
        class row_ref {
        public:
            row_ref(dynamic_sized_matrix* m, size_t i)
                : m_matrix{m}, m_row_index{i}
            {}
            auto size() {
                auto [m, n] = m_matrix->size();
                return n;
            }
            auto& operator[](size_t i) {
                return m_matrix->operator[]({m_row_index,i});
            }
            const auto& operator[](size_t i) const {
                return m_matrix->operator[]({m_row_index,i});
            }
        private:
            dynamic_sized_matrix* m_matrix;
            size_t m_row_index;
        };
        auto row(size_t i) {
            return reference_vector{row_ref{this, i}};
        }
        auto& operator[](index_type index) {
            assert(index.get_column() < m_size.get_column());
            assert(index.get_row() < m_size.get_row());
            return m_elements[index.get_column()*m_size.get_row() + index.get_row()];
        }
        const auto& operator[](index_type index) const {
            assert(index.get_column() < m_size.get_column());
            assert(index.get_row() < m_size.get_row());
            return m_elements[index.get_column()*m_size.get_row() + index.get_row()];
        }
        auto size() {
            return m_size;
        }

        auto& operator=(concept_helper::matrix auto&& A) {
            assert(m_size.get_column() == A.size().get_column());
            assert(m_size.get_row() == A.size().get_row());
            foreach_index(m_size,
                    [this, &A](auto i) {
                        (*this)[i] = A[i];
                    });
            return *this;
        }
        auto& operator=(identity_matrix_type I) {
            foreach_index(m_size,
                    [this](auto i) {
                        (*this)[i] = i.get_column() == i.get_row() ? 1 : 0;
                    });
            return *this;
        }
    private:
        index_type m_size;
        std::vector<T> m_elements;
    };

    template<class T, class U>
    auto transpose(matrix_index<T, U> size) {
        auto transposed_size = size;
        transposed_size.set_column(size.get_row());
        transposed_size.set_row(size.get_column());
        return transposed_size;
    }

    template<class T>
    auto transpose(dynamic_sized_matrix<T> A) {
        auto size = A.size();
        auto transposed_size = transpose(size);
        dynamic_sized_matrix<T> AT(transposed_size);
        foreach_index(transposed_size,
                [&AT, &A](auto i) {
                    AT[i] = A[transpose(i)];
                });
        return AT;
    }
    template<class T>
    auto operator*(dynamic_sized_matrix<T> lhs, dynamic_sized_matrix<T> rhs) {
        assert(lhs.size().get_column() == rhs.size().get_row());
        auto size = lhs.size();
        size.set_column(rhs.size().get_column());
        dynamic_sized_matrix<T> res(size);
        foreach_index(res,
                [&res, &lhs, &rhs](auto index) {
                    res[index] = dot_product(lhs.row(index.get_row()), rhs.column(index.get_column()));
                });
        return res;
    }
}
