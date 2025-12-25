#pragma once
#include <matrix/matrix.hpp>
#include <vector/resizeable_vector.hpp>
#include "type_cast.hpp"

namespace linear_algebra {
    using dynamic_sized_matrix_index = matrix_index<size_t,size_t>;
    template<class T>
    class dynamic_sized_matrix {
    public:
        using element_type = T;
        using index_type = dynamic_sized_matrix_index;
        dynamic_sized_matrix() = default;
        dynamic_sized_matrix(index_type s)
            : m_size{s}, m_elements(m_size.get_column() * m_size.get_row())
        {}
        template<matrix Matrix>
            requires (!std::same_as<std::remove_cvref_t<Matrix>, dynamic_sized_matrix>)
        dynamic_sized_matrix(Matrix& A)
            : m_size{A.size()}, m_elements(m_size.get_column() * m_size.get_row())
        {
            foreach_index(m_size,
                    [this, &A](auto i) {
                        (*this)[i] = A[i];
                    }
                    );
        }
        class column_ref {
        public:
            using element_type = T;
            column_ref(dynamic_sized_matrix* m, size_t i)
                : m_matrix{m}, m_column_index{i}
            {}
            auto size() const {
                auto [m, n] = m_matrix->size();
                return m;
            }
            auto&& operator[](this auto&& self, size_t i) {
                auto& ref = parent_cast<column_ref&>(self);
                return std::forward_like<decltype(self)>(
                        ref.m_matrix->operator[]({i, ref.m_column_index})
                        );
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
            using element_type = T;
            row_ref(dynamic_sized_matrix* m, size_t i)
                : m_matrix{m}, m_row_index{i}
            {}
            auto size() const {
                auto [m, n] = m_matrix->size();
                return n;
            }
            auto&& operator[](this auto&& self, size_t i) {
                auto& ref = parent_cast<row_ref&>(self);
                return std::forward_like<decltype(self)>(
                        ref.m_matrix->operator[]({ref.m_row_index, i})
                        );
            }
        private:
            dynamic_sized_matrix* m_matrix;
            size_t m_row_index;
        };
        auto row(size_t i) {
            return reference_vector{row_ref{this, i}};
        }
        auto&& operator[](this auto&& child, index_type index) {
            auto& self = parent_cast<dynamic_sized_matrix&>(child);
            assert(index.get_column() < self.m_size.get_column());
            assert(index.get_row() < self.m_size.get_row());
            return self.m_elements[index.get_column()*self.m_size.get_row() + index.get_row()];
        }
        auto size() const {
            return m_size;
        }
        auto resize(index_type s) {
            m_size = s;
            m_elements.resize(s.get_row()*s.get_column());
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
}
