#pragma once

#include <algorithm>
#include <vector/fixsized_vector.hpp>

namespace linear_algebra {
	
    using fixsized_matrix_index = matrix_index<size_t,size_t>;

	template<class T, size_t ROW, size_t COLUMN>
	class fixsized_matrix {
	public:
        using element_type = T;
		using index_type = fixsized_matrix_index;

		fixsized_matrix()
			: m_columns{} {}
		fixsized_matrix(std::initializer_list<std::initializer_list<T>> row_list)
			: m_columns{ create_columns_from_rows(row_list) } {}
		fixsized_matrix(std::array<column_vector<T, ROW>, COLUMN> columns) : m_columns{ columns } {}
		fixsized_matrix(std::initializer_list<column_vector<T, ROW>> columns) : m_columns{} {
			size_t i = 0;
			for (auto col : columns) {
				m_columns[i] = col;
				i++;
			}
		}
        fixsized_matrix(identity_matrix_type I) : m_columns{} {
            for (size_t i = 0; i < COLUMN; i++) {
                for (size_t j = 0; j < ROW; j++) {
                    m_columns[i][j] = ((i == j) ? 1 : 0);
                }
            }
        }

        auto&& column(this auto&& self, size_t i) {
            return std::forward_like<decltype(self)>(
                        parent_cast<fixsized_matrix&>(self)
                        .m_columns[i]
                    );
        }
		auto row(size_t i) {
			fixsized_pointer_vector<T, COLUMN> res{};
			for (size_t j = 0; j < COLUMN; j++) {
				res.set(j, &m_columns[j][i]);
			}
			return res;
		}
		const auto row(size_t i) const {
			fixsized_pointer_vector<const T, COLUMN> res{};
			for (size_t j = 0; j < COLUMN; j++) {
				res.set(j, &m_columns[j][i]);
			}
			return res;
		}
        auto&& operator[](this auto&& self, fixsized_matrix_index i) {
            return std::forward_like<decltype(self)>(
                    parent_cast<fixsized_matrix&>(self)
                    .m_columns[i.get_column()][i.get_row()]
                    );
        }
		static constexpr auto size() {
			return fixsized_matrix_index{ ROW, COLUMN };
		}
	private:
		static auto create_columns_by_get_element(auto&& get_element) {
			std::array<fixsized_vector<T, ROW>, COLUMN> columns{};
			for (size_t col = 0; col < COLUMN; col++) {
				for (size_t row = 0; row < ROW; row++) {
					columns[col][row] = get_element(row, col);
				}
			}
			return columns;
		}
		static auto create_columns_from_rows(auto&& rows) {
			std::array<fixsized_vector<T, ROW>, COLUMN> columns{};
			size_t row_index = 0;
			for (auto row : rows) {
				size_t col_index = 0;
				for (auto ele : row) {
					columns[col_index][row_index] = ele;
					col_index++;
				}
				row_index++;
			}
			return columns;
		}
		std::array<fixsized_vector<T, ROW>, COLUMN> m_columns;
	};

	template<class T, size_t COLUMN, size_t ROW>
	constexpr size_t fixsized_matrix_column_count<fixsized_matrix<T, COLUMN, ROW>> = COLUMN;

	template<class T, size_t Rows, size_t Columns>
	auto make_fixsized_matrix_with_columns(std::initializer_list<fixsized_vector<T, Rows>> columns) {
		fixsized_matrix<T, Rows, Columns> res{};
		for (int i = 0; i < columns.size(); i++) {
			res.column(i) = *(columns.begin() + i);
		}
		return res;
	}

	template<class T, size_t m, size_t n, size_t p,
        typename Element_add = std::plus<void>,
        typename Element_multiplies = std::multiplies<void>>
    auto multiplies(
            fixsized_matrix<T, m, n> lhs,
            fixsized_matrix<T, n, p> rhs,
            Element_add&& element_add = std::plus<void>{},
            Element_multiplies&& element_multiplies = std::multiplies<void>{}) {
        fixsized_matrix<T, m, p> res;
        foreach_index(res,
                [&res, &lhs, &rhs, &element_add, &element_multiplies](auto index) {
                    res[index] = dot_product(
                            lhs.row(index.get_row()),
                            rhs.column(index.get_column()),
                            element_add,
                            element_multiplies);
                });
        return res;
    }

	template<class T, size_t m, size_t n, size_t p>
	auto operator*(fixsized_matrix<T, m, n> lhs, fixsized_matrix<T, n, p> rhs) {
        fixsized_matrix<T, m, p> res;
        foreach_index(res,
                [&res, &lhs, &rhs](auto index) {
                    res[index] = dot_product(lhs.row(index.get_row()), rhs.column(index.get_column()));
                });
        return res;
	}
	template<class T, size_t m, size_t n>
	auto operator*(const fixsized_matrix<T, m, n>& lhs, const fixsized_vector<T, n> v) {
		fixsized_vector<T, m> res{};
		foreach_index(res,
			[&res, &lhs, &v](auto i) {
				res[i] = dot_product(lhs.row(i), v);
			});
		return res;
	}

    template<class T, size_t m, size_t n>
    auto transpose(const fixsized_matrix<T, m, n>& A) {
        fixsized_matrix<T, n, m> A_T;
        foreach_index(A_T,
                [&A_T, &A](auto i) {
                    A_T[i] = A[{i.get_column(), i.get_row()}];
                });
        return A_T;
    }

	template<typename T, size_t ROW, size_t COLUMN>
	auto remove_column(const fixsized_matrix<T, ROW, COLUMN> A, size_t column) {
		fixsized_matrix<T, ROW, COLUMN - 1> Q;
		iterate_from_0_to(Q.size().get_column(),
			[&Q, &A, column](auto i) {
				Q.column(i) = A.column(i >= column ? i + 1 : i);
			});
		return Q;
	}

	template<typename T, size_t ROW, size_t COLUMN>
	auto cofactor_matrix(const fixsized_matrix<T, ROW, COLUMN> A) {
		auto res = A;
		foreach(A,
			[&A, &res](auto i) {
				auto E = remove_column(remove_row(A, i.get_row()), i.get_column());
				res[i] = ((i.get_column() + i.get_row()) % 2 == 0 ? 1 : -1) * determinant(E);
			});
		return res;
	}
}
