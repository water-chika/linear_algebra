#pragma once
#include "vector.hpp"
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
			{ m.size() } -> matrix_size;
			m.column(i);
			m.row(i);
		};
	}
	template<class T>
	class matrix_index {
	public:
		constexpr matrix_index() : m_row{}, m_column{} {}
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
	private:
		T m_row;
		T m_column;
	};
	template<class T>
	auto make_matrix_pivot_index(T i) {
		return matrix_index<T>{}.set_row(i).set_column(i);
	}
	template<class T, size_t size>
	class column_vector : public vector<T, size> {
	public:
		column_vector() : vector<T, size>{} {}
		column_vector(vector<T, size> v) : vector<T, size>{ v } {}
	};
	template<class T, size_t Rows, size_t Columns>
	class matrix {
	public:
		static constexpr size_t m = Rows;
		static constexpr size_t n = Columns;
		matrix()
			: m_columns{} {}
		matrix(std::initializer_list<std::initializer_list<T>> row_list)
			: m_columns{ create_columns_from_rows(row_list) } {}
		matrix(std::array<column_vector<T, Rows>, Columns> columns) : m_columns{ columns } {}
		matrix(std::initializer_list<column_vector<T, Rows>> columns) : m_columns{} {
			size_t i = 0;
			for (auto col : columns) {
				m_columns[i] = col;
				i++;
			}
		}
		vector<T, m>& column(size_t i) {
			return m_columns[i];
		}
		vector<T, n> row(size_t i) const {
			assert(i < n);
			vector<T, n> res{};
			for (size_t j = 0; j < n; j++) {
				res[j] = m_columns[j][i];
			}
			return res;
		}
		void row_subtract(size_t i, size_t j, T multi) {
			for (size_t k = 0; k < Columns; k++) {
				m_columns[k][i] -= multi * m_columns[k][j];
			}
		}
		void row_mul(size_t i, T multi) {
			for (size_t k = 0; k < Columns; k++) {
				m_columns[k][i] *= multi;
			}
		}
		static auto create_matrix_by_get_element(auto&& get_element) {
			return matrix{ create_columns_by_get_element(get_element) };
		}
		auto& operator[](concept_helper::matrix_index auto&& i) {
			return m_columns[i.get_column()][i.get_row()];
		}
		static constexpr auto size() {
			return matrix_index<size_t>{}.set_column(Columns).set_row(Rows);
		}
	private:
		static auto create_columns_by_get_element(auto&& get_element) {
			std::array<column_vector<T, m>, n> columns{};
			for (size_t col = 0; col < n; col++) {
				for (size_t row = 0; row < m; row++) {
					columns[col][row] = get_element(row, col);
				}
			}
			return columns;
		}
		static auto create_columns_from_rows(auto&& rows) {
			std::array<column_vector<T, m>, n> columns{};
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
		std::array<column_vector<T, m>, n> m_columns;
	};
	template<class T>
	void iterate_from_0_to(T end, auto&& f) {
		for (T i = 0; i < end; i++) {
			f(i);
		}
	}
	auto& operator<<(std::ostream& out, linear_algebra::concept_helper::matrix auto&& m) {
		out << "{" << std::endl;
		iterate_from_0_to(m.size().get_row(),
			[&out, &m](auto row) {
				out << "{";
				iterate_from_0_to(m.size().get_column(),
					[&out, &m, &row](auto column) {
						out << std::setw(8) << m[linear_algebra::matrix_index<decltype(row)>{}.set_row(row).set_column(column)] << ",";
					});
				out << "}," << std::endl;
			});
		out << "}";
		return out;
	}

	template<class T, size_t Rows, size_t Columns>
	auto make_matrix_with_columns(std::initializer_list<vector<T, Rows>> columns) {
		matrix<T, Rows, Columns> res{};
		for (int i = 0; i < columns.size(); i++) {
			res.column(i) = *(columns.begin() + i);
		}
		return res;
	}
	template<class T, size_t Rows>
	auto identity_matrix() {
		matrix<T, Rows, Rows> res{};
		for (int i = 0; i < Rows; i++) {
			res.column(i)[i] = 1;
		}
		return res;
	}
	auto concatenate_columns(concept_helper::matrix auto&& A, concept_helper::matrix auto&& B) {
		using element_type = std::remove_cvref_t<decltype(A[make_matrix_pivot_index(0)])>;
		matrix<element_type,
			A.size().get_row(), A.size().get_column() + B.size().get_column()> res{};
		static_assert(A.size().get_row() == B.size().get_row());
		for (int i = 0; i < A.size().get_column(); i++) {
			res.column(i) = A.column(i);
		}
		for (int i = A.size().get_column(); i < A.size().get_column() + B.size().get_column(); i++) {
			res.column(i) = B.column(i-A.size().get_column());
		}
		return res;
	}
	template<int... Columns>
	auto select_columns(concept_helper::matrix auto&& A) {
		using element_type = std::remove_cvref_t<decltype(A[make_matrix_pivot_index(0)])>;
		constexpr auto column_indices = std::array{ Columns... };
		matrix<element_type,
			A.size().get_row(), column_indices.size()> res{};
		for (int i = 0; i < column_indices.size(); i++) {
			res.column(i) = A.column(column_indices[i]);
		}
		return res;
	}

	template<class T, size_t m, size_t n, size_t p>
	auto operator*(matrix<T, m, n> lhs, matrix<T, n, p> rhs) {
		return matrix<T, m, p>::create_matrix_by_get_element(
			[&lhs, &rhs](size_t row, size_t col) {
				return dot_product(lhs.row(row), rhs.column(col));
			}
		);
	}

	template<std::integral Int>
	void for_index_range(std::convertible_to<Int> auto&& begin, Int end, std::invocable<Int> auto&& f) {
		for (Int i = begin; i < end; i++) {
			f(i);
		}
	}

	auto eliminate(concept_helper::matrix auto&& A, concept_helper::matrix_index auto&& index) {
		auto column = index.get_column();
		auto row = index.get_row();
		auto pivot = A[make_matrix_pivot_index(column)];
		if (pivot == 0) {
			throw std::runtime_error{ "matrix is not invertible" };
		}
		auto multier = A[index] / pivot;
		auto EA = A;
		EA.row_subtract(row, column, multier);
		return EA;
	}

	template<concept_helper::matrix Matrix>
	Matrix eliminate(Matrix&& A) {
		Matrix res = A;
		for_index_range(0, res.size().get_column(),
			[&res](auto column) {
				for_index_range(column+1, res.size().get_row(),
				[&res, column](auto row) {
						res = eliminate(res, make_matrix_pivot_index(column).set_row(row));
					});
			});
		return res;
	}
	
	template<concept_helper::matrix Matrix>
	Matrix back_substitution(Matrix&& A) {
		Matrix res = A;
		for (auto row = res.size().get_row() - 1; true; row--) {
			for (auto column = res.size().get_row() - 1; true; column--) {
				if (column != row) {
					res.row_subtract(row, column, res[matrix_index<decltype(row)>{}.set_row(row).set_column(column)]);
				}
				else {
					auto pivot = res[make_matrix_pivot_index(row)];
					if (pivot == 0) {
						throw std::runtime_error{ "matrix is not invertible" };
					}
					res.row_mul(row, 1/pivot);
				}
				if (column == row) {
					break;
				}
			}
			if (row == 0) {
				break;
			}
		}
		return res;
	}
}