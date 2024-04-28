#pragma once
#include "vector.hpp"
#include <algorithm>

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
		auto set_row(T row) {
			m_row = row;
			return *this;
		}
		auto set_column(T column) {
			m_column = column;
			return *this;
		}
		auto get_row() const {
			return m_row;
		}
		auto get_column() const {
			return m_column;
		}
	private:
		T m_row;
		T m_column;
	};
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
		template<class... List>
		matrix(List... list) : m_columns{ list... } {}
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
		static auto create_matrix_by_get_element(auto&& get_element) {
			return matrix{ create_columns_by_get_element(get_element) };
		}
		auto& operator[](concept_helper::matrix_index auto&& i) {
			return m_columns[i.get_column()][i.get_row()];
		}
		auto size() const {
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
		std::array<column_vector<T, m>, n> m_columns;
	};

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

	template<concept_helper::matrix Matrix>
	Matrix eliminate(Matrix&& A) {
		Matrix res = A;
		for_index_range(0, res.size().get_column(),
			[&res](auto column) {
				for_index_range(column+1, res.size().get_row(),
				[&res, column](auto row) {
						auto pivot = res[matrix_index<decltype(column)>{}.set_column(column).set_row(column)];
						assert(pivot != 0);
						auto multi = res[matrix_index<decltype(column)>{}.set_column(column).set_row(row)] / pivot;
						res.row_subtract(row, column, multi);
					});
			});
		return res;
	}
}