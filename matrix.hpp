#pragma once
#include "vector.hpp"
#include <algorithm>

namespace linear_algebra {
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
}