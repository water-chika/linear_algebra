#pragma once
#include "vector.hpp"
#include <algorithm>

namespace linear_algebra {
	template<class T, size_t size>
	class column_vector : public vector<T, size> {
	public:
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
		vector<T, m> row(size_t i) const {
			assert(i < m);
			vector<T, m> res{};
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
	private:
		std::array<column_vector<T, m>, n> m_columns;
	};
}