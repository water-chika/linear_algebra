#pragma once
#include "vector.hpp"

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
	private:
		std::array<column_vector<T, m>, n> m_columns;
	};
}