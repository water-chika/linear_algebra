#pragma once
#include <vector>
#include <array>
#include <iostream>

namespace linear_algebra {
	template<class T, size_t Size>
	class vector {
	public:
		vector() : m_data{} {}
		template<class... List>
		vector(List... data) : m_data{ data... } {}
		T& operator[](size_t i) {
			return m_data[i];
		}
	private:
		std::array<T, Size> m_data;
	};
	template<class T, size_t Size>
	vector<T, Size> operator+(vector<T, Size> lhs, vector<T, Size> rhs) {
		vector<T, Size> res{};
		for (size_t i = 0; i < Size; i++) {
			res[i] = lhs[i] + rhs[i];
		}
		return res;
	}
	template<class T, size_t Size>
	vector<T, Size> operator-(vector<T, Size> lhs, vector<T, Size> rhs) {
		vector<T, Size> res{};
		for (size_t i = 0; i < Size; i++) {
			res[i] = lhs[i] - rhs[i];
		}
		return res;
	}
	template<class T, size_t Size>
	T dot_product(vector<T, Size> lhs, vector<T, Size> rhs) {
		T res{};
		for (size_t i = 0; i < Size; i++) {
			res += lhs[i] * rhs[i];
		}
		return res;
	}
	template<class T, size_t Size>
	vector<T, Size> operator*(T lhs, vector<T, Size> rhs) {
		vector<T, Size> res{};
		for (size_t i = 0; i < Size; i++) {
			res[i] = lhs * rhs[i];
		}
		return res;
	}
	template<class T, size_t Size>
	vector<T, Size> operator*(vector<T, Size> lhs, T rhs) {
		vector<T, Size> res{};
		for (size_t i = 0; i < Size; i++) {
			res[i] = lhs[i] * rhs;
		}
		return res;
	}
	template<class T, size_t Size>
	std::ostream& operator<<(std::ostream& out, vector<T, Size> v) {
		out << "[";
		for (size_t i = 0; i+1 < Size; i++) {
			out << v[i] << ", ";
		}
		out << v[Size - 1];
		out << "]";
		return out;
	}
}
