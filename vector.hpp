#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <cassert>

namespace linear_algebra {
	template<class T>
	concept vector = requires (T t1, T t2) {
		t1.size();
		t1[0];
		t1 + t2;
		t1 - t2;
		t1 += t2;
		t1 -= t2;
		t1* t1[0];
		t1[0] * t1;
		t1 / t1[0];
		t1 *= t1[0];
		t1 /= t1[0];
	};
	template<class T, size_t Size>
	class fixsized_vector {
	public:
		fixsized_vector() : m_data{} {}
		fixsized_vector(std::initializer_list<T> data) : m_data{array_from_initializer_list(data)} {}
		T& operator[](size_t i) {
			return m_data[i];
		}
		T operator[](size_t i) const {
			return m_data[i];
		}
		fixsized_vector& operator-=(const fixsized_vector& rhs) {
			for (size_t i = 0; i < Size; i++) {
				m_data[i] -= rhs[i];
			}
			return *this;
		}
	private:
		std::array<T, Size> array_from_initializer_list(std::initializer_list<T> data) {
			assert(data.size() == Size);
			std::array<T, Size> res{};
			std::ranges::copy(data, res.begin());
			return res;
		}
		std::array<T, Size> m_data;
	};
	vector auto operator+(const vector auto& lhs, const vector auto& rhs) {
		auto res{lhs};
		for (size_t i = 0; i < Size; i++) {
			res[i] += rhs[i];
		}
		return res;
	}
	vector auto operator-(const vector auto& lhs, const vector auto& rhs) {
		auto res{ lhs };
		for (size_t i = 0; i < Size; i++) {
			res[i] -= rhs[i];
		}
		return res;
	}
	vector auto dot_product(const vector auto& lhs, const vector auto& rhs) {
		auto res{ lhs };
		for (size_t i = 0; i < Size; i++) {
			res[i] *= rhs[i];
		}
		return res;
	}
	vector auto operator*(const vector auto& lhs, const vector auto& rhs) {
		auto res{ lhs };
		for (size_t i = 0; i < Size; i++) {
			res[i] *= rhs[i];
		}
		return res;
	}
	template<class T, size_t Size>
	fixsized_vector<T, Size> operator*(fixsized_vector<T, Size> lhs, T rhs) {
		fixsized_vector<T, Size> res{};
		for (size_t i = 0; i < Size; i++) {
			res[i] = lhs[i] * rhs;
		}
		return res;
	}
	template<class T, size_t Size>
	fixsized_vector<T, Size> operator/(fixsized_vector<T, Size> lhs, T rhs) {
		fixsized_vector<T, Size> res{};
		for (size_t i = 0; i < Size; i++) {
			res[i] = lhs[i] / rhs;
		}
		return res;
	}
	template<class T, size_t Size>
	std::ostream& operator<<(std::ostream& out, fixsized_vector<T, Size> v) {
		out << "[";
		for (size_t i = 0; i+1 < Size; i++) {
			out << v[i] << ", ";
		}
		out << v[Size - 1];
		out << "]";
		return out;
	}
}
