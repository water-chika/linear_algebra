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
	};
	template<class Vector, class ElementType>
	concept vector_element_type = vector<Vector> && requires(Vector v) {
		{ v[0] } -> std::convertible_to<ElementType>;
	};
	template<class T, size_t Size>
	class fixsized_vector {
	public:
		fixsized_vector() : m_data{} {}
		fixsized_vector(std::initializer_list<T> data) : m_data{array_from_initializer_list(data)} {}
		size_t size() const {
			return Size;
		}
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
		fixsized_vector& operator+=(const fixsized_vector& rhs) {
			for (size_t i = 0; i < Size; i++) {
				m_data[i] += rhs[i];
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
	template<vector Vector, class T>
		requires vector_element_type<Vector, T>
	Vector operator*(const Vector& lhs, const T& rhs) {
		Vector res{ lhs };
		for (size_t i = 0; i < res.size(); i++) {
			res[i] *= rhs;
		}
		return res;
	}
	template<vector Vector, class T>
		requires vector_element_type<Vector, T>
	Vector& operator*=(Vector& lhs, const T& rhs) {
		lhs = lhs * rhs;
		return lhs;
	}
	bool operator==(const vector auto& lhs, const vector auto& rhs) {
		bool equal{ lhs.size() == rhs.size() };
		for (size_t i = 0; equal && i < lhs.size(); i++) {
			equal = equal && (lhs[i] == rhs[i]);
		}
		return equal;
	}
	vector auto operator+(const vector auto& lhs, const vector auto& rhs) {
		if (lhs.size() != rhs.size()) {
			throw std::runtime_error{ "size not equal" };
		}
		auto res{lhs};
		for (size_t i = 0; i < lhs.size(); i++) {
			res[i] += rhs[i];
		}
		return res;
	}
	vector auto operator-(const vector auto& lhs, const vector auto& rhs) {
		if (lhs.size() != rhs.size()) {
			throw std::runtime_error{ "size not equal" };
		}
		auto res{ lhs };
		for (size_t i = 0; i < res.size(); i++) {
			res[i] -= rhs[i];
		}
		return res;
	}
	auto dot_product(const vector auto& lhs, const vector auto& rhs) {
		if (lhs.size() != rhs.size()) {
			throw std::runtime_error{ "size not equal" };
		}
		auto res{ lhs[0]*rhs[0] };
		for (size_t i = 1; i < lhs.size(); i++) {
			res += lhs[i] * rhs[i];
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
