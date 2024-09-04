#pragma once

namespace linear_algebra {
	template<class T, size_t Size>
	class fixsized_vector {
	public:
		fixsized_vector() : m_data{} {}
		fixsized_vector(std::initializer_list<T> data) : m_data{ array_from_initializer_list(data) } {}
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
}