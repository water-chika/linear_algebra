#include "matrix.hpp"
#include <iostream>
#include <measure_duration.hpp>
#include <concepts>

template<class T>
void iterate_from_0_to(T end, auto&& f) {
	for (T i = 0; i < end; i++) {
		f(i);
	}
}

auto& operator<<(std::ostream& out, linear_algebra::concept_helper::matrix auto&& m) {
	out << std::endl;
	iterate_from_0_to(m.size().get_row(),
		[&out, &m](auto row) {
			out << "|";
			iterate_from_0_to(m.size().get_column(),
			[&out, &m, &row](auto column) {
					out << std::setw(4) << m[linear_algebra::matrix_index<decltype(row)>{}.set_row(row).set_column(column)];
				});
			out << "|" << std::endl;
		});
	return out;
}

namespace water {
	namespace concept_helper {
		template<class T>
		concept number = std::integral<T> || std::floating_point<T>;
	}
}

class linear_algebra_test {
public:
	template<water::concept_helper::number Number>
	class set_number {
	public:
		template<std::integral Integral>
		class set_index {
		public:
			static void run() {
				linear_algebra::vector<Number, 2> v{ 4.0f, 1.0f }, u{ 2.0f, -1.0f };
				auto combination = Number{ 2.0 } *v + Number{ -3.0f }*u;
				std::cout << combination << std::endl;
				std::cout << dot_product(Number{ 0.5f } *v * Number{ 2.0f }, Number{ 2.0f } *u * Number{ 0.5f }) << std::endl;
				linear_algebra::matrix<Number, 2, 2> A{ linear_algebra::column_vector{v}, linear_algebra::column_vector{u} };
				std::cout << "A is below: " << A << std::endl;
				std::cout << A.column(0) << std::endl;
				std::cout << A.column(1) << std::endl;
				std::cout << A.row(0) << std::endl;
				std::cout << A.row(1) << std::endl;
				A.row_subtract(1, 0, 2);
				std::cout << "After subtract 2 * row 0 from row 1, A is below: " << A << std::endl;

				auto large_matrix_100_101 = linear_algebra::matrix<Number, 100, 101>::create_matrix_by_get_element(
					[](Integral row, Integral col) {
						return row + col;
					}
				);
				auto large_matrix_101_100 = linear_algebra::matrix<Number, 101, 100>::create_matrix_by_get_element(
					[](Integral row, Integral col) {
						return row + col + 0.1f;
					}
				);
				auto duration = water::measure_duration_average<100>(
					[&large_matrix_100_101, &large_matrix_101_100]() {
						large_matrix_100_101* large_matrix_101_100;
					});
				std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(duration) << std::endl;
			}
		};
	};
};

int main() {
	linear_algebra_test::set_number<double>::set_index<uint8_t>::run();
	return 0;
}