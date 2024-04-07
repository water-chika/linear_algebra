#include "matrix.hpp"
#include <iostream>

int main() {
	linear_algebra::vector<float, 2> v{ 4.0f, 1.0f }, u{2.0f, -1.0f };
	auto combination = 2.0f * v + (-3.0f) * u;
	std::cout << combination << std::endl;
	std::cout << dot_product(0.5f*v*2.0f, 2.0f*u*0.5f) << std::endl;
	linear_algebra::matrix<float, 2, 2> A{ linear_algebra::column_vector{v}, linear_algebra::column_vector{u} };
	std::cout << A.column(0) << std::endl;
	std::cout << A.column(1) << std::endl;
	std::cout << A.row(0) << std::endl;
	std::cout << A.row(1) << std::endl;
	A.row_subtract(1, 0, 2);
	std::cout << A.row(1) << std::endl;
	return 0;
}