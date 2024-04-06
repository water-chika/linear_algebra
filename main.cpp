#include "matrix.hpp"
#include "vector.hpp"
#include <iostream>

int main() {
	linear_algebra::vector<float, 2> v{ 0.4f, 1.0f }, u{ 0.1f, 1.2f };
	std::cout << dot_product(v, u) << std::endl;
	return 0;
}