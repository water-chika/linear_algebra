#include <matrix.hpp>

int main() {
	auto A = linear_algebra::matrix<float, 3, 3>{
		{1, 2, 3},
		{-1, 0, -3},
		{0, -2, 3}
	};
	std::cout << A << std::endl;
	std::cout << "After gram schmidt process:" << std::endl;
	std::cout << gram_schmidt(A) << std::endl;
	return 0;
}