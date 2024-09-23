#include "matrix.hpp"
#include "combined_reference_matrix.hpp"
#include <iostream>
#include <measure_duration.hpp>
#include <concepts>
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
                linear_algebra::fixsized_vector<Number, 2> v{ 4.0f, 1.0f }, u{ 2.0f, -1.0f };
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
                
                auto U = linear_algebra::eliminate(A);
                std::cout << "The elimination U of A is below: " << U << std::endl;

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

                {
                    using namespace linear_algebra;
                    auto A_I = make_matrix_with_columns<double, 3, 6>({fixsized_vector<double, 3>{2,1,1}, fixsized_vector<double, 3>{1,2,1}, fixsized_vector<double, 3>{1,1,2},
                        fixsized_vector<double, 3>{1,0,0}, fixsized_vector<double, 3>{0,1,0}, fixsized_vector<double, 3>{0,0,1}});
                    std::cout << "set A_I = " << A_I << std::endl;
                    auto res = eliminate(A_I);
                    std::cout << "after eliminate, it =" << res << std::endl;
                    res = back_substitution(res);
                    std::cout << "after back substitution, it =" << res << std::endl;
                }
                {
                    auto A = linear_algebra::matrix<int, 3, 3>{
                        {1, 0, 1},
                        {2, 3, 0},
                        {5, 4, 7}
                    };
                    std::cout << "test initializer, A = " << A << std::endl;
                }
                try{
                    using namespace linear_algebra;
                    std::cout << "This is chapter 2, exercise 31." << std::endl;
                    auto A = matrix<Number, 3, 3>{
                        {2, 1, 1},
                        {1, 2, 1},
                        {1, 1, 2}
                    };
                    std::cout << "A = " << A << std::endl;
                    auto A_i = identity_matrix<Number, 3>();
                    auto A_I = combined_reference_matrix{A, A_i};
                    back_substitution(eliminate(A_I));
                    std::cout << "A inverse = " << A_i << std::endl;
                    auto B = matrix<Number, 3, 3>{
                        {2, -1, -1},
                        {-1, 2, -1},
                        {-1, -1, 2}
                    };
                    std::cout << "B = " << B << std::endl;
                    auto B_I = concatenate_columns(B, identity_matrix<Number, 3>());
                    std::cout << "B inverse = " <<
                        select_columns<3, 4, 5>(back_substitution(eliminate(B_I))) << std::endl;
                }
                catch (std::exception& e) {
                    std::cout << e.what() << std::endl;
                }
            }
        };
    };
};

int main() {
    linear_algebra_test::set_number<double>::set_index<uint8_t>::run();
    return 0;
}
