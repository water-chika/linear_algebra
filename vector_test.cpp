#include <vector.hpp>
#include <vector>
#include <functional>

bool test0() {
    linear_algebra::fixsized_vector<float, 2> a{1.0f, 0.0f};
    a *= 2.0f;
    return a == linear_algebra::fixsized_vector<float, 2>{2.0f, 0.0f};
}

bool test1() {
    linear_algebra::fixsized_vector<double, 3> a{ 0.5, 0.5, 0.5 };
    return 0.75 == dot_product(a, a);
}

bool test_add() {
    linear_algebra::fixsized_vector<double, 2> a{ 0.5, -0.5 };
    linear_algebra::fixsized_vector<double, 2> b{ -1.0, 1.5 };
    return a + b == linear_algebra::fixsized_vector<double, 2>{-0.5, 1.0};
}
bool test_sub() {
    linear_algebra::fixsized_vector<double, 2> a{ 0.5, -0.5 };
    linear_algebra::fixsized_vector<double, 2> b{ -1.0, 1.5 };
    return a - b == linear_algebra::fixsized_vector<double, 2>{1.5, -2.0};
}
bool test_element_multi() {
    linear_algebra::fixsized_vector<double, 2> a{ 0.5, -0.5 };
    linear_algebra::fixsized_vector<double, 2> b{ -1.0, 1.5 };
    return element_multi(a, b) == linear_algebra::fixsized_vector<double, 2>{-0.5, -0.5*1.5};
}

int main() {
    std::vector<std::function<bool()>> tests{
        test0,
        test1,
        test_add,
        test_sub,
        test_element_multi
    };
    for (auto& test : tests) {
        if (test()) {
            std::cout << "passed" << std::endl;
        }
        else {
            std::cout << "failed" << std::endl;
        }
    }
    return 0;
}