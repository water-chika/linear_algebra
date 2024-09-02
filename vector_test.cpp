#include <vector.hpp>
#include <vector>
#include <functional>

bool test0() {
    linear_algebra::fixsized_vector<float, 2> a{1.0f, 0.0f};
    a *= 2.0f;
    return a == linear_algebra::fixsized_vector<float, 2>{2.0f, 0.0f};
}

int main() {
    std::vector<std::function<bool()>> tests{
        test0
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