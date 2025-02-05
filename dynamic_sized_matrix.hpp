#pragma once

namespace linear_algebra {
    template<class T>
    class dynamic_sized_matrix {
    public:
        auto column(size_t i) {
        }
        auto row(size_t i) {
        }
        auto operator[](std::pair<size_t,size_t> index) {
        }
        auto size() {
        }
    private:
        std::pari<size_t, size_t> m_size;
        std::vector<T> m_elements;
    };
}
