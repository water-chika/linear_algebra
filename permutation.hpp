#pragma once

#include <algorithm>

namespace permutation {
    class permutation {
    public:
        permutation(std::vector<size_t> data) : m_data{std::move(data)} {}
        auto size() const {
            return m_data.size();
        }
        auto operator()(size_t x) const {
            assert(x < size());
            return m_data[x];
        }
    private:
        std::vector<size_t> m_data;

        friend permutation operator*(const permutation& lhs, const permutation& rhs);
    };

    permutation operator*(const permutation& lhs, const permutation& rhs){
        assert(lhs.size() == rhs.size());
        std::vector<size_t> data(lhs.size());

        std::vector<size_t> indices(lhs.size());
        std::iota(indices.begin(), indices.end(), 0);
        
        std::for_each(indices.begin(), indices.end(),
                [&data, &lhs, &rhs](auto i) {
                    data[i] = lhs(rhs(i));
                }
                );
        return permutation(std::move(data));
    }
}
