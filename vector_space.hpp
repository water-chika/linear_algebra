#pragma once

#include <vector.hpp>

namespace linear_algebra {
	template<vector Vector>
	class orthonormal_bases_vector_space {
	public:
		orthonormal_bases_vector_space(std::vector<Vector> vectors)
			: bases{
				construct_orthonormal_vectors(vectors)
			} {}
		auto get_bases() {
			return bases;
		}
	private:
		std::vector<Vector> bases;
	};
	template<vector Vector>
	auto project(Vector v, orthonormal_bases_vector_space<Vector> space) {
		auto bases = space.get_bases();
		Vector p{};
		for (auto q : bases) {
			p += q * dot_product(v, q);
		}
		return p;
	}
}