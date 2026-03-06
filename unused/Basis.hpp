/*
 * Basis.cpp
 *
 *  Created on: Mar 5, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_BASIS_CPP_
#define INCLUDE_BASIS_CPP_

#include "Matrix.hpp"
#include "Quadrature.hpp"

#include <vector>

template <typename T, Integer O>
auto basis() {
	constexpr auto tmp = gaussLobattoQuadrature<O + 1>();
	constexpr auto x = tmp.first;
	constexpr auto w = tmp.second;
	constexpr Integer M = O;
	constexpr Integer N = O + 1;
	Matrix<T, M, N> n2m;
	Matrix<T, N, M> m2n;
	for (Integer n = 0; n < N; n++) {
		for (Integer m = 0; m < M; m++) {
			auto const p = legendreP(n, x[m]);
			m2n[n][m] = p;
			n2m[m][n] = p * w[m] * 0.5_R * Real(2_I * m + 1_I);
		}
	}
	return std::pair(m2n, n2m);
}

template <typename T, Integer O, Integer D>
void legendreTransform(auto &dst, auto const &src) {
	constexpr Integer M = O;
	constexpr Integer N = O + 1;
	auto const modeCount = binco(M + D - 1, D);
	auto const nodeCount = pow(N, D);
	auto const lambda = [&]<Integer n0>() {
		constexpr auto m2n = basis<T, O - n0>().first;
		auto const L = pow(N, D - 1);
		for (Integer n = 0; n0 + n < N; n++) {
			for (Integer l = 0; l < L; l++) {
				dst[n * L + l] = 0_R;
			}
			for (Integer m = 0; m < M; m++) {
				for (Integer l = 0; l < L; l++) {
					dst[n * L + l] += m2n[n][m] * src[m * L + l];
				}
			}
		}
	};
}

#endif /* INCLUDE_BASIS_CPP_ */
