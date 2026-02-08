/*
 * Legendre.hpp
 *
 *  Created on: Feb 5, 2026
 *      Author: dmarce1
 */

#ifndef LEGENDRE_HPP_
#define LEGENDRE_HPP_

#include "Constants.hpp"

#include <array>
#include <cassert>
#include <type_traits>

constexpr auto legendreP(int l, std::floating_point auto x) {
	using T = decltype(x);
	if (l == 0) {
		return one<T>;
	}
	T Pn = one<T>;
	T Pnm1 = zero<T>;
	T Pnp1;
	for (int n = 0;; n++) {
		Pnp1 = (2 * n + 1) * x * Pn - n * Pnm1;
		Pnp1 /= n + 1;
		if (n + 1 == l) {
			return Pnp1;
		}
		Pnm1 = Pn;
		Pn = Pnp1;
	}
	assert(false);
	return T{};
};

template<int D = 1>
constexpr auto dLegendrePdx(int l, std::floating_point auto x) {
	using T = decltype(x);
	if (l == 0) {
		return zero<T>;
	}
	std::array<T, D + 1> Pn{};
	std::array<T, D + 1> Pnm1{};
	std::array<T, D + 1> Pnp1;
	Pn[0] = one<T>;
	for (int n = 0;; n++) {
		int const b = std::max(D + 1 - l + n, 0);
		for (int m = b; m <= D; m++) {
			Pnp1[m] = x * Pn[m];
			if (m) Pnp1[m] += m * Pn[m - 1];
			Pnp1[m] *= 2 * n + 1;
			Pnp1[m] -= n * Pnm1[m];
			Pnp1[m] /= n + 1;
		}
		if (n + 1 == l) {
			return Pnp1[D];
		}
		Pnm1 = Pn;
		Pn = Pnp1;
	}
	assert(false);
	return T{};
};

#endif /* LEGENDRE_HPP_ */
