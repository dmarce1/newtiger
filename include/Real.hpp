/*
 * Real.hpp
 *
 *  Created on: Feb 20, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_REAL_HPP_
#define INCLUDE_REAL_HPP_

#include <cmath>
#include <limits>

#include "Integer.hpp"

using Real = double;

inline constexpr Real operator"" _R(long double x) {
	return Real(x);
}

inline constexpr Real operator"" _R(long long unsigned x) {
	return Real(x);
}

inline constexpr auto eps_R = std::numeric_limits<Real>::epsilon();
inline constexpr auto huge_R = std::numeric_limits<Real>::max();
inline constexpr auto inf_R = std::numeric_limits<Real>::infinity();
inline constexpr auto NaN_R = std::numeric_limits<Real>::signaling_NaN();
inline constexpr auto tiny_R = std::numeric_limits<Real>::min();

inline constexpr Real inv(Real v) {
	return 1_R / v;
}

inline constexpr Real sqr(Real v) {
	return v * v;
}

inline constexpr Integer pow(Real x, Integer n) {
	if (n == 0_I) return 1_R;
	if (n < 0_I) return inv(pow(x, -n));
	Real z = 1_R;
	Real y = x;
	while (n != 1_I) {
		if ((n & 1_I) != 0_I) {
			z *= y;
		}
		n >>= 1_I;
		y *= y;
	}
	return z * y;
}

inline constexpr Real sgn(Real x) {
	return std::copysign(1_R, x);
}

inline constexpr Real minmod(Real a, Real b) {
	return 0.5_R * (sgn(a) + sgn(b)) * std::min(std::abs(a), std::abs(b));
}

inline constexpr Real minmod(Real a, Real b, Real θ) {
	return minmod(θ * minmod(a, b), 0.5_R * (a + b));
}


#endif /* INCLUDE_REAL_HPP_ */
