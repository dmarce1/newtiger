/*
 * Math.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef MATH_HPP_
#define MATH_HPP_

#include "Debug.hpp"
#include "Integer.hpp"
#include "Real.hpp"
#include "Simd.hpp"

#include <cmath>

template <typename T>
constexpr auto sqr(T const &v) {
	return v * v;
}

template <typename T>
constexpr auto inv(T const &v) {
	ASSERT_NONZERO(v);
	return 1_R / v;
}

template <typename T>
constexpr auto pow(T x, Integer n) {
	x = (n < 0_I) ? inv(x) : x;
	n = (n < 0_I) ? -n : +n;
	T z = T(1_I);
	T y = x;
	while (true) {
		z = bool(n & 1_I) ? z * y : z;
		n >>= 1_I;
		if (!n) return z;
		y *= y;
	}
	return z;
}

constexpr Integer fact(Integer n) {
	Integer y = 1_I;
	for (Integer l = 1_I; l <= n; l++) {
		y *= l;
	}
	return y;
}

constexpr Integer binco(Integer n, Integer k) {
	if (n < 2_I * k) {
		return binco(n, n - k);
	}
	Integer y = 1_I;
	Integer z = 1_I;
	for (Integer l = 0_I; l < k; l++) {
		z *= l + 1_I;
		y *= n - l;
	}
	return y / z;
}

constexpr Real root(Real x, Integer n) {
	int exp;
	Real y, y0, a;
	std::frexp(x, &exp);
	y = 1_I << Integer(std::abs(exp) / Real(n) + 0.5_R);
	if (x < 1_R) y = inv(y);
	a = inv(n);
	do {
		y0 = y;
		y *= 1_I + a * (x * pow(y, -n) - 1_I);
	} while (std::abs(y - y0) > n * std::abs(y * eps_R));
	return y;
}


#endif /* MATH_HPP_ */
