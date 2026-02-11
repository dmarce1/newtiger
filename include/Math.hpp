/*
 * Math.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef MATH_HPP_
#define MATH_HPP_

#include <config.hpp>
#include "Constants.hpp"

#include <cmath>
#include <concepts>

template <typename T>
struct IsScalar : std::false_type {};

template <std::integral T>
struct IsScalar<T> : std::true_type {};

template <std::floating_point T>
struct IsScalar<T> : std::true_type {};

template <typename T>
concept Scalar = IsScalar<T>::value;

template <typename T>
constexpr auto sqr(T const &v) {
	return v * v;
}

template <typename T>
constexpr auto inv(T const &v) {
	return 1 / v;
}

template <typename T>
constexpr auto pow(T const &x, std::integral auto n) {
	if (n == 0) {
		return one<T>;
	}
	if (n < 0) {
		return inv(pow(x, -n));
	}
	T z = one<T>;
	T y = x;
	while (n != 1) {
		if ((n & 1) != 0) {
			z *= y;
		}
		n >>= 1;
		y *= y;
	}
	return z * y;
}

template <std::integral T>
constexpr T fact(T n) {
	T y = one<T>;
	for (T l = 1; l <= n; l++) {
		y *= l;
	}
	return y;
}

template <std::integral T>
constexpr T binco(T n, T k) {
	if (n < 2 * k) {
		return binco(n, n - k);
	}
	T y = one<T>;
	T z = one<T>;
	for (T l = 0; l < k; l++) {
		z *= l + 1;
		y *= n - l;
	}
	return y / z;
}

#endif /* MATH_HPP_ */
