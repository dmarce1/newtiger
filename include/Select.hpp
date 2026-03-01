/*
 * Select.hpp
 *
 *  Created on: Feb 24, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_SELECT_HPP_
#define INCLUDE_SELECT_HPP_

#include "Concepts.hpp"
#include "Simd.hpp"

constexpr auto select(bool f, Scalar auto a, Scalar auto b) {
	return f ? a : b;
}

template<typename T, Integer W>
inline auto select(Simd<bool, W> B, Simd<T, W> const &C, Simd<T, W> const &D) noexcept {
	Simd<T, W> A;
	for (Integer i = 0; i < Simd<bool, W>::size(); i++) {
		A[i] = B[i] ? C[i] : D[i];
	}
	return A;
}

template<typename T, Integer W>
inline auto select(Simd<bool, W> B, T const &c, Simd<T, W> const &D) noexcept {
	Simd<T, W> A;
	for (Integer i = 0; i < Simd<bool, W>::size(); i++) {
		A[i] = B[i] ? c : D[i];
	}
	return A;
}

template<typename T, Integer W>
inline auto select(Simd<bool, W> B, Simd<T, W> const &C, T const &d) noexcept {
	Simd<T, W> A;
	for (Integer i = 0; i < Simd<bool, W>::size(); i++) {
		A[i] = B[i] ? C[i] : d;
	}
	return A;
}

template<Integer W>
inline auto select(Simd<bool, W> f, ArrayType auto const& a, ArrayType auto const&b) {
	using Array = std::remove_cvref_t<decltype(a)>;
	Array c;
	for(Integer w = 0; w < f.size(); w++) {
		for(Integer d = 0; d < a.size(); d++ ) {
			c[d][w] = f[w] ? a[d][w] : b[d][w];
		}
	}
	return c;
}

#endif /* INCLUDE_SELECT_HPP_ */
