/*
 * Limiters.hpp
 *
 *  Created on: Feb 24, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_LIMITERS_HPP_
#define INCLUDE_LIMITERS_HPP_

#include <cmath>

#include "Real.hpp"
#include "Select.hpp"

inline constexpr Real minmod(Real a, Real b) {
	return (std::copysign(0.5_R, a) + std::copysign(0.5_R, b)) * std::min(std::abs(a), std::abs(b));
}

inline constexpr auto minmod(ArrayType auto A, ArrayType auto const& B) {
	constexpr auto size = A.size();
	for(Integer i = 0; i < size; i++) {
		A[i] = minmod(A[i], B[i]);
	}
	return A;
}

#endif /* INCLUDE_LIMITERS_HPP_ */
