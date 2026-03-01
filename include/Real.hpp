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




#endif /* INCLUDE_REAL_HPP_ */
