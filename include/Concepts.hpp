/*
 * Concepts.hpp
 *
 *  Created on: Feb 17, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_CONCEPTS_HPP_
#define INCLUDE_CONCEPTS_HPP_

#include <concepts>
#include <cstdlib>

template <typename T>
concept ArrayType = requires(T const &t, size_t i) {
	{ t[i] };
};

template <typename T>
concept Scalar = std::is_integral_v<T> || std::is_floating_point_v<T>;

#endif /* INCLUDE_CONCEPTS_HPP_ */
