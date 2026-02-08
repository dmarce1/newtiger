/*
 * Constants.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_


template<typename Type>
inline constexpr Type zero = Type(0);

template<typename Type>
inline constexpr Type one = Type(1);

template<typename Type>
inline constexpr Type two = Type(2);

template<typename Type>
inline constexpr Type half = one<Type> / two<Type>;


#endif /* CONSTANTS_HPP_ */
