/*
 * DoubleValue.hpp
 *
 *  Created on: Feb 22, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_DOUBLEVALUE_HPP_
#define INCLUDE_DOUBLEVALUE_HPP_

#include "Real.hpp"
#include "Vector.hpp"

template <typename T = Real>
struct DoubleValue : Vector<2, T> {
	DoubleValue() :
		L(*(new(&(*this)(0)) T{})), R(*(new(&(*this)(1)) T{})) {
	}
	virtual ~DoubleValue() {
		L.~T();
		R.~T();
	}
	T &L;
	T &R;
};

#endif /* INCLUDE_DOUBLEVALUE_HPP_ */
