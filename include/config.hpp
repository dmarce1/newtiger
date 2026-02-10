#pragma once 

#include "Vector.hpp"

namespace config {
	using Real = double;
	
	template<typename T, int N>
	using Simd = Vector<T, N>;
	
	constexpr int dimCount = 1;
}

