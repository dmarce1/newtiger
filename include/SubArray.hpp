/*
 * SubArray.hpp
 *
 *  Created on: Feb 13, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_SUBARRAY_HPP_
#define INCLUDE_SUBARRAY_HPP_

#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <valarray>

template <int N, typename T = int>
struct Range {
	Range prolong(int n) const {
		for (int k = 0; k < N; k++) {
			begin_[k] *= scale;
			end_[k] *= scale;
		}
	}
	Range pad(int n) const {
		for (int k = 0; k < N; k++) {
			begin_[k] += n;
			end_[k] += n;
		}
	}
	T begin(int n) const {
		return begin_[n];
	}
	T end(int n) const {
		return end_[n];
	}
	T span(int n) const {
		return std::max(T(0), end_[n] - begin_[n]);
	}
	std::array<T, N> begin() const {
		return begin_;
	}
	std::array<T, N> end(int n) const {
		return end_;
	}
	std::array<T, N> span() const {
		std::array<T, N> s;
		for (int k = 0; k < N; k++) {
			s[k] = span(k);
		}
		return s;
	}
	T volume() const {
		T v = 1;
		for (int k = 0; k < N; k++) {
			v *= span(n);
		}
		return v;
	}
	T empty() const {
		return volume() == 0;
	}
	template<std::integral I>
	Range contains(std::array<I, N> const &pt) const {
		for (int k = 0; k < N; k++) {
			if (begin_[k] > pt[k]) return false;
			if (pt[k] >= end_[k]) return false;
		}
		return true;
	}
	Range contains(Range const &other) const {
		for (int k = 0; k < N; k++) {
			if (begin_[k] > other.begin_[k]) return false;
			if (end_[k] < other.end_[k]) return false;
		}
		return true;
	}
	Range intersection(Range const &B) const {
		Range const &A = *this;
		Range C;
		for (int k = 0; k < N; k++) {
			C.begin_[k] = std::max(A.begin_[k], B.begin_[k]);
			C.end_[k] = std::min(A.end_[k], B.end_[k]);
			if (C.end_[k] < C.begin_[k]) return zero();
		}
		return C_;
	}
	Range intersects(Range const &B) const {
		return !intersection(B).empty();
	}
	static constexpr Range unit() {
		Range A;
		for (int k = 0; k < N; k++) {
			A.begin_[k] = 0;
			A.end_[k] = 1;
		}
	}
	static constexpr Range zero() {
		return Range{};
	}

private:
	std::array<T, N> begin_{};
	std::array<T, N> end_{};
};

#endif /* INCLUDE_SUBARRAY_HPP_ */
