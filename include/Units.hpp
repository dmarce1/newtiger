/*
 * Units.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef UNITS_HPP_
#define UNITS_HPP_

#include "Integer.hpp"
#include "Matrix.hpp"
#include "Rational.hpp"
#include "Real.hpp"
#include <array>
#include <cmath>

template <Integer L, Integer M, Integer T, Integer K>
struct Units;

template <typename>
struct IsUnits : std::false_type {};

template <Integer L, Integer M, Integer T, Integer K>
struct IsUnits<Units<L, M, T, K>> : std::true_type {};

template <typename T>
concept UnitsType = IsUnits<T>::value;

template <Integer L, Integer M, Integer T, Integer K>
struct Units {
	constexpr inline Units() = default;
	constexpr inline Units(Units const &) = default;
	constexpr inline Units(Units &&) = default;
	constexpr inline Units &operator=(Units const &) = default;
	constexpr inline Units &operator=(Units &&) = default;
	constexpr inline Units &operator+=(Units const &other) {
		v_ += other.v_;
		return *this;
	}
	constexpr inline Units &operator-=(Units const &other) {
		v_ -= other.v_;
		return *this;
	}
	constexpr inline Units &operator*=(Real const &scalar) {
		v_ *= scalar;
		return *this;
	}
	constexpr inline Units &operator/=(Real const &scalar) {
		v_ /= scalar;
		return *this;
	}
	constexpr inline Units const &operator+() const {
		return *this;
	}
	constexpr inline Units operator-() const {
		Units result;
		result.v_ = -v_;
		return result;
	}
	constexpr inline Units operator+(Units const &other) const {
		auto result = *this;
		result += other;
		return result;
	}
	constexpr inline Units operator-(Units const &other) const {
		auto result = *this;
		result -= other;
		return result;
	}
	constexpr inline Units operator*(Real const &scalar) const {
		auto result = *this;
		result *= scalar;
		return result;
	}
	constexpr inline Units operator/(Real const &scalar) const {
		auto result = *this;
		result /= scalar;
		return result;
	}
	template <Integer L2, Integer M2, Integer T2, Integer K2>
	constexpr inline auto operator*(Units<L2, M2, T2, K2> const &other) const {
		constexpr Integer L3 = L + L2;
		constexpr Integer M3 = M + M2;
		constexpr Integer T3 = T + T2;
		constexpr Integer K3 = K + K2;
		auto const res = v_ * other.v_;
		if constexpr ((L3 | M3 | T3 | K3) != 0) {
			Units<L3, M3, T3, K3> result;
			result.v_ = res;
			return result;
		} else {
			return res;
		}
	}
	template <Integer L2, Integer M2, Integer T2, Integer K2>
	constexpr inline auto operator/(Units<L2, M2, T2, K2> const &other) const {
		constexpr Integer L3 = L - L2;
		constexpr Integer M3 = M - M2;
		constexpr Integer T3 = T - T2;
		constexpr Integer K3 = K - K2;
		auto const res = v_ / other.v_;
		if constexpr ((L3 | M3 | T3 | K3) != 0) {
			Units<L3, M3, T3, K3> result;
			result.v_ = res;
			return result;
		} else {
			return res;
		}
	}
	constexpr void set(Real scalar) {
		v_ = scalar;
	}
	constexpr Real get() const {
		return v_;
	}
	static consteval auto powers() {
		Vector<Rational, 4> p;
		p[0] = L;
		p[1] = M;
		p[2] = T;
		p[3] = K;
		return p;
	}
	friend constexpr inline Units operator*(Real const &scalar, Units const units) {
		return units * scalar;
	}
	template <Integer, Integer, Integer, Integer>
	friend class Units;

	Real v_{};
};


#define DEFINE_UNIT(l, m, t, K, name, scale)                                                                                               \
	constexpr auto operator"" _##name(long double value) {                                                                                 \
		Units<l, m, t, K> v;                                                                                                               \
		v.set(value *scale);                                                                                                               \
		return v;                                                                                                                          \
	}                                                                                                                                      \
	constexpr auto operator"" _##name(unsigned long long value) {                                                                          \
		Units<l, m, t, K> v;                                                                                                               \
		v.set(value *scale);                                                                                                               \
		return v;                                                                                                                          \
	}

DEFINE_UNIT(+0, +0, +0, +1, K, 1_R);
DEFINE_UNIT(+0, +0, +1, +0, s, 1_R);
DEFINE_UNIT(+0, +1, +0, +0, g, 1_R);
DEFINE_UNIT(+1, +0, +0, +0, cm, 1_R);
DEFINE_UNIT(+2, +1, -2, +0, erg, 1_R);
DEFINE_UNIT(+0, +0, +0, +4, K4, 1_R);
DEFINE_UNIT(+0, +0, -1, +0, Hz, 1_R);
DEFINE_UNIT(+0, +0, +2, +0, s2, 1_R);
DEFINE_UNIT(+2, +0, +0, +0, cm2, 1_R);
DEFINE_UNIT(+3, +0, +0, +0, cm3, 1_R);
DEFINE_UNIT(+1, +1, -2, +0, dyne, 1_R);

#endif /* UNITS_HPP_ */
