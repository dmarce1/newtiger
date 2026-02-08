/*
 * Units.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef UNITS_HPP_
#define UNITS_HPP_

template <int L, int M, int T, int K, typename S = double>
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
	constexpr inline Units &operator*=(S const &scalar) {
		v_ *= scalar;
		return *this;
	}
	constexpr inline Units &operator/=(S const &scalar) {
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
	constexpr inline Units operator*(S const &scalar) const {
		auto result = *this;
		result *= scalar;
		return result;
	}
	constexpr inline Units operator/(S const &scalar) const {
		auto result = *this;
		result /= scalar;
		return result;
	}
	template <int L2, int M2, int T2, int K2>
	constexpr inline auto operator*(Units<L2, M2, T2, K2, S> const &other) const {
		constexpr int L3 = L + L2;
		constexpr int M3 = M + M2;
		constexpr int T3 = T + T2;
		constexpr int K3 = K + K2;
		auto const res = v_ * other.v_;
		if constexpr ((L3 | M3 | T3 | K3) != 0) {
			Units<L3, M3, T3, K3, S> result;
			result.v_ = res;
			return result;
		} else {
			return res;
		}
	}
	template <int L2, int M2, int T2, int K2>
	constexpr inline auto operator/(Units<L2, M2, T2, K2, S> const &other) const {
		constexpr int L3 = L - L2;
		constexpr int M3 = M - M2;
		constexpr int T3 = T - T2;
		constexpr int K3 = K - K2;
		auto const res = v_ / other.v_;
		if constexpr ((L3 | M3 | T3 | K3) != 0) {
			Units<L3, M3, T3, K3, S> result;
			result.v_ = res;
			return result;
		} else {
			return res;
		}
	}
	constexpr void set(S scalar) {
		v_ = scalar;
	}
	template <int, int, int, int, typename>
	friend class Units;

private:
	S v_{};
};

#define DEFINE_UNIT(l, m, t, K, name, scale)                                                                                               \
	constexpr auto operator"" _##name(long double value) {                                                                                 \
		Units<l, m, t, K, double> v;                                                                                                       \
		v.set(value *scale);                                                                                                               \
		return v;                                                                                                                          \
	}                                                                                                                                      \
	constexpr auto operator"" _##name(unsigned long long value) {                                                                          \
		Units<l, m, t, K, double> v;                                                                                                       \
		v.set(value *scale);                                                                                                               \
		return v;                                                                                                                          \
	}

DEFINE_UNIT(+0, +0, +0, +1, K, 1.0);
DEFINE_UNIT(+0, +0, +0, +4, K4, 1.0);
DEFINE_UNIT(+0, +0, -1, +0, Hz, 1.0);
DEFINE_UNIT(+0, +0, +1, +0, s, 1.0);
DEFINE_UNIT(+0, +0, +2, +0, s2, 1.0);
DEFINE_UNIT(+0, +1, +0, +0, g, 1.0);
DEFINE_UNIT(+1, +0, +0, +0, cm, 1.0);
DEFINE_UNIT(+2, +0, +0, +0, cm2, 1.0);
DEFINE_UNIT(+3, +0, +0, +0, cm3, 1.0);
DEFINE_UNIT(+1, +1, -2, +0, dyne, 1.0);
DEFINE_UNIT(+2, +1, -2, +0, erg, 1.0);

namespace PhysicalConstants {
inline constexpr auto c = 2.99792458e10_cm / 1_s;
inline constexpr auto G = 6.67430e-8_cm3 / (1_g * 1_s2);
inline constexpr auto σ = 5.670374419e-5_erg / (1_cm2 * 1_s * 1_K4);
inline constexpr auto kB = 1.380649e-16_erg / 1_K;
inline constexpr auto NA = 6.02214076e23;
inline constexpr auto π = 3.141592653589793238463;
} // namespace PhysicalConstants


#endif /* UNITS_HPP_ */
