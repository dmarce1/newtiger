/*
 * Quadrature.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef QUADRATURE_HPP_
#define QUADRATURE_HPP_

#include "Constants.hpp"
#include "IO.hpp"
#include "Legendre.hpp"
#include "Math.hpp"
#include "Matrix.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <numeric>
#include <ostream>

template <typename T, int N>
struct GaussLegendreQuadrature {
	static constexpr int size() {
		return N;
	}
	static constexpr auto const exact2degree() {
		return 2 * N - 1;
	}
	constexpr GaussLegendreQuadrature() {
		using HiPrec = long double;
		constexpr auto π = std::numbers::pi_v<HiPrec>;
		constexpr int No2 = N / 2;
		for (int k = 0; k < No2; k++) {
			HiPrec x, p0;
			HiPrec xm = cos(π * HiPrec(2 * k + 2) / HiPrec(2 * N + 1));
			HiPrec xp = cos(π * HiPrec(2 * k + 1) / HiPrec(2 * N + 1));
			HiPrec pp = legendreP(N, xp);
			do {
				x = half<HiPrec> * (xm + xp);
				p0 = legendreP(N, x);
				if ((pp * p0) > zero<HiPrec>) {
					xp = x;
					pp = p0;
				} else {
					xm = x;
				}
			} while (T(xp) > T(xm));
			x_[k] = +T(x);
			x_[N - k - 1] = -T(x);
			w_[k] = w_[N - k - 1] = two<HiPrec> * inv((one<HiPrec> - sqr(x)) * sqr(dLegendrePdx(N, x)));
		}
		if constexpr (N % 2 == 1) {
			x_[N / 2] = zero<T>;
			w_[N / 2] = two<HiPrec> * inv(sqr(dLegendrePdx(N, zero<T>)));
		}
	}
	constexpr T weight(int n) const {
		return w_[n];
	}
	constexpr T position(int n) const {
		return x_[n];
	}
	constexpr auto const &positions() const {
		return x_;
	}
	constexpr auto const &weights() const {
		return w_;
	}

private:
	Vector<T, size()> x_{};
	DiagonalMatrix<T, size()> w_{};
};

// deg = order - 1
// 2*deg = 2*order - 2 <= 2 * n - 1
// deg = order - 1 <= n - 0.5
// 2*deg = 2*order - 2 <= 2 * n - 3
// order + 1 <= n
// order - 0.5 <= n
// 2 * O + 2 = 2 * N - 1
// 2 * O + 1 = 2 * N
// O + 0.5 =  N

template <typename T, int N>
struct GaussLobattoQuadrature {
	static constexpr int size() {
		return N;
	}
	static constexpr auto const exact2degree() {
		return 2 * N - 3;
	}
	constexpr GaussLobattoQuadrature() {
		using HiPrec = long double;
		constexpr int No2 = N / 2;
		constexpr int Nm1 = N - 1;
		constexpr GaussLegendreQuadrature<T, Nm1> gl{};
		for (int k = 1; k < No2; k++) {
			HiPrec x, p0;
			HiPrec xm = gl.position(k - 1);
			HiPrec xp = gl.position(k);
			if (xp < xm) std::swap(xp, xm);
			HiPrec pp = dLegendrePdx<1>(Nm1, xp);
			do {
				x = half<HiPrec> * (xm + xp);
				p0 = dLegendrePdx<1>(Nm1, x);
				if ((pp * p0) > zero<HiPrec>) {
					xp = x;
					pp = p0;
				} else {
					xm = x;
				}
			} while (T(xp) > T(xm));
			x_[k] = +T(x);
			x_[N - k - 1] = -T(x);
			w_[k] = w_[N - k - 1] = two<HiPrec> * inv(HiPrec(N * Nm1) * sqr(legendreP(Nm1, x)));
		}
		x_[0] = +T(1);
		x_[Nm1] = -T(1);
		w_[0] = w_[Nm1] = two<HiPrec> * inv(HiPrec(N * Nm1));
		if constexpr (N % 2 == 1) {
			x_[N / 2] = zero<T>;
			w_[N / 2] = two<HiPrec> * inv(HiPrec(N * Nm1) * sqr(legendreP(Nm1, zero<T>)));
		}
	}
	constexpr T weight(int n) const {
		return w_[n];
	}
	constexpr T position(int n) const {
		return x_[n];
	}
	constexpr auto const &positions() const {
		return x_;
	}
	constexpr auto const &weights() const {
		return w_;
	}

private:
	Vector<T, size()> x_{};
	DiagonalMatrix<T, size()> w_{};
};

template <typename T, int N>
std::ostream &operator<<(std::ostream &os, GaussLegendreQuadrature<T, N> const &q) {
	os << "N = " << N << " Gauss-Legendre quadrature" << std::endl;
	os << "     position                 weight" << std::endl;
	for (int n = 0; n < N; n++) {
		os << print2string("%3i %+24.17e %24.17e\n", n, q.position(n), q.weight(n));
	}
	return os;
}

template <typename T, int N>
std::ostream &operator<<(std::ostream &os, GaussLobattoQuadrature<T, N> const &q) {
	os << "N = " << N << " Gauss-Lobatto quadrature" << std::endl;
	os << "     position                 weight" << std::endl;
	T wsum = zero<T>;
	for (int n = 0; n < N; n++) {
		os << print2string("%3i %+24.17e %24.17e\n", n, q.position(n), q.weight(n));
		wsum += q.weight(n);
	}
	os << "total weight = " << wsum << std::endl;
	return os;
}

template <typename T>
struct IsQuadrature : public std::false_type {};

template <typename T, int N>
struct IsQuadrature<GaussLegendreQuadrature<T, N>> : public std::true_type {};

template <typename T, int N>
struct IsQuadrature<GaussLobattoQuadrature<T, N>> : public std::true_type {};

template <typename T>
concept Quadrature = IsQuadrature<T>::value;

#endif /* QUADRATURE_HPP_ */
