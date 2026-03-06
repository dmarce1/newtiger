/*
 * Quadrature.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef QUADRA2TURE_HPP_
#define QUADRA2TURE_HPP_

#include "Integer.hpp"
#include "Math.hpp"
#include "Real.hpp"
#include "Vector.hpp"

#include <array>
#include <cmath>
#include <numbers>

constexpr HiPrec legendreP(Integer l, HiPrec x) {
	if (l == 0_I) return 1_Hi;
	auto p1 = 1_Hi;
	auto p0 = 0_Hi;
	for (Integer n = 0; n < l; n++) {
		auto const p2 = (HiPrec(2_I * n + 1_I) * x * p1 - HiPrec(n) * p0) / HiPrec(n + 1_I);
		p0 = p1;
		p1 = p2;
	}
	return p1;
};

constexpr HiPrec dLegendrePdX(Integer l, HiPrec x) {
	if (l == 0_I) return 0_Hi;
	auto p1 = 1_Hi;
	auto p0 = 0_Hi, d0 = 0_Hi, d1 = 0_Hi;
	for (Integer n = 0; n < l; n++) {
		auto const p2 = (HiPrec(2_I * n + 1_I) * x * p1 - HiPrec(n) * p0) / HiPrec(n + 1_I);
		auto const d2 = (HiPrec(2_I * n + 1_I) * (p1 + x * d1) - HiPrec(n) * d0) / HiPrec(n + 1_I);
		p0 = p1;
		p1 = p2;
		d0 = d1;
		d1 = d2;
	}
	return d1;
};

template <Integer N>
constexpr auto gaussLegendreQuadrature() {
	constexpr HiPrec π = std::numbers::pi_v<HiPrec>;
	Vector<Real, N> x_, w_;
	for (Integer k = 0; 2 * k < N; k++) {
		HiPrec x, p0;
		HiPrec xm = cos(π * HiPrec(2_I * k + 2_I) / HiPrec(2_I * N + 1_I));
		HiPrec xp = cos(π * HiPrec(2_I * k + 1_I) / HiPrec(2_I * N + 1_I));
		HiPrec p1 = legendreP(N, xp);
		do {
			x = 0.5_Hi * (xm + xp);
			p0 = legendreP(N, x);
			if ((p1 * p0) > 0_Hi) {
				xp = x;
				p1 = p0;
			} else {
				xm = x;
			}
		} while (Real(xp) > Real(xm));
		x_[k] = -std::abs(x);
		x_[N - k - 1] = -x_[k];
		w_[k] = w_[N - k - 1] = 2_Hi * inv((1_Hi - sqr(x)) * sqr(dLegendrePdX(N, x)));
	}
	if constexpr (N % 2 == 1) {
		x_[N / 2] = 0_Hi;
		w_[N / 2] = 2_Hi * inv(sqr(dLegendrePdX(N, 0_Hi)));
	}
	return std::pair(x_, w_);
}

template <Integer N>
constexpr auto gaussLobattoQuadrature() {
	constexpr HiPrec π = std::numbers::pi_v<HiPrec>;
	auto const [glx, glw] = gaussLegendreQuadrature<N - 1>();
	Vector<Real, N> x_, w_;
	for (Integer k = 1; 2 * k < N; k++) {
		HiPrec x, p0;
		HiPrec xm = glx[k - 1];
		HiPrec xp = glx[k];
		HiPrec p1 = dLegendrePdX(N - 1, xp);
		do {
			x = 0.5_Hi * (xm + xp);
			p0 = dLegendrePdX(N - 1_I, x);
			if ((p1 * p0) > 0_Hi) {
				xp = x;
				p1 = p0;
			} else {
				xm = x;
			}
		} while (Real(xp - xm) > eps_R);
		x_[k] = -std::abs(x);
		x_[N - k - 1] = -x;
		w_[k] = w_[N - k - 1] = 2_Hi * inv(HiPrec(N * (N - 1_I)) * sqr(legendreP(N - 1_I, x)));
	}
	x_[0] = -1_Hi;
	x_[N - 1] = +1_Hi;
	w_[0] = w_[N - 1] = 2_Hi / HiPrec(N * (N - 1_I));
	if constexpr ((N & 1_I) == 1_I) {
		x_[N / 2] = 0_Hi;
		w_[N / 2] = 2_Hi * inv(HiPrec(N * (N - 1_I)) * sqr(legendreP(N - 1_I, 0_Hi)));
	}
	return std::pair(x_, w_);
}

#endif /* QUADRATURE_HPP_ */
