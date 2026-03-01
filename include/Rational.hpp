/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#pragma once

#include "Math.hpp"

#include <cmath>
#include <concepts>
#include <cstdint>
#include <numeric>
#include <ostream>

struct Rational {
	using Type = intmax_t;
	constexpr Rational() :
		n_(0), d_(1) {
	}
	constexpr Rational(std::integral auto n) :
		n_(n), d_(1) {
	}
	constexpr Rational(std::integral auto n, std::integral auto d) :
		n_(n), d_(d) {
	}
	constexpr Rational &operator=(std::integral auto i) {
		n_ = i;
		d_ = 1;

		return *this;
	}
	constexpr Rational &operator=(Rational i) {
		n_ = i.n_;
		d_ = i.d_;
		return *this;
	}
	constexpr Rational &operator+=(Rational i) {
		*this = *this + i;
		return *this;
	}
	constexpr Rational &operator-=(Rational i) {
		*this = *this - i;
		return *this;
	}
	constexpr Rational &operator*=(Rational i) {
		*this = *this * i;
		return *this;
	}
	constexpr Rational &operator/=(Rational i) {
		*this = *this / i;
		return *this;
	}
	friend constexpr Rational operator+(Rational b) {
		return b;
	}
	friend constexpr Rational operator-(Rational a) {
		a.n_ = -a.n_;
		return a;
	}
	friend constexpr Rational operator+(Rational b, Rational c) {
		if (b.n_ == 0) {
			return c;
		} else if (c.n_ == 0) {
			return b;
		} else {
			Rational a;
			auto g = std::gcd(b.d_, c.d_);
			auto const bd1 = b.d_ / g;
			auto const cd1 = c.d_ / g;
			a.n_ = b.n_ * cd1 + c.n_ * bd1;
			a.d_ = c.d_ * bd1;
			if (a.n_) {
				g = std::gcd(a.n_, a.d_);
				a.n_ /= g;
				a.d_ /= g;
			} else {
				a.d_ = 1;
			}
			return a;
		}
	}
	friend constexpr Rational operator-(Rational a, Rational b) {
		return a + (-b);
	}
	friend constexpr Rational operator*(Rational b, Rational c) {
		if ((b.n_ == 0) && (c.n_ == 0)) {
			return Rational(0);
		} else {
			Rational a;
			a.n_ = b.n_ * c.n_;
			a.d_ = b.d_ * c.d_;
			if (a.n_) {
				auto const g = std::gcd(a.n_, a.d_);
				a.n_ /= g;
				a.d_ /= g;
			}
			return a;
		}
	}
	friend constexpr Rational operator*(std::integral auto b, Rational c) {
		if ((b == 0) && (c.n_ == 0)) {
			return Rational(0);
		} else {
			Rational a;
			a.n_ = b * c.n_;
			a.d_ = c.d_;
			if (a.n_) {
				auto const g = std::gcd(a.n_, a.d_);
				a.n_ /= g;
				a.d_ /= g;
			}
			return a;
		}
	}
	friend constexpr Rational operator*(Rational b, std::integral auto c) {
		return c * b;
	}
	constexpr Rational reciprocal() const {
		Rational a = *this;
		std::swap(a.n_, a.d_);
		if (a.d_ < 0) {
			a.n_ = -a.n_;
			a.d_ = -a.d_;
		}
		return a;
	}
	friend constexpr Rational operator/(Rational b, Rational c) {
		return b * c.reciprocal();
	}
	friend constexpr Rational operator/(std::integral auto b, Rational c) {
		return b * c.reciprocal();
	}
	friend constexpr Rational operator/(Rational b, std::integral auto c) {
		return (c / b).reciprocal();
	}
	friend constexpr bool operator==(Rational a, Rational b) {
		auto const gn = std::gcd(a.n_, b.n_);
		auto const gd = std::gcd(a.d_, b.d_);
		if (gn) {
			a.n_ /= gn;
			b.n_ /= gn;
		}
		a.d_ /= gd;
		b.d_ /= gd;
		return (a.n_ * b.d_) == (a.d_ * b.n_);
	}
	friend constexpr bool operator<(Rational a, Rational b) {
		auto const gn = std::gcd(a.n_, b.n_);
		auto const gd = std::gcd(a.d_, b.d_);
		if (gn) {
			a.n_ /= gn;
			b.n_ /= gn;
		}
		a.d_ /= gd;
		b.d_ /= gd;
		return (a.n_ * b.d_) < (a.d_ * b.n_);
	}
	friend constexpr bool operator!=(Rational a, Rational b) {
		return !(a == b);
	}
	friend constexpr bool operator>=(Rational a, Rational b) {
		return !(a < b);
	}
	friend constexpr bool operator>(Rational a, Rational b) {
		return b < a;
	}
	friend constexpr bool operator<=(Rational a, Rational b) {
		return b >= a;
	}
	constexpr operator double() const {
		return double(n_) / double(d_);
	}
	constexpr Type denominator() const {
		return d_;
	}
	constexpr Type numerator() const {
		return n_;
	}
	friend constexpr Rational inv(Rational v) {
		std::swap(v.n_, v.d_);
		v.n_ = v.d_ > 0_I ? +v.n_ : -v.n_;
		v.d_ = v.d_ > 0_I ? +v.d_ : -v.d_;
		return v;
	}
	friend constexpr Rational abs(Rational v) {
		v.n_ = (v.n_ >= 0) ? +v.n_ : -v.n_;
		return v;
	}
	friend constexpr auto pow(auto x, Rational r) {
		return pow(root(x, r.d_), r.n_);
	}
	friend std::ostream &operator<<(std::ostream &os, Rational f) {
		os << "(" << f.n_;
		if (f.d_ != 1) {
			os << "/" << f.d_;
		}
		os << ")";
		return os;
	}

private:
	Type n_ = 0;
	Type d_ = 1;
};
