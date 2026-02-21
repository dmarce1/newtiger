/*
 * Vector.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include <array>
#include <cmath>

#include "Real.hpp"

template <Integer D>
class Vector {
	std::array<Real, D> data_{};

public:
	constexpr Vector() = default;
	constexpr Vector(Real ival) {
		data_.fill(ival);
	}
	constexpr Vector(std::initializer_list<Real> &&iList) {
		*this = std::move(iList);
	}
	constexpr Vector(std::array<Real, D> const &arr) {
		data_ = arr;
	}
	constexpr Vector &operator=(std::initializer_list<Real> const &iList) {
		data_.fill(0_R);
		std::copy(iList.begin(), iList.end(), data_.begin());
		return *this;
	}
	constexpr Vector &operator=(Real const &init) {
		data_.fill(init);
		return *this;
	}
	constexpr operator std::array<Real, D>() const {
		return data_;
	}
	constexpr Real const &operator[](Integer i) const {
		return data_[i];
	}
	constexpr Real &operator[](Integer i) {
		return data_[i];
	}
	constexpr Vector &operator+=(Vector const &other) {
		*this = *this + other;
		return *this;
	}
	constexpr Vector &operator-=(Vector const &other) {
		*this = *this - other;
		return *this;
	}
	constexpr Vector &operator*=(Real const &scalar) {
		*this = *this * scalar;
		return *this;
	}
	constexpr Vector &operator/=(Real const &scalar) {
		*this = *this / scalar;
		return *this;
	}
	constexpr Vector const &operator+() const {
		return *this;
	}
	constexpr Vector operator-() const {
		Vector result;
		for (Integer i = 0; i < D; i++) {
			result[i] = -data_[i];
		}
		return result;
	}
	constexpr Vector operator+(Vector const &other) const {
		Vector result;
		for (Integer i = 0; i < D; i++) {
			result[i] = data_[i] + other.data_[i];
		}
		return result;
	}
	constexpr Vector operator-(Vector const &other) const {
		Vector result;
		for (Integer i = 0; i < D; i++) {
			result[i] = data_[i] - other.data_[i];
		}
		return result;
	}
	constexpr auto operator*(Real const &scale) const {
		Vector result;
		for (Integer i = 0; i < D; i++) {
			result[i] = scale * data_[i];
		}
		return result;
	}
	constexpr auto operator/(Real const &scalar) const {
		return (*this) * inv(scalar);
	}
	constexpr auto dot(Vector const &other) const {
		Real result = 0_R;
		for (Integer i = 0; i < D; i++) {
			result += data_[i] * other.data_[i];
		}
		return result;
	}
	constexpr auto data() const {
		return data_.data();
	}
	constexpr auto data() {
		return data_.data();
	}
	constexpr auto begin() const {
		return data_.cbegin();
	}
	constexpr auto end() const {
		return data_.cend();
	}
	constexpr auto begin() {
		return data_.begin();
	}
	constexpr auto end() {
		return data_.end();
	}
	static constexpr auto size() {
		return D;
	}
	static constexpr auto unit(Integer i) {
		Vector u{};
		u[i] = 1_R;
		return u;
	}
	friend constexpr auto operator*(Real const &scalar, Vector const &vector) {
		return vector * scalar;
	}
	friend constexpr auto abs(Vector const &vector) {
		using std::sqrt;
		return sqrt(vector.dot(vector));
	}
	friend constexpr auto normalize(Vector const &vector) {
		return vector / abs(vector);
	}
};


template <Integer D>
inline constexpr auto minmod(Vector<D> const &a, Vector<D> const &b, Real θ) {
	Vector<D> c;
	for (Integer d = 0; d < D; d++) {
		c[d] = minmod(θ * minmod(a[d], b[d]), 0.5_R * (a[d] + b[d]));
	}
	return c;
}


#endif /* VECTOR_HPP_ */
