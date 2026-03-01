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

#include "Math.hpp"
#include "Real.hpp"
#include "Select.hpp"

#define ZERO T(0_I)
#define ONE T(1_I)

template <typename T, Integer D>
class Vector {
	std::array<T, D> data_{};

public:
	constexpr Vector() = default;
	constexpr Vector(T ival) {
		data_.fill(ival);
	}
	constexpr Vector(std::initializer_list<T> &&iList) {
		*this = std::move(iList);
	}
	constexpr Vector(std::array<T, D> const &arr) {
		data_ = arr;
	}
	constexpr Vector &operator=(std::initializer_list<T> const &iList) {
		data_.fill(ZERO);
		std::copy(iList.begin(), iList.end(), data_.begin());
		return *this;
	}
	constexpr Vector &operator=(T const &init) {
		data_.fill(init);
		return *this;
	}
	constexpr operator std::array<T, D>() const {
		return data_;
	}
	constexpr T const &operator[](Integer i) const {
		return data_[i];
	}
	constexpr T &operator[](Integer i) {
		return data_[i];
	}
	template <typename U>
	constexpr Vector &operator+=(Vector<U, D> const &other) {
		*this = *this + other;
		return *this;
	}
	template <typename U>
	constexpr Vector &operator-=(Vector<U, D> const &other) {
		*this = *this - other;
		return *this;
	}
	constexpr Vector &operator*=(T const &scalar) {
		*this = *this * scalar;
		return *this;
	}
	constexpr Vector &operator/=(T const &scalar) {
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
	template <typename U>
	constexpr auto operator+(Vector<U, D> const &other) const {
		using R = decltype(T{} + U{});
		Vector<R, D> result;
		for (Integer i = 0; i < D; i++) {
			result[i] = data_[i] + other[i];
		}
		return result;
	}
	template <typename U>
	constexpr auto operator-(Vector<U, D> const &other) const {
		using R = decltype(T{} - U{});
		Vector<R, D> result;
		for (Integer i = 0; i < D; i++) {
			result[i] = data_[i] - other[i];
		}
		return result;
	}
	constexpr auto operator*(T const &scale) const {
		Vector result;
		for (Integer i = 0; i < D; i++) {
			result[i] = scale * data_[i];
		}
		return result;
	}
	constexpr auto operator/(T const &scalar) const {
		return (*this) * inv(scalar);
	}
	template <typename U>
	constexpr auto dot(Vector<U, D> const &other) const {
		auto result = data_[0] * other.data_[0];
		for (Integer i = 1; i < D; i++) {
			result += data_[i] * other[i];
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
	constexpr auto &front() {
		return data_[0];
	}
	constexpr auto &back() {
		return data_[D - 1];
	}
	constexpr auto const &front() const {
		return data_[0];
	}
	constexpr auto const &back() const {
		return data_[D - 1];
	}
	static constexpr auto size() {
		return D;
	}
	static constexpr auto unit(Integer i) {
		Vector u{};
		u[i] = 1_R;
		return u;
	}
	friend constexpr auto operator*(T const &scalar, Vector const &vector) {
		return vector * scalar;
	}
	friend constexpr auto abs(Vector const &vector) {
		using std::sqrt;
		return sqrt(vector.dot(vector));
	}
	constexpr auto normalize() const {
		auto const v0 = abs(*this);
		auto const den = select(v0 != ZERO, v0, ONE);
		return *this / den;
	}
	template <Integer B, Integer E>
	constexpr auto sub() const {
		Vector<T, E - B> result;
		std::copy(begin() + B, begin() + E, result.begin());
		return result;
	}
	template <Integer D2>
	friend constexpr auto concatenate(Vector const &v1, Vector<T, D2> const &v2) {
		Vector<T, D + D2> result;
		std::copy(v1.begin(), v1.end(), result.begin());
		std::copy(v2.begin(), v2.end(), result.begin() + D);
		return result;
	}
	friend constexpr auto concatenate(T const &s1, Vector<T, D> const &v2) {
		Vector<T, D + 1> result;
		result[0] = s1;
		std::copy(v2.begin(), v2.end(), result.begin() + 1);
		return result;
	}
	friend constexpr auto concatenate(Vector<T, D> const &v1, T const &s2) {
		Vector<T, D + 1> result;
		result[D] = s2;
		std::copy(v1.begin(), v1.end(), result.begin());
		return result;
	}
};

template <typename T, Integer W>
std::ostream &operator<<(std::ostream &os, Vector<T, W> const &v) {
	os << "(";
	os << v[0];
	for (Integer i = 1; i < W; i++) {
		os << std::string(", ") << v[i];
	}
	os << ")";
	return os;
}

#endif /* VECTOR_HPP_ */
