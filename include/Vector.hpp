/*
 * Vector.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "Constants.hpp"
#include "Math.hpp"

#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <numeric>
#include <ostream>

template <typename, int>
struct Vector;

template <typename T>
struct IsVector : std::false_type {};

template <typename T, int N>
struct IsVector<Vector<T, N>> : std::true_type {};

template <typename T, int L>
struct Vector {
	constexpr static int N = L;
	constexpr Vector() = default;
	constexpr Vector(T ival) {
		data_.fill(ival);
	}
	constexpr Vector(std::initializer_list<T> &&iList) {
		*this = std::move(iList);
	}
	constexpr Vector(std::array<T, L> const &arr) {
		data_ = arr;
	}
	constexpr Vector &operator=(std::initializer_list<T> const &iList) {
		data_.fill(zero<T>);
		std::copy(iList.begin(), iList.end(), data_.begin());
		return *this;
	}
	constexpr Vector &operator=(T const &init) {
		data_.fill(init);
		return *this;
	}
	constexpr operator std::array<T, N>() const {
		return data_;
	}
	constexpr T const &operator[](std::integral auto i) const {
		return data_[i];
	}
	constexpr T &operator[](std::integral auto i) {
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
	constexpr Vector &operator*=(Scalar auto const &scalar) {
		*this = *this * scalar;
		return *this;
	}
	constexpr Vector &operator/=(Scalar auto const &scalar) {
		*this = *this / scalar;
		return *this;
	}
	constexpr Vector const &operator+() const {
		return *this;
	}
	constexpr Vector operator-() const {
		Vector result;
		for (int i = 0; i < N; i++) {
			result[i] = -data_[i];
		}
		return result;
	}
	constexpr Vector operator+(Vector const &other) const {
		Vector result;
		for (int i = 0; i < N; i++) {
			result[i] = data_[i] + other.data_[i];
		}
		return result;
	}
	constexpr Vector operator-(Vector const &other) const {
		Vector result;
		for (int i = 0; i < N; i++) {
			result[i] = data_[i] - other.data_[i];
		}
		return result;
	}
	constexpr auto operator*(Scalar auto const &scale) const {
		using R = decltype(scale * T{});
		Vector<R, N> result;
		for (int i = 0; i < N; i++) {
			result[i] = scale * data_[i];
		}
		return result;
	}
	constexpr auto operator/(Scalar auto const &scalar) const {
		return (*this) * inv(scalar);
	}
	template <typename U>
	constexpr auto dot(Vector<U, N> const &other) const {
		auto result = data_[0] * other.data_[0];
		for (int i = 1; i < N; i++) {
			result += data_[i] * other.data_[i];
		}
		return result;
	}
	constexpr auto data()  const{
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
		return L;
	}
	static constexpr auto unit(int i) {
		Vector u{};
		u[i] = one<T>;
		return u;
	}
	friend constexpr auto operator*(Scalar auto const &scalar, Vector const &vector) {
		return vector * scalar;
	}
	friend constexpr auto abs(Vector const &vector) {
		using std::sqrt;
		return sqrt(vector.dot(vector));
	}
	friend constexpr auto normalize(Vector const &vector) {
		return vector / abs(vector);
	}
	template<std::integral Int>
	T &operator()(std::array<Int, 1> const &xyz) {
		return data_[xyz[0]];
	}
	template<std::integral Int>
	T operator()(std::array<Int, 1> const &xyz) const {
		return data_[xyz[0]];
	}

private:
	std::array<T, N> data_{};
};

#endif /* VECTOR_HPP_ */
