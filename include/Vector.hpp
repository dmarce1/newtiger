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

template <typename, int, int = 1>
struct Vector;

template <typename T>
struct IsVector : std::false_type {};

template <typename T, int N, int D>
struct IsVector<Vector<T, N, D>> : std::true_type {};

template <typename T, int L>
struct Vector<T, L, 1> {
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

template<typename T, int N, int D>
constexpr Vector<T, N, D> fma(Scalar auto const& a, Vector<T, N, D> const& b, Vector<T, N, D> c) {
	for(int i = 0; i < N; i++) {
		c[i] += a * b[i];
	}
	return c;
}


template <typename T, int L, int D>
struct Vector : public Vector<Vector<T, L, D - 1>, L> {
	using base_type = Vector<Vector<T, L, D - 1>, L>;
	constexpr Vector() = default;
	constexpr Vector(base_type const &other) :
		base_type(other) {
	}
	template<std::integral Int>
	T &operator()(std::array<Int, D> const &nml) {
		T *ptr = reinterpret_cast<T *>(this);
		return ptr[flatIndex(nml)];
	}
	template<std::integral Int>
	T operator()(std::array<Int, D> const &nml) const {
		T const *ptr = reinterpret_cast<T const *>(this);
		return ptr[flatIndex(nml)];
	}

private:
	static constexpr int flatIndex(std::array<int, D> const &xyz) {
		int index = 0;
		for (int d = 0; d < D; d++) {
			index = L * index + xyz[d];
		}
		return index;
	}
};

//template <typename T, int L, int D>
//struct TriangularVector : public Vector<T, binco(L + D - 1, D)> {
//	static constexpr int size() {
//		return binco(L + D - 1, D);
//	}
//	using base_type = Vector<Vector<T, size()>, L>;
//	constexpr TriangularVector() = default;
//	constexpr TriangularVector(base_type const &other) :
//		base_type(other) {
//	}
//	T &operator()(std::array<int, D> const &nml) {
//		T *ptr = reinterpret_cast<T *>(this);
//		return ptr[flatIndex(nml)];
//	}
//	T operator()(std::array<int, D> const &nml) const {
//		T const *ptr = reinterpret_cast<T const *>(this);
//		return ptr[flatIndex(nml)];
//	}
//	template<int I>
//	auto getSubvector() {
//		
//	}
//
//private:
//	static constexpr int flatIndex(std::array<int, D> nml) {
//		std::reverse(nml.begin(), nml.end());
//		for (int d = 0; d + 1 < D; d++) {
//			nml[d + 1] += nml[d];
//		}
//		int index = 0;
//		for (int d = 0; d < D; d++) {
//			index += binco(nml[d] + d, d + 1);
//		}
//		return index;
//	}
//};
//
template <typename TA, typename TB>
constexpr auto cross(Vector<TA, 3> const &A, Vector<TB, 3> const &B) {
	using TC = decltype(TA{} * TB{});
	Vector<TC, 3> C;
	C[0] = A[1] * B[2] - A[2] * B[1];
	C[1] = A[2] * B[0] - A[0] * B[2];
	C[2] = A[0] * B[1] - A[1] * B[0];
	return C;
}

template <typename T, int N>
std::ostream &operator<<(std::ostream &os, Vector<T, N> const &A) {
	os << "(";
	for (int n = 0; n < N; n++) {
		if (n) {
			os << ", ";
		}
		os << A[n];
	}
	os << ")";
	return os;
}

#endif /* VECTOR_HPP_ */
