/*
 * Matrix.hpp
 *
 *  Created on: Feb 4, 2026
 *      Author: dmarce1
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cassert>
#include <iostream>
#include <ostream>

#include "Math.hpp"
#include "Vector.hpp"

template <typename T, int N, int M = N>
struct Matrix;

template <typename T>
struct IsMatrix : std::false_type {};

template <typename T, int N, int M>
struct IsMatrix<Matrix<T, N, M>> : std::true_type {};

template <typename T, int N, int M>
struct Matrix {
	constexpr Matrix() = default;
	constexpr Matrix(T const &init) {
		std::fill(data_.begin(), data_.end(), init);
	}
	constexpr Matrix(std::initializer_list<std::initializer_list<T>> &&iLists) {
		*this = std::move(iLists);
	}
	constexpr Matrix &operator=(std::initializer_list<std::initializer_list<T>> &&iLists) {
		data_.fill(zero<T>);
		std::copy(iLists.begin(), iLists.end(), data_.begin());
		return *this;
	}
	constexpr Matrix const &operator+() const {
		return *this;
	}
	constexpr auto const &operator[](std::integral auto i) const {
		return data_[i];
	}
	constexpr auto &operator[](std::integral auto i) {
		return data_[i];
	}
	constexpr Matrix operator-() const {
		Matrix result;
		for (int i = 0; i < N; i++) {
			result[i] = -data_[i];
		}
		return result;
	}
	constexpr Matrix operator+(Matrix const &other) const {
		Matrix result;
		for (int i = 0; i < N; i++) {
			result[i] = data_[i] + other.data_[i];
		}
		return result;
	}
	constexpr Matrix operator-(Matrix const &other) const {
		Matrix result;
		for (int i = 0; i < N; i++) {
			result[i] = data_[i] - other.data_[i];
		}
		return result;
	}
	constexpr Matrix &operator+=(Matrix const &other) {
		*this = *this + other;
		return *this;
	}
	constexpr Matrix &operator-=(Matrix const &other) {
		*this = *this - other;
		return *this;
	}
	template <std::enable_if_t<std::is_same_v<decltype(T{1} * T{1}), T>, int> = 0>
	Matrix &operator*=(Matrix const &other) {
		static_assert(N == M);
		*this = *this * other;
		return *this;
	}
	constexpr Matrix &operator*=(Scalar auto &scalar) {
		*this = *this * scalar;
		return *this;
	}
	constexpr Matrix &operator/=(Scalar auto &scalar) {
		*this = *this / scalar;
		return *this;
	}
	template <Scalar S>
	constexpr Matrix<decltype(S{} * T{}), N, M> operator*(S const &scale) const {
		using R = decltype(S{} * T{});
		Matrix<R, N, M> result;
		for (int i = 0; i < N; i++) {
			result[i] = data_[i] * scale;
		}
		return result;
	}
	constexpr auto operator/(Scalar auto &scalar) const {
		return (*this) * inv(scalar);
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
	template <typename U, int L>
	constexpr auto operator*(Matrix<U, M, L> const &other) const {
		using R = decltype(U{} * T{});
		Matrix<R, N, L> result;
		for (int n = 0; n < N; n++) {
			for (int l = 0; l < L; l++) {
				result[n][l] = zero<R>;
				for (int m = 0; m < M; m++) {
					result[n][l] += data_[n][m] * other[m][l];
				}
			}
		}
		return result;
	}
	friend constexpr Matrix<T, M, N> transpose(Matrix const &other) {
		Matrix<T, M, N> result;
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < N; n++) {
				result[m][n] = other[n][m];
			}
		}
		return result;
	}
	template <Scalar S>
	friend constexpr Matrix<decltype(S{} * T{}), N, M> operator*(S &scalar, Matrix const &matrix) {
		return matrix * scalar;
	}
	static constexpr auto size() {
		return N * M;
	}
	static constexpr auto rowCount() {
		return N;
	}
	static constexpr auto colCount() {
		return M;
	}
	static constexpr auto identity() {
		Matrix<T, N, M> I{zero<T>};
		constexpr int L = std::min(N, M);
		for (int n = 0; n < L; n++) {
			I.data_[n][n] = one<T>;
		}
		return I;
	}
	static constexpr bool isZero(T const &value) {
		if constexpr (std::is_floating_point<T>::value) {
			constexpr T ε = std::sqrt(std::numeric_limits<T>::epsilon());
			return ((-ε < value) && (value < +ε));
		} else {
			return value == zero<T>;
		}
	}
	template <int N2, int M2>
	constexpr auto getBlock(std::integral auto j, std::integral auto k) const {
		Matrix<T, N2, M2> B;
		for (int n = 0; n < N2; n++) {
			for (int m = 0; m < M2; m++) {
				B[n][m] = (*this)[n + j][m + k];
			}
		}
		return B;
	}
	template <int N2, int M2>
	constexpr auto &setBlock(std::integral auto j, std::integral auto k, Matrix<T, N2, M2> const &B) {
		for (int n = 0; n < N2; n++) {
			for (int m = 0; m < M2; m++) {
				(*this)[n + j][m + k] = B[n][m];
			}
		}
		return *this;
	}
	friend constexpr auto inv(Matrix A) {
		Matrix<T, M, N> iA;
		if constexpr (N == M) {
			iA = identity();
			for (int n = 0; n < N; n++) {
				if (isZero(A[n][n])) {
					for (int j = n; j < N; j++) {
						T const Aⱼₙ = A[j][n];
						if (!isZero(Aⱼₙ)) {
							std::swap(iA[j], iA[n]);
							std::swap(A[j], A[n]);
							break;
						}
						assert(j != N - 1);
					}
				}
				T const iAₙₙ = inv(A[n][n]);
				iA[n] *= iAₙₙ;
				A[n] *= iAₙₙ;
				for (int j = n + 1; j < N; j++) {
					T const Aⱼₙ = A[j][n];
					if (!isZero(Aⱼₙ)) {
						iA[j] -= Aⱼₙ * iA[n];
						A[j] -= Aⱼₙ * A[n];
					}
				}
			}
			for (int n = N - 1; n >= 0; n--) {
				for (int j = 0; j < n; j++) {
					T const Aⱼₙ = A[j][n];
					if (!isZero(Aⱼₙ)) {
						iA[j] -= Aⱼₙ * iA[n];
						A[j] -= Aⱼₙ * A[n];
					}
				}
			}
		} else {
			constexpr int L = N + M;
			constexpr auto Iₙₙ = Matrix<T, M>::identity();
			constexpr auto Oₘₘ = Matrix<T, M>(zero<T>);
			Matrix<T, L> K;
			auto const trA = transpose(A);
			K.setBlock(0, 0, Iₙₙ);
			K.setBlock(0, N, A);
			K.setBlock(N, 0, trA);
			K.setBlock(N, N, Oₘₘ);
			iA = inv(K).template getBlock<M, N>(N, 0);
		}
		return iA;
	}
	friend constexpr auto det(Matrix A) {
		static_assert(N == M);
		T det = one<T>;
		for (int n = 0; n < N; n++) {
			T Aₙₙ = A[n][n];
			if (isZero(Aₙₙ)) {
				for (int j = n; j < N; j++) {
					if (j == N) {
						return zero<T>;
					}
					T const Aⱼₙ = A[j][n];
					if (!isZero(Aⱼₙ)) {
						std::swap(A[j], A[n]);
						Aₙₙ = A[n][n];
						det = -det;
						break;
					}
				}
			}
			A[n] *= inv(Aₙₙ);
			det *= Aₙₙ;
			for (int j = n + 1; j < N; j++) {
				T const Aⱼₙ = A[j][n];
				if (!isZero(Aⱼₙ)) {
					A[j] -= Aⱼₙ * A[n];
				}
			}
		}
		return det;
	}
	friend constexpr auto luDecomposition(Matrix LU) {
		static_assert(N == M);
		for (int n = 0; n < N; n++) {
			assert(!isZero(LU[n][n]));
			auto const iLUnn = inv(LU[n][n]);
			for (int j = n + 1; j < N; j++) {
				LU[j][n] *= iLUnn;
				for (int k = n + 1; k < N; k++) {
					LU[j][k] -= LU[j][n] * LU[n][k];
				}
			}
		}
		return LU;
	}

private:
	std::array<Vector<T, M>, N> data_{};
};

template <typename T, typename U, int N1, int N2>
constexpr auto operator*(Matrix<T, N1, N2> const &M, Vector<U, N2> const &V) {
	using R = decltype(T{} * U{});
	Vector<R, N1> Y;
	for (int n = 0; n < N1; n++) {
		Y[n] = M[n].dot(V);
	}
	return Y;
}

template <typename T, typename U, int N1, int N2>
constexpr auto operator*(Vector<U, N2> const &V, Matrix<T, N2, N1> const &M) {
	return transpose(M) * V;
}

template <typename T, int N>
struct DiagonalMatrix : public Vector<T, N> {
	operator Matrix<T, N>() const {
		Matrix<T, N> D{};
		for(int n = 0; n < N; n++) {
			D[n][n] = (*this)[n];
		}
		return D;
	}
	auto &operator*=(DiagonalMatrix<T, N> const &other) {
		*this = *this * other;
		return *this;
	}
	template <typename U>
	auto operator*(DiagonalMatrix<U, N> const &other) const {
		using R = decltype(U{} * T{});
		DiagonalMatrix<R, N> result;
		for (int n = 0; n < N; n++) {
			result[n] = (this)[n] * other[n];
		}
		return result;
	}
	template <typename U, int M>
	friend auto operator*(Matrix<U, M, N> const &A, DiagonalMatrix const &D) {
		using R = decltype(U{} * T{});
		Matrix<R, M, N> AD;
		for (int n = 0; n < N; n++) {
			T const d = D[n];
			for (int m = 0; m < M; m++) {
				AD[m][n] = A[m][n] * d;
			}
		}
		return AD;
	}
	template <typename U, int M>
	friend auto operator*(DiagonalMatrix const &D, Matrix<U, N, M> const &A) {
		using R = decltype(T{} * U{});
		Matrix<R, N, M> DA;
		for (int n = 0; n < N; n++) {
			T const d = D[n];
			for (int m = 0; m < M; m++) {
				DA[n][m] = d * A[n][m] ;
			}
		}
		return DA;
	}
	friend auto inv(DiagonalMatrix const &D) {
		using R = decltype(inv(T{1}));
		DiagonalMatrix<R, N> iD;
		for (int n = 0; n < N; n++) {
			iD[n] = inv(D[n]);
		}
		return iD;
	}
};

template <typename T, int N, int M>
std::ostream &operator<<(std::ostream &os, Matrix<T, N, M> const &A) {
	using MatrixType = Matrix<T, N, M>;
	std::array<std::array<std::string, M>, N> strs;
	size_t maxLen = 0;
	for (int n = 0; n < N; n++) {
		for (int m = 0; m < M; m++) {
			auto &str = strs[n][m];
			auto const &Anm = A[n][m];
			if (!MatrixType::isZero(Anm)) {
				str = std::to_string(Anm);
				maxLen = std::max(maxLen, str.size());
			}
		}
	}
	size_t const fieldLen = maxLen + 2;
	auto const hlineOut = [fieldLen, &os](char const *left, char const *mid, char const *right) {
		os << left;
		for (size_t j = 0; j < M; ++j) {
			for (size_t k = 0; k < fieldLen; k++) {
				os << "─";
			}
			if (j + 1 < M) {
				os << mid;
			}
		}
		os << right << std::endl;
	};
	hlineOut("┌", "┬", "┐");
	for (int n = 0; n < N; n++) {
		os << "|";
		for (int m = 0; m < M; m++) {
			auto const &str = strs[n][m];
			auto const leftPad = (maxLen - str.size() + 1) / 2 + 1;
			auto const rightPad = (maxLen - str.size()) / 2 + 1;
			os << std::string(leftPad, ' ') + strs[n][m] + std::string(rightPad, ' ') << '|';
		}
		os << std::endl;
		if (n + 1 < N) {
			hlineOut("├", "┼", "┤");
		}
	}
	hlineOut("└", "┴", "┘");
	return os;
}

#endif /* MATRIX_HPP_ */
