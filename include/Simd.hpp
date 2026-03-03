/*
 * Simd.hpp
 *
 *  Created on: Feb 23, 2026
 *      Author: dmarce1
 */

#ifndef SIMD_HPP_
#define SIMD_HPP_

#include "Concepts.hpp"
#include "Integer.hpp"
#include "Real.hpp"

#include <climits>
#include <functional>
#include <ostream>

inline constexpr Integer maxSimdBitCount = 256;

template <typename T>
inline constexpr Integer maxSimdSize() {
	return maxSimdBitCount / (sizeof(T) * CHAR_BIT);
}

template <typename T>
concept SimdScalar = std::is_same_v<T, Integer> || std::is_same_v<T, Real> || std::is_same_v<T, bool>;

#define SIMD_ASSIGNMENT_OP(op)                                                                                                             \
	template <typename OtherType>                                                                                                          \
	inline Simd &operator op(Simd<OtherType, W> const &other) noexcept {                                                                   \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			vec_[i] op other[i];                                                                                                           \
		}                                                                                                                                  \
		return *this;                                                                                                                      \
	}                                                                                                                                      \
	inline Simd &operator op(SimdScalar auto const &a) noexcept {                                                                          \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			vec_[i] op a;                                                                                                                  \
		}                                                                                                                                  \
		return *this;                                                                                                                      \
	}

#define SIMD_UNARY_OP(op)                                                                                                                  \
	friend inline Simd operator op(Simd const &A) noexcept {                                                                               \
		Simd C;                                                                                                                            \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C.vec_[i] = op A.vec_[i];                                                                                                      \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}

#define SIMD_BINARY_OP(op)                                                                                                                 \
	friend inline Simd operator op(Simd const &A, Simd const &B) noexcept {                                                                \
		Simd C;                                                                                                                            \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C.vec_[i] = A.vec_[i] op B.vec_[i];                                                                                            \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}                                                                                                                                      \
	friend inline Simd operator op(T const &a, Simd const &B) noexcept {                                                                   \
		Simd C;                                                                                                                            \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C.vec_[i] = a op B.vec_[i];                                                                                                    \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}                                                                                                                                      \
	friend inline Simd operator op(Simd const &A, T const &b) noexcept {                                                                   \
		Simd C;                                                                                                                            \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C.vec_[i] = A.vec_[i] op b;                                                                                                    \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}

#define SIMD_COMPARE_OP(op)                                                                                                                \
	friend inline Simd<bool, W> operator op(Simd const &A, Simd const &B) noexcept {                                                       \
		Simd<bool, W> C;                                                                                                                   \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C[i] = A.vec_[i] op B.vec_[i];                                                                                                 \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}                                                                                                                                      \
	friend inline Simd<bool, W> operator op(T const &a, Simd const &B) noexcept {                                                          \
		Simd<bool, W> C;                                                                                                                   \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C[i] = a op B.vec_[i];                                                                                                         \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}                                                                                                                                      \
	friend inline Simd<bool, W> operator op(Simd const &A, T const &b) noexcept {                                                          \
		Simd<bool, W> C;                                                                                                                   \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C[i] = A.vec_[i] op b;                                                                                                         \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}

#define SIMD_UNARY_FUNCTION(func)                                                                                                          \
	friend inline Simd func(Simd const &A) noexcept {                                                                                      \
		using namespace std;                                                                                                               \
		Simd C;                                                                                                                            \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C.vec_[i] = func(A.vec_[i]);                                                                                                   \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}

#define SIMD_BINARY_FUNCTION(func)                                                                                                         \
	friend inline Simd func(Simd const &A, Simd const &B) noexcept {                                                                       \
		using namespace std;                                                                                                               \
		Simd C;                                                                                                                            \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C.vec_[i] = func(A.vec_[i], B.vec_[i]);                                                                                        \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}                                                                                                                                      \
	friend inline Simd func(Simd const &A, T const &b) noexcept {                                                                          \
		using namespace std;                                                                                                               \
		Simd C;                                                                                                                            \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C.vec_[i] = func(A.vec_[i], b);                                                                                                \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}                                                                                                                                      \
	friend inline Simd func(T const &a, Simd const &B) noexcept {                                                                          \
		using namespace std;                                                                                                               \
		Simd C;                                                                                                                            \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			C.vec_[i] = func(a, B.vec_[i]);                                                                                                \
		}                                                                                                                                  \
		return C;                                                                                                                          \
	}

#define SIMD_REDUCTION1(func)                                                                                                              \
	friend inline T func(Simd const &simd) noexcept {                                                                                      \
		using namespace std;                                                                                                               \
		T result = simd.vec_[0];                                                                                                           \
		for (Integer i = 1; i < W; i++) {                                                                                                  \
			result = func(result, simd.vec_[i]);                                                                                           \
		}                                                                                                                                  \
		return result;                                                                                                                     \
	}

#define SIMD_REDUCTION2(name, op)                                                                                                          \
	friend inline T name(Simd const &simd) noexcept {                                                                                      \
		T result = simd.vec_[0];                                                                                                           \
		for (Integer i = 1; i < W; i++) {                                                                                                  \
			result = result op simd.vec_[i];                                                                                               \
		}                                                                                                                                  \
		return result;                                                                                                                     \
	}

#define SIMD_METHODS()                                                                                                                     \
	static_assert(W <= maxSimdBitCount / (CHAR_BIT * sizeof(Real)));                                                                       \
	static consteval Integer size() noexcept {                                                                                             \
		return W;                                                                                                                          \
	}                                                                                                                                      \
	inline constexpr Simd() = default;                                                                                                     \
	template <typename OtherType>                                                                                                          \
	inline Simd(Simd<OtherType, W> const &other) noexcept {                                                                                \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			vec_[i] = other[i];                                                                                                            \
		}                                                                                                                                  \
	}                                                                                                                                      \
	inline Simd(SimdScalar auto const &a) noexcept {                                                                                       \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			vec_[i] = a;                                                                                                                   \
		}                                                                                                                                  \
	}                                                                                                                                      \
	inline T &operator[](Integer i) noexcept {                                                                                             \
		return vec_[i];                                                                                                                    \
	}                                                                                                                                      \
	inline T operator[](Integer i) const noexcept {                                                                                        \
		return vec_[i];                                                                                                                    \
	}                                                                                                                                      \
	inline T *begin() noexcept {                                                                                                           \
		return vec_;                                                                                                                       \
	}                                                                                                                                      \
	inline T *end() noexcept {                                                                                                             \
		return vec_ + W;                                                                                                                   \
	}                                                                                                                                      \
	inline T const *begin() const noexcept {                                                                                               \
		return vec_;                                                                                                                       \
	}                                                                                                                                      \
	inline T const *end() const noexcept {                                                                                                 \
		return vec_ + W;                                                                                                                   \
	}                                                                                                                                      \
	inline void load(T const *ptr) {                                                                                                       \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			vec_[i] = ptr[i];                                                                                                              \
		}                                                                                                                                  \
	}                                                                                                                                      \
	inline void store(T *ptr) const {                                                                                                      \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			ptr[i] = vec_[i];                                                                                                              \
		}                                                                                                                                  \
	}                                                                                                                                      \
	inline void serialize(auto &arc, unsigned) {                                                                                           \
		for (Integer i = 0; i < W; i++) {                                                                                                  \
			arc &vec_[i];                                                                                                                  \
		}                                                                                                                                  \
	}

#define SIMD_COMPATIBILITY_REDUCTION(Type, func)                                                                                           \
	inline constexpr Type func(Type a) noexcept {                                                                                                    \
		return a;                                                                                                                          \
	}                                                                                                                                      \
	template <typename T, Integer W>                                                                                                       \
	inline constexpr Type func(Simd<T, W> a) noexcept {                                                                                              \
		return a.func();                                                                                                                   \
	}

template <typename T, Integer W>
struct Simd;

template <Integer W>
struct Simd<bool, W> {
	using T = bool;
	SIMD_METHODS();
	SIMD_ASSIGNMENT_OP(=);
	SIMD_UNARY_OP(!);
	SIMD_BINARY_OP(&&);
	SIMD_BINARY_OP(||);
	SIMD_COMPARE_OP(==);
	SIMD_COMPARE_OP(!=);
	SIMD_COMPARE_OP(>=);
	SIMD_COMPARE_OP(<=);
	SIMD_COMPARE_OP(<);
	SIMD_COMPARE_OP(>);
	SIMD_REDUCTION2(any, ||);
	SIMD_REDUCTION2(all, &&);

private:
	bool vec_[W];
};

template <Integer W>
struct Simd<Integer, W> {
	using T = Integer;
	SIMD_METHODS();
	SIMD_ASSIGNMENT_OP(=);
	SIMD_ASSIGNMENT_OP(+=);
	SIMD_ASSIGNMENT_OP(-=);
	SIMD_ASSIGNMENT_OP(*=);
	SIMD_ASSIGNMENT_OP(/=);
	SIMD_ASSIGNMENT_OP(%=);
	SIMD_ASSIGNMENT_OP(|=);
	SIMD_ASSIGNMENT_OP(&=);
	SIMD_ASSIGNMENT_OP(^=);
	SIMD_ASSIGNMENT_OP(>>=);
	SIMD_ASSIGNMENT_OP(<<=);
	SIMD_UNARY_OP(+);
	SIMD_UNARY_OP(-);
	SIMD_UNARY_OP(~);
	SIMD_BINARY_OP(+);
	SIMD_BINARY_OP(-);
	SIMD_BINARY_OP(*);
	SIMD_BINARY_OP(/);
	SIMD_BINARY_OP(%);
	SIMD_BINARY_OP(|);
	SIMD_BINARY_OP(&);
	SIMD_BINARY_OP(^);
	SIMD_BINARY_OP(>>);
	SIMD_BINARY_OP(<<);
	SIMD_UNARY_FUNCTION(abs);
	SIMD_BINARY_FUNCTION(copysign);
	SIMD_BINARY_FUNCTION(max);
	SIMD_BINARY_FUNCTION(min);
	SIMD_COMPARE_OP(==);
	SIMD_COMPARE_OP(!=);
	SIMD_COMPARE_OP(>=);
	SIMD_COMPARE_OP(<=);
	SIMD_COMPARE_OP(<);
	SIMD_COMPARE_OP(>);
	SIMD_REDUCTION1(max);
	SIMD_REDUCTION1(min);
	SIMD_REDUCTION2(sum, +);
	SIMD_REDUCTION2(product, *);

private:
	T vec_[W];
};

template <Integer W>
struct Simd<Real, W> {
	using T = Real;
	SIMD_METHODS();
	SIMD_ASSIGNMENT_OP(=);
	SIMD_ASSIGNMENT_OP(+=);
	SIMD_ASSIGNMENT_OP(-=);
	SIMD_ASSIGNMENT_OP(*=);
	SIMD_ASSIGNMENT_OP(/=);
	SIMD_UNARY_OP(+);
	SIMD_UNARY_OP(-);
	SIMD_BINARY_OP(+);
	SIMD_BINARY_OP(-);
	SIMD_BINARY_OP(*);
	SIMD_BINARY_OP(/);
	SIMD_UNARY_FUNCTION(abs);
	SIMD_UNARY_FUNCTION(sqrt);
	SIMD_BINARY_FUNCTION(copysign);
	SIMD_BINARY_FUNCTION(max);
	SIMD_BINARY_FUNCTION(min);
	SIMD_BINARY_FUNCTION(pow);
	SIMD_COMPARE_OP(==);
	SIMD_COMPARE_OP(!=);
	SIMD_COMPARE_OP(>=);
	SIMD_COMPARE_OP(<=);
	SIMD_COMPARE_OP(<);
	SIMD_COMPARE_OP(>);
	SIMD_REDUCTION1(max);
	SIMD_REDUCTION1(min);
	SIMD_REDUCTION2(sum, +);
	SIMD_REDUCTION2(product, *);

private:
	T vec_[W];
};

SIMD_COMPATIBILITY_REDUCTION(bool, any);
SIMD_COMPATIBILITY_REDUCTION(bool, all);
SIMD_COMPATIBILITY_REDUCTION(Integer, max);
SIMD_COMPATIBILITY_REDUCTION(Integer, min);
SIMD_COMPATIBILITY_REDUCTION(Integer, sum);
SIMD_COMPATIBILITY_REDUCTION(Integer, product);
SIMD_COMPATIBILITY_REDUCTION(Real, max);
SIMD_COMPATIBILITY_REDUCTION(Real, min);
SIMD_COMPATIBILITY_REDUCTION(Real, sum);
SIMD_COMPATIBILITY_REDUCTION(Real, product);

template <typename T, Integer W>
std::ostream &operator<<(std::ostream &os, Simd<T, W> const &simd) {
	os << "(";
	os << std::to_string(simd[0]);
	for (Integer i = 1; i < W; i++) {
		os << std::string(", ") << simd[i];
	}
	os << ")";
	return os;
}

#endif /* SIMD_HPP_ */
