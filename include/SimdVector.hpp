/*
 * SimdVector.hpp
 *
 *  Created on: Feb 10, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_SIMDVECTOR_HPP_
#define INCLUDE_SIMDVECTOR_HPP_

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#include <array>
#include <bitset>
#include <cmath>
#include <concepts>
#include <type_traits>

// template <typename C>
// concept SimdExpressionType = requires(C &c, C const &cc, std::size_t i) {
//	typename C::value_type;
//	{ cc.size() } -> std::integral;
//	requires(!std::same_as<decltype(c[i]), void>);
//	requires(!std::same_as<decltype(cc[i]), void>);
// };
//
template <typename, size_t, typename>
struct SimdExpression;

template <typename S, size_t W>
struct SimdStorage {
	using value_type = S;
	constexpr size_t size() const {
		return W;
	}
	S operator()(int i) const {
		return vec_[i];
	}
	S &operator()(int i) {
		return vec_[i];
	}

private:
	std::array<S, W> vec_;
};

#define SIMD_UNARY_OPERATOR(op)                                                                                                            \
	auto operator op() const {                                                                                                             \
		auto lambda = [this](int i) {                                                                                                      \
			return op f_(i);                                                                                                               \
		};                                                                                                                                 \
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));                                                                  \
	}

#define SIMD_BINARY_OPERATOR(op)                                                                                                           \
	template <typename F_>                                                                                                                 \
	auto operator op(SimdExpression<S, W, F_> const &other) const {                                                                        \
		auto lambda = [this, other](int i) {                                                                                               \
			return f_(i) op other(i);                                                                                                      \
		};                                                                                                                                 \
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));                                                                  \
	}                                                                                                                                      \
	auto operator op(S s) const {                                                                                                          \
		auto lambda = [this, s](int i) {                                                                                                   \
			return f_(i) op s;                                                                                                             \
		};                                                                                                                                 \
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));                                                                  \
	}                                                                                                                                      \
	friend auto operator op(S s, SimdExpression const &expr) {                                                                             \
		auto lambda = [expr, s](int i) {                                                                                                   \
			return s op expr(i);                                                                                                           \
		};                                                                                                                                 \
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));                                                                  \
	}

#define SIMD_UNARY_FUNCTION(func)                                                                                                          \
	friend auto func(SimdExpression const &expr) {                                                                                         \
		using std::func;                                                                                                                   \
		auto lambda = [expr](int i) {                                                                                                      \
			return func(expr(i));                                                                                                          \
		};                                                                                                                                 \
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));                                                                  \
	}

#define SIMD_BINARY_FUNCTION(func)                                                                                                         \
	template <typename F_>                                                                                                                 \
	friend auto func(SimdExpression const &expr1, SimdExpression<S, W, F_> const &expr2) {                                                 \
		auto lambda = [expr1, expr2](int i) {                                                                                              \
			return func(expr1(i), expr2(i));                                                                                               \
		};                                                                                                                                 \
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));                                                                  \
	}                                                                                                                                      \
	friend auto func(SimdExpression const &expr1, S s) {                                                                                   \
		auto lambda = [expr1, s](int i) {                                                                                                  \
			return func(expr1(i), s);                                                                                                      \
		};                                                                                                                                 \
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));                                                                  \
	}                                                                                                                                      \
	friend auto func(S s, SimdExpression const &expr2) {                                                                                   \
		auto lambda = [s, expr2](int i) {                                                                                                  \
			return func(s, expr2(i));                                                                                                      \
		};                                                                                                                                 \
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));                                                                  \
	}

#define SIMD_COMPARE_OPERATOR(op)                                                                                                          \
	template <typename F_>                                                                                                                 \
	auto operator op(SimdExpression<S, W, F_> const &other) const {                                                                        \
		auto lambda = [this, other](int i) {                                                                                               \
			return (*this)(i)op other(i);                                                                                                  \
		};                                                                                                                                 \
		return SimdExpression<bool, W, decltype(lambda)>(std::move(lambda));                                                               \
	}

template <size_t W>
using SimdMask = std::bitset<W>;

template <typename S, size_t W, typename F>
struct SimdExpression {
	using value_type = S;
	constexpr size_t size() const {
		return W;
	}
	SimdExpression(F &&f) :
		f_{std::move(f)} {
	}
	S operator()(int i) const {
		return f_(i);
	}
	S &operator()(int i) {
		return f_(i);
	}
	auto cshift(int di) const {
		auto lambda = [di, this](int i) {
			int j = i - di;
			if (j < 0) {
				do {
					j += W;
				} while (j < 0);
			} else if (j >= W) {
				do {
					j -= W;
				} while (j >= W);
			}
			return this->operator()(j);
		};
		return SimdExpression<S, W, decltype(lambda)>(std::move(lambda));
	}
	S max() const {
		S m = (*this)[0];
		for (int w = 1; w < W; w++) {
			m = std::max(m, (*this)[w]);
		}
	}
	S min() const {
		S m = (*this)[0];
		for (int w = 1; w < W; w++) {
			m = std::min(m, (*this)[w]);
		}
	}
	SIMD_BINARY_OPERATOR(+);
	SIMD_BINARY_OPERATOR(-);
	SIMD_BINARY_OPERATOR(*);
	SIMD_BINARY_OPERATOR(/);
	SIMD_BINARY_FUNCTION(max);
	SIMD_BINARY_FUNCTION(min);
	SIMD_BINARY_FUNCTION(pow);
	SIMD_COMPARE_OPERATOR(<);
	SIMD_COMPARE_OPERATOR(==);
	SIMD_COMPARE_OPERATOR(<=);
	SIMD_COMPARE_OPERATOR(>=);
	SIMD_COMPARE_OPERATOR(>);
	SIMD_COMPARE_OPERATOR(!=);
	SIMD_UNARY_FUNCTION(abs);
	SIMD_UNARY_FUNCTION(exp);
	SIMD_UNARY_FUNCTION(log);
	SIMD_UNARY_FUNCTION(sqrt);
	SIMD_UNARY_FUNCTION(cbrt);
	SIMD_UNARY_FUNCTION(sin);
	SIMD_UNARY_FUNCTION(cos);
	SIMD_UNARY_FUNCTION(tan);
	SIMD_UNARY_FUNCTION(asin);
	SIMD_UNARY_FUNCTION(acos);
	SIMD_UNARY_FUNCTION(atan);
	SIMD_UNARY_OPERATOR(+);
	SIMD_UNARY_OPERATOR(-);

private:
	F f_{};
};

#define SIMDVECTOR_ASSIGN_OPERATION(op)                                                                                                    \
	template <typename F_>                                                                                                                 \
	SimdVector &operator op(SimdExpression<S, W, F_> const &other) {                                                                       \
		base_type &self = static_cast<base_type &>(*this);                                                                                 \
		for (int i = 0; i < W; i++) {                                                                                                      \
			self(i) op other(i);                                                                                                           \
		}                                                                                                                                  \
		return *this;                                                                                                                      \
	}

#define SIMDVECTOR_MASKED_ASSIGN_OPERATION(op)                                                                                             \
	template <typename F_>                                                                                                                 \
	Expression &operator op(SimdExpression<S, W, F_> const &other) {                                                                       \
		for (int i = 0; i < W; i++) {                                                                                                      \
			if (mask_(i)) {                                                                                                                \
				self_(i) op other(i);                                                                                                      \
			}                                                                                                                              \
		}                                                                                                                                  \
		return *this;                                                                                                                      \
	}                                                                                                                                      \
	Expression &operator op(S s) {                                                                                                         \
		for (int i = 0; i < W; i++) {                                                                                                      \
			if (mask_(i)) {                                                                                                                \
				self_(i) op s;                                                                                                             \
			}                                                                                                                              \
		}                                                                                                                                  \
		return *this;                                                                                                                      \
	}

template <typename S, size_t W>
struct SimdVector : SimdExpression<S, W, SimdStorage<S, W>> {
	using base_type = SimdExpression<S, W, SimdStorage<S, W>>;
	SimdVector() :
		base_type(SimdStorage<S, W>{}) {
	}
	template <typename F>
	struct Expression {
		Expression(base_type &self, SimdExpression<bool, W, F> const &mask) :
			self_{self}, mask_{mask} {
		}
		SIMDVECTOR_MASKED_ASSIGN_OPERATION(=);
		SIMDVECTOR_MASKED_ASSIGN_OPERATION(+=);
		SIMDVECTOR_MASKED_ASSIGN_OPERATION(-=);
		SIMDVECTOR_MASKED_ASSIGN_OPERATION(*=);
		SIMDVECTOR_MASKED_ASSIGN_OPERATION(/=);

	private:
		base_type &self_;
		SimdExpression<bool, W, F> const &mask_;
	};

	SIMDVECTOR_ASSIGN_OPERATION(=);
	SIMDVECTOR_ASSIGN_OPERATION(+=);
	SIMDVECTOR_ASSIGN_OPERATION(-=);
	SIMDVECTOR_ASSIGN_OPERATION(*=);
	SIMDVECTOR_ASSIGN_OPERATION(/=);

	template <typename F>
	auto operator[](SimdExpression<bool, W, F> const &mask) {
		return Expression(*this, mask);
	}
	S operator[](int i) const {
		base_type const &self = static_cast<base_type const &>(*this);
		return self.operator()(i);
	}
	S &operator[](int i) {
		base_type &self = static_cast<base_type &>(*this);
		return self.operator()(i);
	}
};
//#define BINARY_SIMD_EXPRESSION(name, op)                                                                                                   \
//	template <typename T, SimdExpressionType H1, SimdExpressionType H2>                                                                    \
//	struct Lazy##name##1 {                                                                                                                 \
//		using value_type = T;                                                                                                              \
//		constexpr Lazy##name##1(H1 const &h1, H2 const &h2) :                                                                              \
//			h1_(&h1), h2_(&h2) {                                                                                                             \
//		}                                                                                                                                  \
//		constexpr auto operator[](unsigned i) const {                                                                                      \
//			return (*h1_)[i] op (*h2_)[i];                                                                                                       \
//		}                                                                                                                                  \
//		inline static constexpr unsigned width() {                                                                                         \
//			return H1::width();                                                                                                            \
//		}                                                                                                                                  \
//                                                                                                                                           \
//	private:                                                                                                                               \
//		H1 const *h1_{};                                                                                                                     \
//		H2 const *h2_{};                                                                                                                     \
//	};                                                                                                                                     \
//	template <typename T, SimdExpressionType H>                                                                                            \
//	struct Lazy##name##2 {                                                                                                                 \
//		using value_type = T;                                                                                                              \
//		constexpr Lazy##name##2(H const &h, T const &s) :                                                                                  \
//			h_(&h), s_(s) {                                                                                                                 \
//		}                                                                                                                                  \
//		constexpr auto operator[](unsigned i) const {                                                                                      \
//			return (*h_)[i] op s_;                                                                                                            \
//		}                                                                                                                                  \
//		inline static constexpr unsigned width() {                                                                                         \
//			return H::width();                                                                                                             \
//		}                                                                                                                                  \
//                                                                                                                                           \
//	private:                                                                                                                               \
//		H const *h_{};                                                                                                                       \
//		T s_{};                                                                                                                              \
//	};                                                                                                                                     \
//	template <typename T, SimdExpressionType H>                                                                                            \
//	struct Lazy##name##3 {                                                                                                                 \
//		using value_type = T;                                                                                                              \
//		constexpr Lazy##name##3(T const &s, H const &h) :                                                                                  \
//			h_(&h), s_(s) {                                                                                                                 \
//		}                                                                                                                                  \
//		constexpr auto operator[](unsigned i) const {                                                                                      \
//			return s_ op (*h_)[i];                                                                                                            \
//		}                                                                                                                                  \
//		inline static constexpr unsigned width() {                                                                                         \
//			return H::width();                                                                                                             \
//		}                                                                                                                                  \
//                                                                                                                                           \
//	private:                                                                                                                               \
//		H const *h_{};                                                                                                                       \
//		T s_{};                                                                                                                              \
//	};                                                                                                                                     \
//	template <SimdExpressionType H1, SimdExpressionType H2>                                                                                \
//	constexpr auto operator op(H1 const &h1, H2 const &h2) {                                                                               \
//		using T = typename H1::value_type;                                                                                                 \
//		return Lazy##name##1 < T, H1, H2 > (h1, h2);                                                                                       \
//	}                                                                                                                                      \
//	template <SimdExpressionType H>                                                                                                        \
//	constexpr auto operator op(H const &h, typename H::value_type const &s) {                                                              \
//		return Lazy##name##2 < typename H::value_type, H > (h, s);                                                                         \
//	}                                                                                                                                      \
//	template <SimdExpressionType H>                                                                                                        \
//	constexpr auto operator op(typename H::value_type const &s, H const &h) {                                                              \
//		return Lazy##name##3 < typename H::value_type, H > (s, h);                                                                         \
//	}                                                                                                                                      \
//	template <SimdVectorType V, SimdExpressionType H>                                                                                      \
//	constexpr V &operator op##=(V & v, H const &h) {                                                                                       \
//		constexpr unsigned N = v.width();                                                                                                  \
//		for (unsigned i = 0; i < N; i++) {                                                                                                 \
//			v[i] op## = h[i];                                                                                                              \
//		}                                                                                                                                  \
//	}                                                                                                                                      \
//	template <typename T, SimdExpressionType H1, SimdExpressionType H2>                                                                    \
//		struct IsSimdExpression<Lazy##name##1 < T, H1, H2> > : std::true_type {};                                                          \
//	template <typename T, SimdExpressionType H>                                                                                            \
//		struct IsSimdExpression<Lazy##name##2 < T, H> > : std::true_type {};                                                               \
//	template <typename T, SimdExpressionType H>                                                                                            \
//		struct IsSimdExpression<Lazy##name##3 < T, H> > : std::true_type {};
//
// #define BINARY_SIMD_FUNCTION(name, function) \
//	template <typename T, SimdExpressionType H1, SimdExpressionType H2>                                                                    \
//	struct Lazy##name##1 {                                                                                                                 \
//		using value_type = T;                                                                                                              \
//		constexpr Lazy##name##1(H1 const &h1, H2 const &h2) :                                                                              \
//			h1_(&h1), h2_(&h2) { \
//		}                                                                                                                                  \
//		constexpr auto operator[](unsigned i) const {                                                                                      \
//			return function((*h1_), (*h2_));                                                                                               \
//		}                                                                                                                                  \
//		inline static constexpr unsigned width() {                                                                                         \
//			return H1::width();                                                                                                            \
//		}                                                                                                                                  \
//                                                                                                                                           \
//	private:                                                                                                                               \
//		H1 const *h1_{}; \
//		H2 const *h2_{}; \
//	};                                                                                                                                     \
//	template <typename T, SimdExpressionType H>                                                                                            \
//	struct Lazy##name##2 {                                                                                                                 \
//		using value_type = T;                                                                                                              \
//		constexpr Lazy##name##2(H const &h, T const &s) :                                                                                  \
//			h_(&h), s_(s) { \
//		}                                                                                                                                  \
//		constexpr auto operator[](unsigned i) const {                                                                                      \
//			return function((*h_)[i], s_); \
//		}                                                                                                                                  \
//		inline static constexpr unsigned width() {                                                                                         \
//			return H::width();                                                                                                             \
//		}                                                                                                                                  \
//                                                                                                                                           \
//	private:                                                                                                                               \
//		H const *h_{}; \
//		T s_{}; \
//	};                                                                                                                                     \
//	template <typename T, SimdExpressionType H>                                                                                            \
//	struct Lazy##name##3 {                                                                                                                 \
//		using value_type = T;                                                                                                              \
//		constexpr Lazy##name##3(T const &s, H const &h) :                                                                                  \
//			h_(&h), s_(s) { \
//		}                                                                                                                                  \
//		constexpr auto operator[](unsigned i) const {                                                                                      \
//			return function(s_, (*h_)[i]); \
//		}                                                                                                                                  \
//		inline static constexpr unsigned width() {                                                                                         \
//			return H::width();                                                                                                             \
//		}                                                                                                                                  \
//                                                                                                                                           \
//	private:                                                                                                                               \
//		H const *h_{}; \
//		T s_{}; \
//	};                                                                                                                                     \
//	template <SimdExpressionType H1, SimdExpressionType H2>                                                                                \
//	constexpr auto name(H1 const &h1, H2 const &h2) {                                                                                      \
//		using T = typename H1::value_type;                                                                                                 \
//		return Lazy##name##1 < T, H1, H2 > (h1, h2);                                                                                       \
//	}                                                                                                                                      \
//	template <SimdExpressionType H>                                                                                                        \
//	constexpr auto name(H const &h, typename H::value_type const &s) {                                                                     \
//		return Lazy##name##2 < typename H::value_type, H > (h, s);                                                                         \
//	}                                                                                                                                      \
//	template <SimdExpressionType H>                                                                                                        \
//	constexpr auto name(typename H::value_type const &s, H const &h) {                                                                     \
//		return Lazy##name##3 < typename H::value_type, H > (s, h);                                                                         \
//	}                                                                                                                                      \
//	template <typename T, SimdExpressionType H1, SimdExpressionType H2>                                                                    \
//		struct IsSimdExpression<Lazy##name##1 < T, H1, H2> > : std::true_type {};                                                          \
//	template <typename T, SimdExpressionType H>                                                                                            \
//		struct IsSimdExpression<Lazy##name##2 < T, H> > : std::true_type {};                                                               \
//	template <typename T, SimdExpressionType H>                                                                                            \
//		struct IsSimdExpression<Lazy##name##3 < T, H> > : std::true_type {};
//
// #define UNARY_SIMD_EXPRESSION(name, op) \
//	template <typename T, SimdExpressionType H>                                                                                            \
//	struct Lazy##name {                                                                                                                    \
//		using value_type = T;                                                                                                              \
//		constexpr Lazy##name(H const &h) :                                                                                                 \
//			h_(&h) { \
//		}                                                                                                                                  \
//		constexpr auto operator[](unsigned i) const {                                                                                      \
//			return op (*h_)[i]; \
//		}                                                                                                                                  \
//		inline static constexpr unsigned width() {                                                                                         \
//			return H::width();                                                                                                             \
//		}                                                                                                                                  \
//                                                                                                                                           \
//	private:                                                                                                                               \
//		H const *h_{}; \
//	};                                                                                                                                     \
//	template <SimdExpressionType H>                                                                                                        \
//	constexpr auto operator op(H const &h) {                                                                                               \
//		using T = typename H::value_type;                                                                                                  \
//		return Lazy##name<T, H>(h);                                                                                                        \
//	}                                                                                                                                      \
//	template <typename T, SimdExpressionType H>                                                                                            \
//	struct IsSimdExpression<Lazy##name<T, H>> : std::true_type {};
//
// #define UNARY_SIMD_FUNCTION(name, function) \
//	template <typename T, SimdExpressionType H>                                                                                            \
//	struct Lazy##name {                                                                                                                    \
//		using value_type = T;                                                                                                              \
//		constexpr Lazy##name(H const &h) :                                                                                                 \
//			h_(&h) { \
//		}                                                                                                                                  \
//		constexpr auto operator[](unsigned i) const {                                                                                      \
//			return function((*h_)[i]); \
//		}                                                                                                                                  \
//		inline static constexpr unsigned width() {                                                                                         \
//			return H::width();                                                                                                             \
//		}                                                                                                                                  \
//                                                                                                                                           \
//	private:                                                                                                                               \
//		H const *h_{}; \
//	};                                                                                                                                     \
//	template <SimdExpressionType H>                                                                                                        \
//	constexpr auto name(H const &h) {                                                                                                      \
//		using T = typename H::value_type;                                                                                                  \
//		return Lazy##name<T, H>(h);                                                                                                        \
//	}                                                                                                                                      \
//	template <typename T, SimdExpressionType H>                                                                                            \
//	struct IsSimdExpression<Lazy##name<T, H>> : std::true_type {};
//
// template <typename, unsigned>
// struct SimdVector;
//
// template <typename T>
// struct IsSimdVector : std::false_type {};
//
// template <typename T, unsigned W>
// struct IsSimdVector<SimdVector<T, W>> : std::true_type {};
//
// template <typename T>
// struct IsSimdExpression : std::false_type {};
//
// template <typename T, unsigned W>
// struct IsSimdExpression<SimdVector<T, W>> : std::true_type {};
//
// template <typename T>
// concept SimdVectorType = IsSimdVector<T>::value;
//
// template <typename T>
// concept SimdExpressionType = IsSimdExpression<T>::value;
//
// template <typename, SimdExpressionType, SimdExpressionType, SimdExpressionType>
// struct LazyFMA;
//
// template <typename T, SimdExpressionType H1, SimdExpressionType H2, SimdExpressionType H3>
// struct IsSimdExpression<LazyFMA<T, H1, H2, H3>> : std::true_type {};
//
// template <typename T, unsigned W, typename F>
// struct LazyExpression {
//	using value_type = T;
//	constexpr LazyExpression(F const &f) :
//		f_(f) {
//	}
//	constexpr auto operator[](unsigned i) const {
//		return f_(i);
//	}
//	inline static constexpr unsigned width() {
//		return W;
//	}
//
// private:
//	F const &f_{};
//};
//
// template <typename T, unsigned W>
// struct SimdVector {
//	using value_type = T;
//	inline SimdVector() = default;
//	inline SimdVector(SimdVector const &) = default;
//	inline SimdVector(SimdVector &&) = default;
//	inline SimdVector(T const &init) {
//		vec_.fill(init);
//	};
//	inline SimdVector(T &&init) {
//		vec_.fill(std::move(init));
//	};
//	template <SimdExpressionType Expression>
//	inline SimdVector(Expression const &x) {
//		for (unsigned i = 0; i < W; i++) {
//			vec_[i] = x[i];
//		}
//	};
//	inline SimdVector &operator=(SimdVector const &) = default;
//	inline SimdVector &operator=(SimdVector &&) = default;
//	inline SimdVector &operator=(T const &init) {
//		*this = SimdVector(init);
//		return *this;
//	}
//	inline SimdVector &operator=(T &&init) {
//		*this = SimdVector(std::move(init));
//		return *this;
//	}
//	template <SimdExpressionType Expression>
//	inline SimdVector &operator=(Expression const &x) {
//		for (unsigned i = 0; i < W; i++) {
//			vec_[i] = x[i];
//		}
//		return *this;
//	};
//	inline T operator[](unsigned i) const {
//		return vec_[i];
//	}
//	inline T &operator[](unsigned i) {
//		return vec_[i];
//	}
//	inline T max() const {
//		T a = (*this)[0];
//		for (int i = 1; i < W; i++) {
//			a = std::max(a, (*this)[i]);
//		}
//		return a;
//	}
//	inline T min() const {
//		T a = (*this)[0];
//		for (int i = 1; i < W; i++) {
//			a = std::min(a, (*this)[i]);
//		}
//		return a;
//	}
//	inline auto cshift(int di) const {
//		auto lambda = [this, di](int i) {
//			int j = i - di;
//			if (j < 0) {
//				j += W;
//			} else if (j >= W) {
//				j -= W;
//			}
//			return (*this)[j];
//		};
//		return LazyExpression<T, W, decltype(lambda)>(std::move(lambda));
//	}
//	inline static constexpr unsigned width() {
//		return W;
//	}
//
// private:
//	std::array<T, W> vec_{};
//};
//
// template <typename T, SimdExpressionType H1, SimdExpressionType H2, SimdExpressionType H3>
// struct LazyFMA {
//	using value_type = T;
//	constexpr LazyFMA(H1 const &h1, H2 const &h2, H3 const &h3) :
//		h1_(&h1), h2_(&h2), h3_(&h3) {
//	}
//	constexpr T operator[](unsigned i) const {
//		using std::fma;
//		return fma(h1_[i], h2_[i], h3_[i]);
//	}
//	inline static constexpr unsigned width() {
//		return H1::width();
//	}
//
// private:
//	H1 const *h1_{};
//	H2 const *h2_{};
//	H3 const *h3_{};
//};
//
// template <SimdExpressionType H1, SimdExpressionType H2, SimdExpressionType H3>
// constexpr auto fma(H1 const &h1, H2 const &h2, H3 const &h3) {
//	using T = typename H1::value_type;
//	return LazyFMA<T, H1, H2, H3>(h1, h2, h3);
//}
//
// BINARY_SIMD_EXPRESSION(Addition, +);
// BINARY_SIMD_EXPRESSION(Subtraction, -);
// BINARY_SIMD_EXPRESSION(Multiplication, *);
// BINARY_SIMD_EXPRESSION(Division, /);
// BINARY_SIMD_FUNCTION(pow, pow);
// UNARY_SIMD_EXPRESSION(Plus, +);
// UNARY_SIMD_EXPRESSION(Minus, -);
// UNARY_SIMD_FUNCTION(abs, fabs);
// UNARY_SIMD_FUNCTION(sqrt, sqrt);

#endif /* INCLUDE_SIMDVECTOR_HPP_ */
