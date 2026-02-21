/*
 * FourierLegendre.hpp
 *
 *  Created on: Feb 7, 2026
 *      Author: dmarce1
 */

#ifndef BASIS_HPP_
#define BASIS_HPP_

#include "Matrix.hpp"
#include "Quadrature.hpp"
#include <concepts>
#include <iostream>
#include <span>
#include <utility>

template <typename, int, int>
struct AnalysisVector;

template <typename T, int N>
struct AnalysisVector<T, N, 1> : public Vector<T, N> {};

template <typename T, int N, int D>
struct AnalysisVector {

	template <std::integral auto I>
	auto &get() {
		if constexpr (I == 0) {
			return this_;
		} else {
			return next_.template get<I - 1>();
		}
	}

	template <std::integral auto I>
	auto const &get() const {
		if constexpr (I == 0) {
			return this_;
		} else {
			return next_.template get<I - 1>();
		}
	}

	template <std::integral Int>
	T &operator()(std::array<Int, D> nml) {
		if (nml[0] == 0) {
			return this_(dimReduce(nml));
		} else if constexpr (N > 1) {
			nml[0]--;
			return next_(nml);
		} else {
			assert(false);
		}
	}

	template <std::integral Int>
	T operator()(std::array<Int, D> nml) const {
		return const_cast<AnalysisVector const *>(this)->operator()(nml);
	}

private:
	template <std::integral Int>
	static constexpr std::array<Int, D - 1> dimReduce(std::array<Int, D> const &idx1) {
		constexpr int Dm1 = D - 1;
		std::array<Int, D - 1> idx2;
		for (int d = 0; d < Dm1; d++) {
			idx2[d] = idx1[d + 1];
		}
		return idx2;
	}
	struct NullType {};
	std::conditional_t<(D > 1), AnalysisVector<T, N, D - 1>, T> this_;
	std::conditional_t<(N > 1), AnalysisVector<T, N - 1, D>, NullType> next_;
};

template <typename T, int N, int D>
using SynthesisVector = Vector<T, N, D>;

template <typename T, int N, int D>
constexpr AnalysisVector<T, N, D> fma(Scalar auto a, AnalysisVector<T, N, D> const &B, AnalysisVector<T, N, D> const &C) {
	auto const lambda = [a, B]<size_t... I>(AnalysisVector<T, N, D> C, std::index_sequence<I...>) {
		(([a, B](AnalysisVector<T, N, D> &C) {
			 C.template get<I>() += a * B.template get<I>();
		 }),
		 ...);
	};
	lambda(C, std::make_index_sequence<N>());
	return C;
}
/*
----------------------------------------------
A    α                 Α                 alpha
B    β                 Β                 beta
G    γ                 Γ                 gamma
D    δ                 Δ                 delta
E    ε                 Ε                 epsilon
Z    ζ                 Ζ                 zeta
H    η                 Η                 eta
Q    θ                 Θ                 theta
I    ι                 Ι                 iota
K    κ                 Κ                 kappa
L    λ                 Λ                 lambda
M    μ                 Μ                 mu
N    ν                 Ν                 nu
C    ξ                 Ξ                 xi
O    ο                 Ο                 omicron
P    π                 Π                 pi
R    ρ                 Ρ                 rho
S    σ / ς             Σ                 sigma
T    τ                 Τ                 tau
U    υ                 Υ                 upsilon
F    φ                 Φ                 phi
X    χ                 Χ                 chi
Y    ψ                 Ψ                 psi
W    ω                 Ω                 omega
----------------------------------------------
₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ
----------------------------------------------
⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ
*/

enum class BasisType : int { legendre, taylor };
enum class QuatratureType : int { legendre, lobatto };

template <typename Type, int modeCount, QuatratureType quadratureType = QuatratureType::legendre, BasisType testType = BasisType::legendre,
		  BasisType trialType = BasisType::legendre>
struct Basis {
	constexpr static int nodeCount = (quadratureType == QuatratureType::legendre) ? modeCount : (modeCount + 1);
	constexpr Basis() {
		constexpr auto x = quadrature().positions();
		constexpr auto w = quadrature().weights();
		for (int n = 0; n < modeCount; n++) {
			for (int m = 0; m < nodeCount; m++) {
				test[n][m] = w[m] * testFunction(n, x[m]);
				dTestDx[n][m] = w[m] * testDerivative(n, x[m]);
				trial[m][n] = trialFunction(n, x[m]);
			}
		}
		inverseMass = inv(test * trial);
	}
	template <int dimensionCount, int derivativeOrder>
	AnalysisVector<Type, modeCount, dimensionCount> analyze(SynthesisVector<Type, nodeCount, dimensionCount> const &x) {
		auto const &A = derivativeOrder ? test : testDerivative;
		if constexpr (dimensionCount == 1) {
			AnalysisVector<Type, nodeCount, 1> y;
			for (int n = 0; n < modeCount; n++) {
				for (int m = 0; m < nodeCount; m++) {
					y[n] = fma(A[n][m], x[m], y[n]);
				}
			}
			return y;
		} else if constexpr (dimensionCount >= 2) {
			AnalysisVector<Type, nodeCount, dimensionCount> z;
			std::array<SynthesisVector<Type, nodeCount, dimensionCount - 1>, modeCount> y;
			for (int n = 0; n < modeCount; n++) {
				for (int m = 0; m < nodeCount; m++) {
					y[n] = fma(A[n][m], x[m], y[n]);
				}
			}
			auto const nextDimension = [&y, &z]<int... I>(std::integer_sequence<int, I...>) {
				(([&y, &z]() {
					 z.template get<I>() = F(y[I]);
				 }()),
				 ...);
			};
			nextDimension(std::make_index_sequence<modeCount>{});
			return z;
		}
	}
	template <int dimensionCount>
	SynthesisVector<Type, modeCount, dimensionCount> synthesize(AnalysisVector<Type, nodeCount, dimensionCount> const &x) {
		auto const &A = trial;
		if constexpr (dimensionCount == 1) {
			SynthesisVector<Type, nodeCount, 1> y;
			for (int n = 0; n < modeCount; n++) {
				for (int m = 0; m < nodeCount; m++) {
					y[n] = fma(A[n][m], x[m], y[n]);
				}
			}
			return y;
		} else if constexpr (dimensionCount >= 2) {
			SynthesisVector<Type, nodeCount, dimensionCount> z;
			std::array<SynthesisVector<Type, nodeCount, dimensionCount - 1>, modeCount> y;
			auto const nextDimension = [&y, &x]<int... I>(std::integer_sequence<int, I...>) {
				(([&y, &x]() {
					 y[I] = F(x.template get<I>());
				 }()),
				 ...);
			};
			nextDimension(std::make_index_sequence<nodeCount>{});
			for (int m = 0; m < nodeCount; m++) {
				for (int n = 0; n < modeCount; n++) {
					z[m] = fma(A[m][n], y[n], z[m]);
				}
			}
			return z;
		}
	}

private:
	constexpr static auto quadrature() {
		if constexpr (quadratureType == QuatratureType::legendre) {
			return GaussLobattoQuadrature<Type, nodeCount>{};
		} else if constexpr (quadratureType == QuatratureType::lobatto) {
			return GaussLegendreQuadrature<Type, nodeCount>{};
		}
	}
	template <BasisType basisType>
	constexpr static Type basisFunction(int n, Type x) {
		if constexpr (testType == BasisType::legendre) {
			return legendreP(n, x);
		} else if constexpr (testType == BasisType::taylor) {
			return pow(x, n);
		}
	}
	template <BasisType basisType>
	constexpr static Type basisDerivative(int n, Type x) {
		if constexpr (testType == BasisType::legendre) {
			return dLegendrePdx<1>(n, x);
		} else if constexpr (testType == BasisType::taylor) {
			return Type(n) * pow(x, n - 1);
		}
	}
	constexpr static Type testFunction(int n, Type x) {
		return basisFunction<testType>(n, x);
	}
	constexpr static Type testDerivative(int n, Type x) {
		return basisDerivative<testType>(n, x);
	}
	constexpr static Type trialFunction(int n, Type x) {
		return basisFunction<trialType>(n, x);
	}
	Matrix<Type, nodeCount, nodeCount> inverseMass;
	Matrix<Type, nodeCount, modeCount> dTestDx;
	Matrix<Type, modeCount, nodeCount> test;
	Matrix<Type, nodeCount, modeCount> trial;
};

// template <int modeCount, typename Test, typename Trial, Quadrature Quadrature, int nodeCount = Quadrature::size()>
// auto createBasis(Test const &testFunction, Trial const &trialFunction) {
//	using Type = double;
//	auto const x = q.positions();
//	auto const w = q.weights();
//	Matrix<Type, modeCount, nodeCount> ψₙₘ{};
//	Matrix<Type, nodeCount, modeCount> φₘₖ{};
//	Matrix<Type, modeCount, nodeCount> dₓψₙₘ{};
//	for (int n = 0; n < modeCount; n++) {
//		for (int m = 0; m < nodeCount; m++) {
//			ψₙₘ[n][m] = w[m] * v(n, x[m]);
//			φₘₖ[m][n] = u(n, x[m]);
//			dₓψₙₘ[n][m] = w[m] * v(n, x[m], 1);
//		}
//	}
//	auto const iMₙₖ = inv(ψₙₘ * φₘₖ);
//	auto const Sₙₖ = dₓψₙₘ * φₘₖ;
// }
//  template <typename T, int N, int D, Quadrature Q>
//  constexpr auto createFourierLegendreTransform(Q) {
//	static_assert(Q::exact2degree() >= 2 * N - 2);
//	using std::fma;
//	constexpr Q q{};
//	constexpr int M = Q::size();
//	constexpr auto x = q.positions();
//	constexpr auto w = q.weights();
//	using Analysis = AnalysisVector<T, N, D>;
//	using Synthesis = SynthesisVector<T, M, D>;
//	Matrix<T, N, M> A{};
//	for (int n = 0; n < N; n++) {
//		for (int m = 0; m < M; m++) {
//			A[n][m] = T(2 * n + 1) / T(2) * legendreP(n, x[m]) * w[m];
//		}
//	}
//	if constexpr (D == 1) {
//		return [A](Synthesis const &fx) {
//			Analysis fn{};
//			for (int n = 0; n < N; n++) {
//				for (int m = 0; m < M; m++) {
//					fn[n] = fma(A[n][m], fx[m], fn[n]);
//				}
//			}
//			return fn;
//		};
//	} else {
//		return [A](Synthesis const &fxy) {
//			constexpr std::make_integer_sequence<int, N> seq{};
//			Analysis fnk;
//			std::array<SynthesisVector<T, M, D - 1>, N> fnx{};
//			for (int n = 0; n < N; n++) {
//				for (int m = 0; m < M; m++) {
//					fnx[n] = fma(A[n][m], fxy[m], fnx[n]);
//				}
//			}
//			auto const nextDimension = [&fnk, &fnx]<int... I>(std::integer_sequence<int, I...>) {
//				(([&fnk, &fnx]() {
//					 constexpr auto F = createFourierLegendreTransform<T, N - I, D - 1, Q>(q);
//					 fnk.template get<I>() = F(fnx[I]);
//				 }()),
//				 ...);
//			};
//			nextDimension(seq);
//			return fnk;
//		};
//	}
//  }
//
//  template <typename T, int N, int D, Quadrature Q>
//  constexpr auto createInverseFourierLegendreTransform(Q) {
//	static_assert(Q::exact2degree() >= 2 * N - 2);
//	constexpr Q q{};
//	constexpr int M = Q::size();
//	constexpr auto x = q.positions();
//	using Analysis = AnalysisVector<T, N, D>;
//	using Synthesis = SynthesisVector<T, M, D>;
//	Matrix<T, M, N> S{};
//	for (int m = 0; m < M; m++) {
//		for (int n = 0; n < N; n++) {
//			S[m][n] = legendreP(n, x[m]);
//		}
//	}
//	if constexpr (D == 1) {
//		return [S](Analysis const &fn) {
//			Synthesis fx{};
//			for (int m = 0; m < M; m++) {
//				for (int n = 0; n < N; n++) {
//					fx[m] = fma(S[m][n], fn[n], fx[m]);
//				}
//			}
//			return fx;
//		};
//	} else {
//		return [S](Analysis const &fnk) {
//			constexpr std::make_integer_sequence<int, N> seq{};
//			Synthesis fxy{};
//			std::array<SynthesisVector<T, M, D - 1>, N> fny;
//			auto const nextDimension = [&fny, &fnk]<int... I>(std::integer_sequence<int, I...>) {
//				(([&fny, &fnk]() {
//					 constexpr auto iF = createInverseFourierLegendreTransform<T, N - I, D - 1, Q>(q);
//					 fny[I] = iF(fnk.template get<I>());
//				 }()),
//				 ...);
//			};
//			nextDimension(seq);
//			for (int m = 0; m < M; m++) {
//				for (int n = 0; n < N; n++) {
//					fxy[m] = fma(S[m][n], fny[n], fxy[m]);
//				}
//			}
//			return fxy;
//		};
//	}
//  }
//
//   template <std::integral Int, std::integral auto B, std::integral auto E, std::integral auto... Is>
//	struct MakeIntegerSequence
//	: public std::conditional_t < B<E, MakeIntegerSequence<Int, B, E - 1, E, Is...>, MakeIntegerSequence<Int, B + 1, E, B, Is...>> {};
//
//   template <std::integral Int, std::integral auto I0, std::integral auto... Is>
//   struct MakeIntegerSequence<Int, I0, I0, Is...> {
//	using type = std::integer_sequence<Int, I0, Is...>;
//   };

#endif /* BASIS_HPP_ */
