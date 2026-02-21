
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

// enum class ForLoopShape : int { square, triangular };
//
// template <int D, ForLoopShape S = ForLoopShape::square, typename F>
// constexpr void mdForLoop(int n, F const &f) {
//	Vector<int, D> idx;
//	auto const lambda = [&]<int depth>(auto const &self) {
//		int &i = idx[depth];
//		for (i = 0; i < n; i++) {
//			if constexpr (depth == 0) {
//				f(idx);
//			} else {
//				if constexpr (S == ForLoopShape::triangular) {
//					n += i;
//				}
//				self.template operator()<depth - 1>(self);
//				if constexpr (S == ForLoopShape::triangular) {
//					n -= i;
//				}
//			}
//		}
//	};
//	lambda.template operator()<D - 1>(lambda);
// }
//template <int N, int D, int I>
//constexpr void swapDimensions(auto *u) {
//	constexpr int size = pow(N, D);
//	if constexpr (I > 0) {
//		constexpr int M = pow(N, I - 1);
//		constexpr int L = pow(N, D - I - 1);
//		constexpr int NM = N * M;
//		constexpr int NMm1 = N * M - 1;
//		constexpr int NMp1 = N * M + 1;
//		for (int l = 0; l < L; l++) {
//			for (int di = 1; di < N; di++) {
//				int const j0 = NM * (N * l + di);
//				for (int m = 0; m < M; m++) {
//					int const j1 = j0 + N * m;
//					for (int i = 0; i + di < N; i++) {
//						int const k1 = j1 + NMp1 * i;
//						int const k2 = k1 - NMm1 * di;
//						assert(k1 < size);
//						assert(k2 < size);
//						std::swap(u[k1], u[k2]);
//					}
//				}
//			}
//		}
//	}
//}
//
