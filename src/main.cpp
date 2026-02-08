
#include "Constants.hpp"
#include "FourierLegendre.hpp"
#include "IO.hpp"
#include "Legendre.hpp"
#include "Math.hpp"
#include "Matrix.hpp"
#include "Units.hpp"
#include "Vector.hpp"
#include "config.hpp"
#include <iostream>
#include <iterator>
#include <span>

template <int N>
constexpr void foo() {
}

int main(int argc, char **argv) {
	using T = double;
	constexpr int N = 4;
//	constexpr int M = N ;
//	constexpr GaussLegendreQuadrature<double, M> q{};
	constexpr int M = N + 1;
	constexpr GaussLobattoQuadrature<double, M> q{};
	std::cout << q;
	constexpr int D = 4;
	using VecA = AnalysisVector<T, N, D>;
	using VecS = SynthesisVector<T, M, D>;
	constexpr auto F = createFourierLegendreTransform<T, N, D>(q);
	constexpr auto iF = createInverseFourierLegendreTransform<T, N, D>(q);
	VecS fxyz{};
	VecA fnkl{};
	for (int n = 0; n < pow(M, D); n++) {
		std::array<int, D> i;
		int m = n;
		for (int d = 0; d < D; d++) {
			i[d] = m % M;
			printf("%c = %i ", 'x' + d, i[d]);
			m /= M;
		}
//		fxyz(i) = 1;
		fxyz(i) = (rand() & 0xf) - 0x7;
		printf(" %e\n", fxyz(i));
	}
	printf("\n");
	fnkl = F(fxyz);
	fxyz = iF(fnkl);
	fnkl = F(fxyz);
	auto fxyz0 = iF(fnkl);
	for (int n = 0; n < pow(N, D); n++) {
		std::array<int, D> i;
		int m = n;
		for (int d = 0; d < D; d++) {
			i[d] = m % N;
			m /= N;
		}
		if (std::accumulate(i.begin(), i.end(), 0) >= N) continue;
		for (int d = 0; d < D; d++) {
			printf("%c = %i ", 'x' + d, i[d]);
		}
		printf(" %e\n", fnkl(i));
	}
	printf("\n");
	double l2 = 0.0;
	for (int n = 0; n < pow(M, D); n++) {
		std::array<int, D> i;
		int m = n;
		for (int d = 0; d < D; d++) {
			i[d] = m % M;
			printf("%c = %i ", 'x' + d, i[d]);
			m /= M;
		}
		auto const err = fxyz0(i) - fxyz(i);
		l2 += sqr(err);
		printf(" %e\n", fxyz0(i) - fxyz(i));
	}
	l2 = sqrt(l2);
	printf("\n");
	printf("err %e\n", l2);
	return 0;
}
