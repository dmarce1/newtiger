#include "Constants.hpp"
#include "CoupledState.hpp"
#include "Grid.hpp"
#include "Integer.hpp"
#include "RadiationState.hpp"
#include "Rational.hpp"
#include "Silo.hpp"
#include "Simd.hpp"
#include "Vector.hpp"

#include <hpx/init.hpp>

constexpr Real cfl = 0.4_R / 3_R;

void driver(Real tmax) {
	constexpr Integer D = 2;
	constexpr Integer O = 3;
	constexpr Integer N = 128;
	Grid<O, D, N, CoupledState> grid;
	//	grid.initialize([](auto const &X) {
	//		GasState<Real, D> u;
	//		auto x = inv(sqrt(2_R)) * (X[0] - X[1]);
	//		std::fill(u.begin(), u.end(), 0_R);
	//		if (x < 0.0_R) {
	//			u.setDensity(1_R);
	//			u.setPressure(1_R);
	//		} else {
	//			u.setDensity(0.125_R);
	//			u.setPressure(0.1_R);
	//		}
	//		u.setEntropy();
	//		return u;
	//	});
	//	grid.initialize([N](auto const &X) {
	//		GasState<Real, D> u;
	//		auto x = abs(X);
	//		std::fill(u.begin(), u.end(), 0_R);
	//		auto const y = exp(-x * N * inv(3_R));
	//		u.setDensity(y);
	//		u.setPressure(y);
	//		u.setEntropy();
	//		return u;
	//	});
	grid.initialize([N](auto const &X) {
		GasState<Real, D> um{};
		RadiationState<Real, D> ur{};
		auto x = abs(X);
		std::fill(um.begin(), um.end(), 0_R);
		std::fill(ur.begin(), ur.end(), 0_R);
		um.setDensity(1_R);
		um.setPressure(1_R);
		um.setEnergy();
		auto const Er = x > 0.1_R ? 0.1_R : 1_R;
		auto const Fr = 0.99_R * Er * Vector<Real, D>::unit(0);
		ur.setEnergy(Er);
		ur.setFlux(Fr);
		CoupledState<Real, D> u;
		u.matter() = um;
		u.radiation() = ur;
		return u;
	});
	constexpr Real beta[4][3] = {{}, {1_R}, {1_R, 0.5_R}, {1_R, 0.25_R, 2_R / 3_R}};
	constexpr Real dx = grid.cellWidth;
	Real t = 0_R;
	Real a, dt;
	Integer iter = 0_I;
	auto const output = [&]() {
		printf("i=%i t=%e dt=%e\n", int(iter), t, dt);
		Vector<Real, D> const origin(-2_R / Real(N));
		std::string const fname = "X." + std::to_string(iter) + ".silo";
		Silo<D> silo(fname);
		grid.output(silo);
	};
	while (std::nextafter(t, tmax) < tmax) {
		for (Integer rk = 0; rk < O; rk++) {
			grid.reconstruct();
			a = grid.fluxes();
			if (rk == 0) {
				dt = std::min(cfl * dx / a, tmax - t);
				output();
			}
			grid.update(dt, beta[O][rk]);
			grid.boundaries();
		}
		grid.store();
		t += dt;
		iter++;
	}
	output();
}

#include <iostream>
int hpx_main(int argc, char *argv[]) {
	driver(1_R / 8_R);
	//	auto const genLegendre = []<Integer L>() {
	//		Vector<Rational, L + 1> cn{}, cnp1{}, cnm1{};
	//		Vector<Rational, L / 2 + 1> cn2{};
	//		cn = Vector<Rational, N>::unit(0);
	//		for (Integer n = 0; n < L; n++) {
	//			cnp1 = Rational(2_I * n + 1_I, n + 1_I) * cn.shiftHi() - Rational(n, n + 1_I) * cnm1;
	//			cnm1 = cn;
	//			cn = cnp1;
	//		}
	//		Integer const odd = ~L & 1_I;
	//		for (Integer n = odd; n <= l; n += 2) {
	//			cn2[n >> 1] = cn[n];
	//		}
	//		return [cn2](Real x) {
	//			Real x2 = sqr(x);
	//			Real y = cn.back();
	//			for (Integer k = cn.size() - 2_I; k >= 0_I; k--) {
	//				y = x * y + cn[k];
	//			}
	//			if (odd) y *= x;
	//		};
	//	};
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(22)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
