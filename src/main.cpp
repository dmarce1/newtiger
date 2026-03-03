#include "AutoDiff.hpp"
#include "Constants.hpp"
#include "CoupledState.hpp"
#include "Grid.hpp"
#include "Integer.hpp"
#include "RadiationState.hpp"
#include "Silo.hpp"
#include "Simd.hpp"

#include <hpx/init.hpp>

constexpr Real cfl = 0.4_R;

void driver(Real tmax) {
	constexpr Integer D = 2;
	constexpr Integer O = 1;
	constexpr Integer N = 64;
	Grid<O, D, N, GasState> grid;
	grid.initialize([](auto const &x) {
		GasState<Real, D> u;
		std::fill(u.begin(), u.end(), 0_R);
		if (std::abs(x[0] - 0.5) < 0.25_R) {
			u.setDensity(1_R);
			u.setPressure(1_R);
		} else {
			u.setDensity(0.125_R);
			u.setPressure(0.1_R);
		}
		u.setEntropy();
		return u;
	});
	constexpr Real beta[3][2] = {{}, {1_R}, {1_R, 0.5_R}};
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

int hpx_main(int argc, char *argv[]) {
	//	std::cout << unitConversion(1_cm, 1_g, 1_s, 1_K) << std::endl;
	//	std::cout << unitConversion(Constants::c, Constants::G, Constants::σ, Constants::kB / Constants::mᵤ) << std::endl;
	driver(1_R / 8_R);

//	constexpr Integer D = 2;
//	RadiationState<Real, D> rad;
//	GasState<Real, D> gas;
//	gas.ρ = 1_R;
//	rad.E = 3.4_R;
//	gas.e = 1_R;
//	rad.F = 0_R;
//	gas.m = 0_R;
//	rad.F[0] = 0.99_R * rad.E;
//	gas.setInternalEnergy();
//	CoupledState<Real, D> coupled(gas, rad);
//	coupled.solveImplicit(1_R, 1_R, 1_R, 1_R);
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(22)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
