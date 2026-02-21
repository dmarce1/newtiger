#include "Grid.hpp"

#include <hpx/init.hpp>

constexpr Real cfl = 0.4_R;

void driver(Real tmax) {
	Real t = 0_R;
	constexpr Integer D = 1;
	constexpr Integer O = 2;
	constexpr Integer N = 128;
	Grid<D, N> grid;
	grid.initialize([](std::array<Integer, D> const &idx) {
		State<D> u{};
		if (2 * idx[0] < N) {
			u.setDensity(1_R);
			u.setPressure(1_R);
		} else {
			u.setDensity(0.125_R);
			u.setPressure(0.1_R);
		}
		return u;
	});
	constexpr std::array<Real, O> beta = {1_R, 0.5_R};
	constexpr Real dx = grid.cellWidth;
	Real a, dt;
	while (std::nextafter(t, tmax) < tmax) {
		for (Integer rk = 0; rk < O; rk++) {
			grid.reconstruct();
			a = grid.fluxes();
			dt = (rk == 0) ? (cfl * dx / a) : dt;
			grid.update(dt, beta[rk]);
			grid.boundaries();
		}
		grid.store();
		t += dt;
	}
}

int hpx_main(int argc, char *argv[]) {
	driver(1.0_R);
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(22)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
