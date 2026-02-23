#include "Grid.hpp"
#include "Silo.hpp"

#include <hpx/init.hpp>

constexpr Real cfl = 0.4_R;

void driver(Real tmax) {
	constexpr Integer D = 1;
	constexpr Integer O = 2;
	constexpr Integer N = 128;
	Grid<O, D, N> grid;
	grid.initialize([](auto const &x) {
		State<D> u{};
		if (x[0] < 0.5_R) {
			u.setDensity(1_R);
			u.setPressure(1_R);
		} else {
			u.setDensity(0.125_R);
			u.setPressure(0.1_R);
		}
		return u;
	});
	constexpr Real beta[3][2] = {{}, {1_R}, {1_R, 0.5_R}};
	constexpr Real dx = grid.cellWidth;
	Real t = 0_R;
	Real a, dt;
	Integer iter = 0_I;
	auto const output = [&]() {
		printf("i=%i t=%e dt=%e\n", int(iter), t, dt);
		Vector<D> const origin(-2_R / Real(N));
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
	driver(1_R / 8_R);
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(22)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
