#include <config.hpp>
#include <hpx/init.hpp>

#include <numeric>
#include <silo.h>

#include "Gas.hpp"
#include "Interval.hpp"

using ConType = ConservedGasState<cfg::Real, cfg::dimCount>;
using PrimType = PrimitiveGasState<cfg::Real, cfg::dimCount>;

using cfg::dimCount;
using cfg::Real;

constexpr int dimLength = 100;
constexpr int boundWidth = 1;
constexpr int gridSize = dimLength + 2 * boundWidth;
constexpr Real cflFactor = 0.4;
constexpr Real domainBegin = 0.0;
constexpr Real domainEnd = 1.0;
constexpr Real cellWidth = (domainEnd - domainBegin) / dimLength;
constexpr Real maxTime = 1.0;

std::vector<ConType> U(gridSize);
std::vector<ConType> F(gridSize);

void enforceBoundaries() {
	for (unsigned i = 0; i < boundWidth; i++) {
		U[i] = U[U.size() - 2 * boundWidth + i];
		U[U.size() - boundWidth + i] = U[boundWidth + i];
	}
}

Real computeTimestep() {
	Real speed = zero<Real>;
	for (unsigned i = boundWidth; i < U.size() - boundWidth + 1; i++) {
		speed = std::max(speed, maxSignalSpeed(con2prim(U[i]), 0));
	}
	return cflFactor * cellWidth / speed;
}

void computeFluxes() {
	for (unsigned i = boundWidth; i < U.size() - boundWidth + 1; i++) {
		F[i] = con2flux(U[i - 1], U[i], 0);
	}
}

void applyUpdate(Real timeStepSize) {
	auto const λ = timeStepSize / cellWidth;
	for (unsigned i = boundWidth; i < U.size() - boundWidth; i++) {
		U[i] -= λ * (F[i + 1] - F[i]);
	}
}

template <typename T, int D>
struct Silo {
	Silo(std::string filename) :
		domain_{T(1)}, hasMesh{false}, db_{DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5)} {
	}
	void write(std::span<T>) {
		assert(hasMesh);
	}
	void writeMesh(int interiorLength, int ghostWidth) {
		assert(!hasMesh);
		std::array<int, D> dims;
		std::array<std::vector<T>, D> coords;
		std::array<std::string, D> coordnames;
		std::array<void const *, D> coordPointers;
		std::array<char const *, D> namePointers;
		int const exteriorLength = interiorLength + 2 * ghostWidth;
		for (int d = 0; d < D; d++) {
			coords[d].resize(exteriorLength + 1);
			coordnames[d] = std::string(1, 'x' + d);
			dims[d] = exteriorLength + 1;
			auto const dx = domain_.span(d) / interiorLength;
			for (int n = 0; n <= exteriorLength; n++) {
				coords[d][n] = (n - ghostWidth) * dx;
			}
			coordPointers[d] = coords[d].data();
			namePointers[d] = coordnames[d].c_str();
		}
		DBPutQuadmesh(db_, "mesh", namePointers.data(), coordPointers.data(), dims.data(), D, siloDatatype(), DB_COLLINEAR, NULL);
	}
	~Silo() {
		DBClose(db_);
	}

private:
	static constexpr int siloDatatype() {
		if constexpr (std::is_same<T, double>::value) {
			return DB_DOUBLE;
		} else if constexpr (std::is_same<T, float>::value) {
			return DB_FLOAT;
		} else {
			static_assert(false);
		}
	}
	Interval<T, D> domain_{};
	bool hasMesh{};
	DBfile *db_{};
};

void output(std::string name, int number) {
	name += "." + std::to_string(number) + ".silo";
	Silo<cfg::Real, cfg::dimCount> silo(name);
	silo.writeMesh(dimLength, boundWidth);
}

void initialize() {
	for (unsigned i = boundWidth; i < U.size() - boundWidth; i++) {
		PrimType v;
		if (i < U.size() / 2) {
			v.ρ = 1.0;
			v.v = 0.0;
			v.p = 1.0;
		} else {
			v.ρ = 0.1;
			v.v = 0.0;
			v.p = 0.1;
		}
		U[i] = prim2con(v);
	}
}

void test() {
	Silo<cfg::Real, 3> silo("X.silo");
	silo.writeMesh(100, 2);
	return;
	initialize();
	Real currentTime = zero<Real>;
	int stepNumber = 0;
	output("X", 0);
	while (currentTime < maxTime) {
		enforceBoundaries();
		computeFluxes();
		Real const dt = std::min(maxTime - currentTime, computeTimestep());
		applyUpdate(dt);
		printf("i = %i dt = %e\n", stepNumber, dt);
		currentTime += dt;
		stepNumber++;
	}
}

int hpx_main(int argc, char *argv[]) {
	test();
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(22)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
