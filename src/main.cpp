#include "Gas.hpp"
#include <config.hpp>
#include <hpx/init.hpp>

#include <numeric>
#include <silo.h>
#include <string>
#include <vector>

#include "Gas.hpp"
#include "Interval.hpp"

using namespace config;

constexpr int interiorLength = 256;
constexpr int ghostWidth = 1;
constexpr int exteriorLength = interiorLength + 2 * ghostWidth;
constexpr Real cflFactor = 0.4;
constexpr Real domainBegin = 0.0;
constexpr Real domainEnd = 1.0;
constexpr Real cellWidth = (domainEnd - domainBegin) / interiorLength;
constexpr Real timeEnd = 1.0;

using ConType = ConservedGasState;
using PrimType = PrimitiveGasState;

Gas::AoS U(exteriorLength);
Gas::AoS F(exteriorLength);

void enforceBoundaries() {
	for (unsigned i = 0; i < ghostWidth; i++) {
		U[i] = U[U.size() - 2 * ghostWidth + i];
		U[U.size() - ghostWidth + i] = U[ghostWidth + i];
	}
}

Real computeTimestep() {
	Real speed = zero<Real>;
	for (unsigned i = ghostWidth; i < U.size() - ghostWidth + 1; i++) {
		speed = std::max(speed, maxSignalSpeed(con2prim(U[i]), 0));
	}
	return cflFactor * cellWidth / speed;
}

void computeFluxes() {
	for (unsigned i = ghostWidth; i < U.size() - ghostWidth + 1; i++) {
		F[i] = con2flux(U[i - 1], U[i], 0);
	}
}

void applyUpdate(Real timeStepSize) {
	auto const λ = timeStepSize / cellWidth;
	for (unsigned i = ghostWidth; i < U.size() - ghostWidth; i++) {
		U[i] -= λ * (F[i + 1] - F[i]);
	}
}

struct Silo {
	static constexpr int siloDimCount = std::max(dimCount, 2);
	Silo(std::string filename) :
		domain_(Real(1)), hasMesh_{false}, db_{DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5)} {
	}
	void writeVariable(std::string const &name, Real const *var) {
		assert(hasMesh_);
		int const exteriorLength = interiorLength_ + 2 * ghostWidth_;
		std::array<int, siloDimCount> dims;
		for (int d = 0; d < dimCount; d++) {
			dims[d] = exteriorLength;
		}
		for (int d = dimCount; d < siloDimCount; d++) {
			dims[d] = 1;
		}
		DBPutQuadvar1(db_, name.c_str(), "mesh", var, dims.data(), siloDimCount, NULL, 0, siloDatatype(), DB_ZONECENT, NULL);
	}
	void writeMesh(int interiorLength, int ghostWidth) {
		assert(!hasMesh_);
		interiorLength_ = interiorLength;
		ghostWidth_ = ghostWidth;
		std::array<int, siloDimCount> dims;
		std::array<std::vector<Real>, siloDimCount> coords;
		std::array<std::string, siloDimCount> coordnames;
		std::array<void const *, siloDimCount> coordPointers;
		std::array<char const *, siloDimCount> namePointers;
		int const exteriorLength = interiorLength + 2 * ghostWidth;
		for (int d = 0; d < dimCount; d++) {
			coords[d].resize(exteriorLength + 1);
			coordnames[d] = std::string(1, 'x' + d);
			dims[d] = exteriorLength + 1;
			auto const dx = domain_.span(d) / interiorLength;
			for (int n = 0; n <= exteriorLength; n++) {
				coords[d][n] = (n - ghostWidth) * dx;
			}
		}
		for (int d = dimCount; d < siloDimCount; d++) {
			coords[d].resize(2);
			coordnames[d] = std::string(1, '_' + d);
			dims[d] = 2;
			coords[d][0] = Real(0.0);
			coords[d][1] = Real(1.0);
		}
		for (int d = 0; d < siloDimCount; d++) {
			coordPointers[d] = coords[d].data();
			namePointers[d] = coordnames[d].c_str();
		}
		int const dataType = siloDatatype();
		DBPutQuadmesh(db_, "mesh", namePointers.data(), coordPointers.data(), dims.data(), siloDimCount, dataType, DB_COLLINEAR, NULL);
		hasMesh_ = true;
	}
	~Silo() {
		DBClose(db_);
	}

private:
	static constexpr int siloDatatype() {
		if constexpr (std::is_same<Real, double>::value) {
			return DB_DOUBLE;
		} else if constexpr (std::is_same<Real, float>::value) {
			return DB_FLOAT;
		} else {
			assert(false);
		}
	}
	Interval<Real, dimCount> domain_{};
	bool hasMesh_{};
	int interiorLength_;
	int ghostWidth_;
	DBfile *db_{};
};

void output(std::string name, int number) {
	name += "." + std::to_string(number) + ".silo";
	Silo silo(name);
	silo.writeMesh(interiorLength, ghostWidth);
	auto const data = aos2soa(U);
	silo.writeVariable("rho", &(data.ρ[0]));
	silo.writeVariable("tau", &(data.τ[0]));
	silo.writeVariable("E", &(data.E[0]));
	for (int d = 0; d < dimCount; d++) {
		std::string const nm = "S" + std::string(1, 'x' + d);
		silo.writeVariable(nm.c_str(), &(data.s[d][0]));
	}
}

void initialize() {
	for (unsigned i = ghostWidth; i < U.size() - ghostWidth; i++) {
		PrimType v;
		if (i < U.size() / 2) {
			v.ρ = 1.0;
			v.v = 0.0;
			v.p = 1.0;
		} else {
			v.ρ = 0.125;
			v.v = 0.0;
			v.p = 0.1;
		}
		U[i] = prim2con(v);
	}
}

void test() {
	initialize();
	Real timeNow = zero<Real>;
	int stepNumber = 0;
	while (timeNow < timeEnd) {
		output("X", stepNumber);
		enforceBoundaries();
		computeFluxes();
		Real const dt = std::min(timeEnd - timeNow, computeTimestep());
		applyUpdate(dt);
		printf("i = %i dt = %e\n", stepNumber, dt);
		timeNow += dt;
		stepNumber++;
	}
	output("X", stepNumber);
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
