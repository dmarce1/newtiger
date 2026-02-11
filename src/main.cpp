#include <config.hpp>
#include <hpx/init.hpp>
// #include "Gas.hpp"
#include "SimdVector.hpp"
//
// #include <numeric>
// #include <silo.h>
// #include <string>
// #include <vector>
//
// #include "Gas.hpp"
// #include "Interval.hpp"
//
// using namespace config;
//
// constexpr int interiorLength = 1024;
// constexpr int ghostWidth = 1;
// constexpr int exteriorLength = interiorLength + 2 * ghostWidth;
// constexpr Real cflFactor = 0.4;
// constexpr Real domainBegin = 0.0;
// constexpr Real domainEnd = 1.0;
// constexpr Real cellWidth = (domainEnd - domainBegin) / interiorLength;
// constexpr Real timeEnd = 0.125;
//
// using ConType = Gas::State;
//
// Gas::SoA<exteriorLength> Usoa;
// Gas::SoA<exteriorLength> Fsoa;
//
// void enforceBoundaries() {
//	auto Uaos = soa2aos(Usoa);
//	for (unsigned i = 0; i < ghostWidth; i++) {
//		Uaos[i] = Uaos[exteriorLength - 2 * ghostWidth + i];
//		Uaos[exteriorLength - ghostWidth + i] = Uaos[ghostWidth + i];
//	}
//	Usoa = aos2soa(Uaos);
//}
//
// Real computeTimestep() {
//	return cflFactor * cellWidth / maxSignalSpeed(Usoa, 0).max();
//}
//
////void computeFluxes() {
////	auto const &Ul = Usoa;
////	auto const Ur = Usoa.cshift(1);
////	Fsoa = con2flux(Ul, Ur, 0);
////}
////
////void applyUpdate(Real timeStepSize) {
////	auto const λ = timeStepSize / cellWidth;
////	auto const &Fl = Usoa;
////	auto const Fr = Fsoa.cshift(-1);
////	Usoa -= λ * (Fr - Fl);
////	dualEnergyUpdate(Usoa);
////}
////
////struct Silo {
////	static constexpr int siloDimCount = std::max(dimCount, 2);
////	Silo(std::string filename) :
////		domain_(Real(1)), hasMesh_{false}, db_{DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5)} {
////	}
////	void writeVariable(std::string const &name, Real const *var) {
////		assert(hasMesh_);
////		int const exteriorLength = interiorLength_ + 2 * ghostWidth_;
////		std::array<int, siloDimCount> dims;
////		for (int d = 0; d < dimCount; d++) {
////			dims[d] = exteriorLength;
////		}
////		for (int d = dimCount; d < siloDimCount; d++) {
////			dims[d] = 1;
////		}
////		DBPutQuadvar1(db_, name.c_str(), "mesh", var, dims.data(), siloDimCount, NULL, 0, siloDatatype(), DB_ZONECENT, NULL);
////	}
////	void writeMesh(int interiorLength, int ghostWidth) {
////		assert(!hasMesh_);
////		interiorLength_ = interiorLength;
////		ghostWidth_ = ghostWidth;
////		std::array<int, siloDimCount> dims;
////		std::array<std::vector<Real>, siloDimCount> coords;
////		std::array<std::string, siloDimCount> coordnames;
////		std::array<void const *, siloDimCount> coordPointers;
////		std::array<char const *, siloDimCount> namePointers;
////		int const exteriorLength = interiorLength + 2 * ghostWidth;
////		for (int d = 0; d < dimCount; d++) {
////			coords[d].resize(exteriorLength + 1);
////			coordnames[d] = std::string(1, 'x' + d);
////			dims[d] = exteriorLength + 1;
////			auto const dx = domain_.span(d) / interiorLength;
////			for (int n = 0; n <= exteriorLength; n++) {
////				coords[d][n] = (n - ghostWidth) * dx;
////			}
////		}
////		for (int d = dimCount; d < siloDimCount; d++) {
////			coords[d].resize(2);
////			coordnames[d] = std::string(1, '_' + d);
////			dims[d] = 2;
////			coords[d][0] = Real(0.0);
////			coords[d][1] = Real(1.0);
////		}
////		for (int d = 0; d < siloDimCount; d++) {
////			coordPointers[d] = coords[d].data();
////			namePointers[d] = coordnames[d].c_str();
////		}
////		int const dataType = siloDatatype();
////		DBPutQuadmesh(db_, "mesh", namePointers.data(), coordPointers.data(), dims.data(), siloDimCount, dataType, DB_COLLINEAR, NULL);
////		hasMesh_ = true;
////	}
////	~Silo() {
////		DBClose(db_);
////	}
////
////private:
////	static constexpr int siloDatatype() {
////		if constexpr (std::is_same<Real, double>::value) {
////			return DB_DOUBLE;
////		} else if constexpr (std::is_same<Real, float>::value) {
////			return DB_FLOAT;
////		} else {
////			assert(false);
////		}
////	}
////	Interval<Real, dimCount> domain_{};
////	bool hasMesh_{};
////	int interiorLength_;
////	int ghostWidth_;
////	DBfile *db_{};
////};
////
////void output(std::string name, int number) {
////	name += "." + std::to_string(number) + ".silo";
////	Silo silo(name);
////	silo.writeMesh(interiorLength, ghostWidth);
////	silo.writeVariable("rho", &(Usoa.ρ[0]));
////	silo.writeVariable("tau", &(Usoa.τ[0]));
////	silo.writeVariable("E", &(Usoa.E[0]));
////	for (int d = 0; d < dimCount; d++) {
////		std::string const nm = "S" + std::string(1, 'x' + d);
////		silo.writeVariable(nm.c_str(), &(Usoa.s[d][0]));
////	}
////}
////
////void initialize() {
////	auto Uaos = soa2aos(Usoa);
////	for (unsigned i = ghostWidth; i < exteriorLength - ghostWidth; i++) {
////		bool const f = i < exteriorLength / 2;
////		auto const ρ = f ? 1.0 : 0.125;
////		auto const p = f ? 1.0 : 0.1;
////		Uaos[i].ρ = ρ;
////		Uaos[i].s = 0;
////		Uaos[i].E = p / (Gas::Γ - 1);
////		Uaos[i].τ = pow(Uaos[i].E, inv(Gas::Γ));
////	}
////	Usoa = aos2soa(Uaos);
////}
////
////void test() {
////	initialize();
////	Real timeNow = zero<Real>;
////	int stepNumber = 0;
////	while (timeNow < timeEnd) {
////		output("X", stepNumber);
////		enforceBoundaries();
////		computeFluxes();
////		Real const dt = std::min(timeEnd - timeNow, computeTimestep());
////		applyUpdate(dt);
////		printf("i = %i dt = %e\n", stepNumber, dt);
////		timeNow += dt;
////		stepNumber++;
////	}
////	output("X", stepNumber);
////}
////
int hpx_main(int argc, char *argv[]) {
	//	test();
	using T = double;
	constexpr size_t W = 56;
	SimdVector<T, W> v1, v2, v3;
	v1[v2 < v3] = 1.0;
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	hpx::init_params init_params;
	std::string sizeString = "hpx.stacks.small_size=" + std::to_string(size_t(1) << size_t(22)) + "\n";
	init_params.cfg.push_back("hpx.commandline.allow_unknown=1");
	init_params.cfg.push_back(sizeString);
	return hpx::init(argc, argv, init_params);
}
