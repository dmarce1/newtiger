/*
 * Silo.hpp
 *
 *  Created on: Feb 21, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_SILO_HPP_
#define INCLUDE_SILO_HPP_

#include <array>
#include <string>
#include <vector>

#include <silo.h>

#include "Integer.hpp"
#include "Real.hpp"
#include "Vector.hpp"

template <Integer D>
struct Silo {
	constexpr static int mode = DB_CLOBBER;
	constexpr static int target = DB_LOCAL;
	constexpr static int filetype = DB_HDF5;
	constexpr static int coordtype = DB_COLLINEAR;
	constexpr static int centering = (D == 1) ? DB_NODECENT : DB_ZONECENT;
	constexpr static int datatype = std::is_same_v<Real, double> ? DB_DOUBLE : DB_FLOAT;
	constexpr static int majororder = DB_COLMAJOR;
	constexpr static int maxopts = 8;
	constexpr static char const *meshname = "mesh";
	constexpr static char const *fileinfo = NULL;
	constexpr static char const *const coordnames[] = {"x", "y", "z", "w"};
	Silo(std::string const &fname) {
		char const *pathname = fname.c_str();
		db_ = DBCreate(pathname, mode, target, fileinfo, filetype);
		optlist_ = DBMakeOptlist(maxopts);
		int opt = majororder;
		DBAddOption(optlist_, DBOPT_MAJORORDER, &opt);
	}
	void writeCoordinates(auto const &origin, Real dx, Integer gridWidth) {
		std::array<std::vector<Real>, D> x;
		nels_ = 1;
		Integer const meshWidth = gridWidth + ((D == 1) ? 0 : 1);
		for (Integer d = 0; d < D; d++) {
			x[d].resize(meshWidth);
			for (Integer n = 0; n <= gridWidth; n++) {
				x[d][n] = origin[d] + Real(n) * dx;
			}
		}
		std::reverse(x.begin(), x.end());
		for (Integer d = 0; d < D; d++) {
			coords_[d] = static_cast<void const *>(x[d].data());
			meshDims_[d] = meshWidth;
			varDims_[d] = gridWidth;
			nels_ *= gridWidth;
		}
		DBPutQuadmesh(db_, meshname, coordnames, coords_.data(), meshDims_.data(), D, datatype, coordtype, optlist_);
	}
	void writeData(auto begin, auto const &names) {
		std::vector<Real> var(nels_);
		void const *varptr = static_cast<void const *>(var.data());
		Integer const nfields = names.size();
		for (Integer f = 0; f < nfields; f++) {
			for (auto n = 0, it = begin; n < nels_; ++n, ++it) {
				var[n] = (*it)[f];
			}
			DBPutQuadvar1(db_, names[f].c_str(), meshname, varptr, varDims_.data(), D, NULL, 0, datatype, centering, optlist_);
		}
	}
	~Silo() {
		DBFreeOptlist(optlist_);
		DBClose(db_);
	}

private:
	DBfile *db_;
	DBoptlist *optlist_;
	int nels_;
	std::array<void const *, D> coords_;
	std::array<int, D> meshDims_;
	std::array<int, D> varDims_;
};

#endif /* INCLUDE_SILO_HPP_ */
