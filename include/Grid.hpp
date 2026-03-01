/*
 * Grid.hpp
 *
 *  Created on: Feb 21, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_GRID_HPP_
#define INCLUDE_GRID_HPP_

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#include "Allocator.hpp"
#include "Box.hpp"
#include "GasState.hpp"
#include "Limiters.hpp"
#include "RadiationState.hpp"
#include "Silo.hpp"
#include "Simd.hpp"

#include <functional>
#include <iostream>

template <Integer order, Integer dimCount, Integer intWidth>
struct Grid : Allocator {
	static constexpr auto simdWidth = 4;
	static constexpr auto fieldCount = GasState<Real, dimCount>::size();
	static constexpr auto faceCount = 2_I * dimCount;
	static constexpr auto ghostWidth = 4_I;
	static constexpr auto extWidth = intWidth + 2_I * ghostWidth;
	static constexpr auto intBox = Box<dimCount>(intWidth);
	static constexpr auto extBox = intBox.pad(ghostWidth);
	static constexpr auto fluxBox = []() {
		std::array<Box<dimCount>, dimCount> fbox;
		auto const lb = intBox.begin();
		auto ub = intBox.end();
		for (Integer d = 0; d < dimCount; d++) {
			ub[d]++;
			fbox[d] = Box<dimCount>(lb, ub);
			ub[d]--;
		}
		return fbox;
	}();
	static constexpr Integer toRightFace(Integer dim) {
		return 2 * dim + 1;
	}
	static constexpr Integer toLeftFace(Integer dim) {
		return 2 * dim;
	}
	static constexpr Integer faceDim(Integer f) {
		return f >> 1_I;
	}
	static constexpr Real faceSign(Integer f) {
		return Real(2_I * (f & 1_I) - 1_I);
	}
	static constexpr auto bndBox = []() {
		std::array<Box<dimCount>, faceCount> bbox;
		std::array<Integer, dimCount> lbl, ubl, lbr, ubr;
		for (Integer d = 0; d < dimCount; d++) {
			lbr = lbl = extBox.begin();
			ubr = ubl = extBox.end();
			lbl[d] = -ghostWidth;
			ubl[d] = 0_I;
			lbr[d] = intWidth;
			ubr[d] = lbr[d] + ghostWidth;
			bbox[toLeftFace(d)] = Box<dimCount>(lbl, ubl);
			bbox[toRightFace(d)] = Box<dimCount>(lbr, ubr);
		}
		return bbox;
	}();
	static constexpr auto intVolume = intBox.volume();
	static constexpr auto extVolume = extBox.volume();
	static constexpr auto cellWidth = inv(intWidth);
	static constexpr auto gridStride = []() {
		std::array<Integer, dimCount> str;
		str.back() = 1;
		for (Integer d = dimCount - 1; d > 0; d--) {
			str[d - 1] = str[d] * extWidth;
		}
		return str;
	}();
	static constexpr auto cellCenters = []() {
		std::array<Vector<Real, dimCount>, extVolume> x;
		forEach(extBox, [&](auto const &idx) {
			auto const i = extBox.flatten(idx);
			for (Integer d = 0; d < dimCount; d++) {
				x[i][d] = cellWidth * (0.5_R + Real(idx[d]));
			}
		});
		return x;
	}();
	void output(Silo<dimCount> &silo) const {
		Vector<Real, dimCount> const origin = -ghostWidth * cellWidth;
		silo.writeCoordinates(origin, cellWidth, extWidth);
		auto const &names = GasState<Real, dimCount>::getFieldNames();
		for (Integer f = 0; f < fieldCount; f++) {
			silo.writeData(un_[f].begin(), names[f]);
		}
	}
	void reconstruct() {
		constexpr auto θ = 1.3_R;
		for (Integer d = 0; d < dimCount; d++) {
			forEachSimd<maxSimdSize<Real>()>(intBox.pad(d, std::pair(1, 1)), [&]<Integer W>(auto const &idx) {
				using SimdType = Simd<Real, W>;
				GasState<SimdType, dimCount> u0, ur, ul;
				Integer const i = extBox.flatten(idx);
				for (Integer f = 0; f < fieldCount; f++) {
					u0[f].load(&un_[f][i]);
				}
				if constexpr (order == 1) {
					ur = ul = u0;
				} else if constexpr (order == 2) {
					GasState<SimdType, dimCount> up, um, Δ;
					auto const di = gridStride[d];
					for (Integer f = 0; f < fieldCount; f++) {
						up[f].load(&un_[f][i + di]);
						um[f].load(&un_[f][i - di]);
					}
					auto Δp = up - u0;
					auto Δm = u0 - um;
					GasState<SimdType, dimCount> GasState(u0);
					auto const [_, R, L] = GasState.eigenstructure(d);
					Δp = L * Δp;
					Δm = L * Δm;
					auto const Δc = 0.5_R * (Δp + Δm);
					for (Integer f = 0; f < fieldCount; f++) {
						Δ[f] = minmod(Δc[f], θ * minmod(Δp[f], Δm[f]));
					}
					Δ = R * Δ;
					ur = u0 + 0.5_R * Δ;
					ul = u0 - 0.5_R * Δ;
				} else {
					static_assert(false);
				}
				for (Integer f = 0; f < fieldCount; f++) {
					ur[f].store(&uf_[toRightFace(d)][f][i]);
					ul[f].store(&uf_[toLeftFace(d)][f][i]);
				}
			});
		}
	}
	Real fluxes() {
		Real λ_max = 0_R;
		for (Integer d = 0; d < dimCount; d++) {
			using std::max;
			auto const di = gridStride[d];
			forEachSimd<maxSimdSize<Real>()>(intBox.pad(d, std::pair(0, 1)), [&]<Integer W>(auto const &idx) {
				using SimdType = Simd<Real, W>;
				Integer const i = extBox.flatten(idx);
				GasState<SimdType, dimCount> ur, ul, Δ;
				for (Integer f = 0; f < fieldCount; f++) {
					ur[f].load(&uf_[toLeftFace(d)][f][i]);
					ul[f].load(&uf_[toRightFace(d)][f][i - di]);
				}
				auto const [flux, λ_i] = riemannFlux(ul, ur, d);
				for (Integer f = 0; f < fieldCount; f++) {
					flux[f].store(&f_[d][f][i]);
				}
				λ_max = max(λ_max, max(λ_i));
			});
		}
		return λ_max;
	}
	void update(Real dt, Real β) {
		constexpr auto dx = cellWidth;
		auto const λ = dt * inv(dx);
		forEachSimd<maxSimdSize<Real>()>(intBox, [&]<Integer W>(auto const &idx) {
			using SimdType = Simd<Real, W>;
			auto const i = extBox.flatten(idx);
			GasState<SimdType, dimCount> u0, un;
			Vector<GasState<SimdType, dimCount>, dimCount> fp, fm;
			for (Integer f = 0; f < fieldCount; f++) {
				un[f].load(&un_[f][i]);
				u0[f].load(&u0_[f][i]);
				for (Integer d = 0; d < dimCount; d++) {
					auto const di = gridStride[d];
					fp[d][f].load(&f_[d][f][i + di]);
					fm[d][f].load(&f_[d][f][i]);
				}
			}
			for (Integer d = 0; d < dimCount; d++) {
				un -= λ * (fp[d] - fm[d]);
			}
			un = β * un + (1_R - β) * u0;
			un = un.updateInternalEnergy();
			for (Integer f = 0; f < fieldCount; f++) {
				un[f].store(&un_[f][i]);
			}
		});
	}
	void boundaries() {
		for (Integer d = 0; d < dimCount; d++) {
			auto const extL = bndBox[toLeftFace(d)];
			auto const extR = bndBox[toRightFace(d)];
			auto const intR = extL.shift(d, +intWidth);
			auto const intL = extR.shift(d, -intWidth);
			forEach(extL, intR, [&](auto const &extIdx, auto const &intIdx) {
				Integer const j = extBox.flatten(extIdx);
				Integer const k = extBox.flatten(intIdx);
				for (Integer f = 0; f < fieldCount; f++) {
					un_[f][j] = un_[f][k];
				}
			});
			forEach(extR, intL, [&](auto const &extIdx, auto const &intIdx) {
				Integer const j = extBox.flatten(extIdx);
				Integer const k = extBox.flatten(intIdx);
				for (Integer f = 0; f < fieldCount; f++) {
					un_[f][j] = un_[f][k];
				}
			});
		}
	}
	void store() {
		u0_ = un_;
	}
	template <typename F>
	void initialize(F const &foo) {
		forEach(intBox, [&](auto const &idx) {
			Integer const i = extBox.flatten(idx);
			auto const u = foo(cellCenters[i]);
			for (Integer f = 0; f < fieldCount; f++) {
				un_[f][i] = u[f];
			}
		});
		boundaries();
		store();
	}
	Grid() :
		Allocator(2 * sizeof(GasStateArray) + sizeof(FluxArray) + sizeof(ReconArray)),
		f_(get<FluxArray>()),
		uf_(get<ReconArray>()),
		un_(get<GasStateArray>()),
		u0_(get<GasStateArray>()) {
	}
	virtual ~Grid() {
	}

private:
	using GasStateArray = std::array<std::array<Real, extVolume>, fieldCount>;
	using ReconArray = std::array<GasStateArray, faceCount>;
	using FluxArray = std::array<GasStateArray, dimCount>;
	FluxArray &f_;
	ReconArray &uf_;
	GasStateArray &un_;
	GasStateArray &u0_;
};

#endif /* INCLUDE_GRID_HPP_ */
