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
#include "Silo.hpp"
#include "State.hpp"

#include <functional>

template <Integer order, Integer dimCount, Integer intWidth>
struct Grid : Allocator {
	static constexpr auto faceCount = 2_I * dimCount;
	static constexpr auto ghostWidth = 2_I;
	static constexpr auto extWidth = intWidth + 2_I * ghostWidth;
	static constexpr auto intBox = Box<dimCount>(intWidth);
	static constexpr auto extBox = intBox.expand(ghostWidth);
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
		std::array<Vector<dimCount>, extVolume> x;
		forEach(extBox, [&](auto const &idx) {
			auto const i = extBox.flatten(idx);
			for (Integer d = 0; d < dimCount; d++) {
				x[i][d] = cellWidth * (0.5_R + Real(idx[d]));
			}
		});
		return x;
	}();
	void output(Silo<dimCount> &silo) const {
		Vector<dimCount> const origin = -ghostWidth * cellWidth;
		silo.writeCoordinates(origin, cellWidth, extWidth);
		silo.writeData(un_.begin(), State<dimCount>::getFieldNames());
	}
	void reconstruct() {
		constexpr auto θ = 1.3_R;
		for (Integer d = 0; d < dimCount; d++) {
			auto const box = intBox.expand(1);
			forEach(box, [&](auto const &idx) {
				auto const i = extBox.flatten(idx);
				auto &rf = uf_[toRightFace(d)][i];
				auto &lf = uf_[toLeftFace(d)][i];
				if constexpr (order == 1) {
					rf = lf = un_[i];
				} else if constexpr (order == 2) {
					auto const di = gridStride[d];
					auto const up = un_[i + di];
					auto const u0 = un_[i];
					auto const um = un_[i - di];
					auto const ux = minmod(up - u0, u0 - um, θ);
					rf = un_[i] + 0.5_R * ux;
					lf = un_[i] - 0.5_R * ux;
				} else {
					static_assert(false);
				}
			});
		}
	}
	Real fluxes() {
		Real λ_max = 0_R;
		for (Integer d = 0; d < dimCount; d++) {
			auto const di = gridStride[d];
			forEach(fluxBox[d], [&](auto const &idx) {
				auto const i = extBox.flatten(idx);
				auto const &ul = uf_[toRightFace(d)][i - di];
				auto const &ur = uf_[toLeftFace(d)][i];
				auto const [f, λ_i] = riemannFlux(ul, ur, d);
				f_[d][i] = f;
				λ_max = std::max(λ_max, λ_i);
			});
		}
		return λ_max;
	}
	void update(Real dt, Real β) {
		constexpr auto dx = cellWidth;
		auto const λ = dt * inv(dx);
		forEach(intBox, [&](auto const &idx) {
			auto const i = extBox.flatten(idx);
			for (Integer d = 0; d < dimCount; d++) {
				auto const &fp = f_[d][i + gridStride[d]];
				auto const &fm = f_[d][i];
				un_[i] -= λ * (fp - fm);
			}
			un_[i] = β * un_[i] + (1_R - β) * u0_[i];
		});
	}
	void boundaries() {
		for (Integer d = 0; d < dimCount; d++) {
			auto const extL = bndBox[toLeftFace(d)];
			auto const extR = bndBox[toRightFace(d)];
			auto const intR = extL.shift(d, +intWidth);
			auto const intL = extR.shift(d, -intWidth);
			forEach(extL, intR, [&](auto const &extIdx, auto const &intIdx) {
				un_[extBox.flatten(extIdx)] = un_[extBox.flatten(intIdx)];
			});
			forEach(extR, intL, [&](auto const &extIdx, auto const &intIdx) {
				un_[extBox.flatten(extIdx)] = un_[extBox.flatten(intIdx)];
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
			un_[i] = foo(cellCenters[i]);
		});
		boundaries();
		store();
	}
	Grid() :
		Allocator(2 * sizeof(StateArray) + sizeof(FluxArray) + sizeof(ReconArray)),
		f_(get<FluxArray>()),
		uf_(get<ReconArray>()),
		un_(get<StateArray>()),
		u0_(get<StateArray>()) {
	}
	virtual ~Grid() {
	}

private:
	using StateArray = std::array<State<dimCount>, extVolume>;
	using ReconArray = std::array<StateArray, faceCount>;
	using FluxArray = std::array<StateArray, dimCount>;
	std::array<StateArray, dimCount> &f_;
	std::array<StateArray, faceCount> &uf_;
	StateArray &un_;
	StateArray &u0_;
};

#endif /* INCLUDE_GRID_HPP_ */
