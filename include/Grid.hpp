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

#include "Box.hpp"
#include "Silo.hpp"
#include "State.hpp"

template <Integer dimCount, Integer intWidth>
struct Grid {
	constexpr static auto faceCount = 2_I * dimCount;
	constexpr static auto ghostWidth = 2_I;
	constexpr static auto intBox = Box<dimCount>(intWidth);
	constexpr static auto intVolume = intBox.volume();
	constexpr static auto extWidth = intWidth + 2_I * ghostWidth;
	constexpr static auto extBox = Box<dimCount>(extWidth);
	constexpr static auto extVolume = extBox.volume();
	constexpr static auto cellWidth = Real(intWidth);
	constexpr static auto fluxBox = []() {
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
	constexpr static auto gridStride = []() {
		std::array<Integer, dimCount> str;
		str.back() = 1;
		for (Integer d = dimCount - 1; d > 0; d--) {
			str[d - 1] = str[d] * extWidth;
		}
		return str;
	}();
	constexpr static auto cellCenters = []() {
		std::array<Vector<dimCount>, extVolume> x;
		forEach(extBox, [&](auto const &idx) {
			auto const i = extBox.flatten(idx);
			for (Integer d = 0; d < dimCount; d++) {
				x[i][d] = cellWidth * (0.5_R + idx[d] - ghostWidth);
			}
		});
		return x;
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
	void output(Silo<dimCount> &silo) const {
		Vector<dimCount> const origin = -ghostWidth * cellWidth;
		silo.writeCoordinates(origin, cellWidth, extWidth);
		silo.writeData(un_.begin(), State<dimCount>::getFieldNames());
	}
	void reconstruct() {
		constexpr auto θ = 1.3_R;
		for (Integer d = 0; d < dimCount; d++) {
			auto const stride = gridStride[d];
			forEach(intBox, [&](auto const &idx) {
				auto const i = extBox.flatten(idx);
				auto const ux = minmod(un_[i + stride] - un_[i], un_[i] - un_[i - stride], θ);
				uf_[toRightFace(d)][i] = un_[i] + 0.5_R * ux;
				uf_[toLeftFace(d)][i] = un_[i] - 0.5_R * ux;
			});
		}
	}
	Real fluxes() {
		Real λ_max = 0_R;
		for (Integer d = 0; d < dimCount; d++) {
			forEach(fluxBox[d], [&](auto const &idx) {
				auto const i = extBox.flatten(idx);
				auto const &ul = uf_[toRightFace(d)][i];
				auto const &ur = uf_[toLeftFace(d)][i - gridStride[d]];
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
			for (Integer n = 0; n < ghostWidth; n++) {
				auto const rightBox = extBox.slice(d, intWidth + n);
				auto const leftBox = extBox.slice(d, n - ghostWidth);
				forEach(rightBox, [&](auto const &idx) {
					Integer const i = extBox.flatten(idx);
					un_[i] = un_[i + intWidth * gridStride[d]];
				});
				forEach(leftBox, [&](auto const &idx) {
					Integer const i = extBox.flatten(idx);
					un_[i] = un_[i - intWidth * gridStride[d]];
				});
			}
		}
	}
	void store() {
		u0_ = un_;
	}
	template <typename F>
	void initialize(F const &foo) {
		forEach(intBox, [&](auto const &idx) {
			Integer const i = extBox.flatten(idx);
			un_[i] = foo(idx);
		});
		boundaries();
		store();
	}

private:
	using StateArray = std::array<State<dimCount>, extVolume>;
	std::array<StateArray, faceCount> f_;
	std::array<StateArray, faceCount> uf_;
	StateArray un_;
	StateArray u0_;
};

#endif /* INCLUDE_GRID_HPP_ */
