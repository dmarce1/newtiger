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
#include "GasState.hpp"
#include "Limiters.hpp"
#include "RadiationState.hpp"
#include "Silo.hpp"
#include "Simd.hpp"

#include <functional>
#include <iostream>

template <Integer order, Integer dimCount, Integer intWidth, template <typename, Integer> typename State>
struct Grid {
	static constexpr auto fieldCount = State<Real, dimCount>::size();
	static constexpr auto faceCount = 2_I * dimCount;
	static constexpr auto ghostWidth = 2_I;
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
				x[i][d] = cellWidth * (0.5_R + Real(idx[d]) - Real(intWidth) * 0.5_R);
			}
		});
		return x;
	}();
	void output(Silo<dimCount> &silo) const {
		Vector<Real, dimCount> const origin = -ghostWidth * cellWidth;
		silo.writeCoordinates(origin, cellWidth, extWidth);
		auto const &names = State<Real, dimCount>::getFieldNames();
		for (Integer f = 0; f < fieldCount; f++) {
			silo.writeData(u_[f].begin(), names[f]);
		}
	}
	void reconstruct() {
		if constexpr (order == 1) return;
		constexpr auto θ = (order <= 2) ? 1.333_R : 2_R;
		for (Integer d = 0; d < dimCount; d++) {
			forEachSimd<maxSimdSize<Real>()>(intBox.pad(d, std::pair(1, 1)), [&]<Integer W>(auto const &idx) {
				using SimdType = Simd<Real, W>;
				State<SimdType, dimCount> u0, ur, ul;
				Integer const i = extBox.flatten(idx);
				for (Integer f = 0; f < fieldCount; f++) {
					u0[f].load(&u_[f][i]);
				}
				State<SimdType, dimCount> up, um, Δ, Δc;
				auto const di = gridStride[d];
				for (Integer f = 0; f < fieldCount; f++) {
					up[f].load(&u_[f][i + di]);
					um[f].load(&u_[f][i - di]);
					if constexpr (order >= 2) {
						Δc[f].load(&du_[f][d][i - di]);
					}
				}
				auto Δp = 0.5_R * (up - u0);
				auto Δm = 0.5_R * (u0 - um);
				State<SimdType, dimCount> state(u0);
				auto const [_, R, L] = state.eigenstructure(d);
				Δp = L(Δp);
				Δm = L(Δm);
				if constexpr (order < 2) {
					Δc = 0.5_R * (Δp + Δm);
				}
				for (Integer f = 0; f < fieldCount; f++) {
					Δ[f] = minmod(Δc[f], θ * minmod(Δp[f], Δm[f]));
				}
				Δ = R(Δ);
				for (Integer f = 0; f < fieldCount; f++) {
					Δ[f].store(&du_[f][d][i]);
				}
			});
		}
	}
	Real fluxes() {
		Real λ_max = 0_R;
		for (Integer d = 0; d < dimCount; d++) {
			using std::max;
			auto const di = gridStride[d];
			auto const thisBox = intBox.pad(d, std::pair(0_I, 1_I));
			forEachSimd<maxSimdSize<Real>()>(thisBox, [&]<Integer W>(auto const &idx) {
				using SimdType = Simd<Real, W>;
				Integer const i = extBox.flatten(idx);
				State<SimdType, dimCount> ur, ul;
				if constexpr (order == 1) {
					for (Integer f = 0; f < fieldCount; f++) {
						ul[f].load(&u_[f][i - di]);
						ur[f].load(&u_[f][i]);
					}
					auto const [flux, λ_i] = riemannFlux(ul, ur, d);
					for (Integer f = 0; f < fieldCount; f++) {
						flux[f].store(&f_[f][d][i]);
					}
					λ_max = max(λ_max, max(λ_i));
				} else if constexpr (order == 2) {
					State<SimdType, dimCount> ur0, ur1, ul0, ul1;
					for (Integer f = 0; f < fieldCount; f++) {
						ul0[f].load(&u_[f][i - di]);
						ul1[f].load(&du_[f][d][i - di]);
						ur0[f].load(&u_[f][i]);
						ur1[f].load(&du_[f][d][i]);
					}
					ul = ul0 + ul1;
					ur = ur0 - ur1;
					auto const [flux, λ_i] = riemannFlux(ul, ur, d);
					for (Integer f = 0; f < fieldCount; f++) {
						flux[f].store(&f_[f][d][i]);
					}
					λ_max = max(λ_max, max(λ_i));
				} else if constexpr (order == 3) {
					constexpr Integer quadCount = 1_I << dimCount;
					constexpr Real dx = inv(sqrt(3_R));
					constexpr Real wt0 = inv(quadCount / 2);
					constexpr Real wt1 = sqrt(3_R) * wt0;
					State<SimdType, dimCount> ur0, ul0;
					State<SimdType, dimCount> flux0{};
					std::array<State<SimdType, dimCount>, dimCount> flux1{};
					std::array<State<SimdType, dimCount>, dimCount> ur1, ul1;
					for (Integer f = 0; f < fieldCount; f++) {
						ul0[f].load(&u_[f][i - di]);
						ur0[f].load(&u_[f][i]);
						for (Integer d1 = 0; d1 < dimCount; d1++) {
							ul1[d1][f].load(&du_[f][d1][i - di]);
							ur1[d1][f].load(&du_[f][d1][i]);
						}
					}
					for (Integer j = 0; j < quadCount; j++) {
						if ((j >> d) & 1_I) continue;
						std::array<Real, dimCount> sgn;
						for (Integer d1 = 0; d1 < dimCount; d1++) {
							sgn[d1] = Real(2_I * ((j >> d1) & 1_I) - 1_I);
						}
						ul = ul0 + ul1[d];
						ur = ur0 - ur1[d];
						for (Integer d1 = 0; d1 < dimCount; d1++) {
							if (d1 == d) continue;
							ul += ul1[d1] * sgn[d1] * dx;
							ur += ur1[d1] * sgn[d1] * dx;
						}
						auto const [thisFlux, λ_i] = riemannFlux(ul, ur, d);
						for (Integer f = 0; f < fieldCount; f++) {
							flux0[f] += thisFlux[f] * wt0;
							for (Integer d1 = 0; d1 < dimCount; d1++) {
								if (d1 == d) continue;
								flux1[d1][f] += thisFlux[f] * sgn[d1] * wt1;
							}
						}
						λ_max = max(λ_max, max(λ_i));
					}
					for (Integer f = 0; f < fieldCount; f++) {
						flux0[f].store(&f_[f][d][i]);
						flux1[d][f].store(&df_[d][f][d][i]);
						for (Integer d1 = 0; d1 < dimCount; d1++) {
							if (d1 == d) continue;
							flux1[d1][f].store(&df_[d1][f][d][i]);
						}
					}
				} else {
					static_assert(false);
				}
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
			if constexpr (order < 3) {
				State<SimdType, dimCount> u0, un;
				Vector<State<SimdType, dimCount>, dimCount> fp, fm;
				for (Integer f = 0; f < fieldCount; f++) {
					un[f].load(&u_[f][i]);
					u0[f].load(&u0_[f][i]);
					for (Integer d = 0; d < dimCount; d++) {
						auto const di = gridStride[d];
						fp[d][f].load(&f_[f][d][i + di]);
						fm[d][f].load(&f_[f][d][i]);
					}
				}
				for (Integer d = 0; d < dimCount; d++) {
					un -= λ * (fp[d] - fm[d]);
				}
				un = β * un + (1_R - β) * u0;
				un = un.updateEntropy();
				for (Integer f = 0; f < fieldCount; f++) {
					un[f].store(&u_[f][i]);
				}
			} else if constexpr (order == 3) {
				constexpr Integer quadCount = 1_I << dimCount;
				constexpr Real dx = sqrt(3_R) / 3_R;
				constexpr Real wt0 = inv(quadCount);
				constexpr Real wt1 = 3_R * wt0;
				State<SimdType, dimCount> u0, un;
				Vector<State<SimdType, dimCount>, dimCount> du0, dun;
				Vector<State<SimdType, dimCount>, dimCount> fp, fm;
				Vector<Vector<State<SimdType, dimCount>, dimCount>, dimCount> dfp, dfm;
				for (Integer f = 0; f < fieldCount; f++) {
					u0[f].load(&u0_[f][i]);
					un[f].load(&u_[f][i]);
					for (Integer d1 = 0; d1 < dimCount; d1++) {
						auto const di = gridStride[d1];
						fp[d1][f].load(&f_[f][d1][i + di]);
						fm[d1][f].load(&f_[f][d1][i]);
						dun[d1][f].load(&du_[f][d1][i]);
						du0[d1][f].load(&du0_[f][d1][i]);
						for (Integer d2 = 0; d2 < dimCount; d2++) {
							dfp[d2][d1][f].load(&df_[f][d2][d1][i + di]);
							dfm[d2][d1][f].load(&df_[f][d2][d1][i]);
						}
					}
				}
				auto const dun_ = dun;
				for (Integer j = 0; j < quadCount; j++) {
					State<SimdType, dimCount> u = un;
					for (Integer d1 = 0; d1 < dimCount; d1++) {
						auto const sgn = Real(2_I * ((j >> d1) & 1_I) - 1_I);
						u += sgn * dx * dun_[d1];
					}
					for (Integer d2 = 0; d2 < dimCount; d2++) {
						State<SimdType, dimCount> const f = u.flux(d2);
						dun[d2] += λ * wt1 * f;
					}
				}
				for (Integer d1 = 0; d1 < dimCount; d1++) {
					un -= λ * (fp[d1] - fm[d1]);
					for (Integer d2 = 0; d2 < dimCount; d2++) {
						Real const sgn = (d1 == d2) ? +1_R : -1_R;
						dun[d2] -= λ * (dfp[d2][d1] + sgn * dfm[d2][d1]);
					}
				}
				un = β * un + (1_R - β) * u0;
				un = un.updateEntropy();
				for (Integer f = 0; f < fieldCount; f++) {
					un[f].store(&u_[f][i]);
					for (Integer d2 = 0; d2 < dimCount; d2++) {
						dun[d2][f].store(&du_[f][d2][i]);
					}
				}
			} else {
				static_assert(false);
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
					u_[f][j] = u_[f][k];
					if constexpr (order > 2) {
						for (Integer d = 0; d < dimCount; d++) {
							du_[f][d][j] = du_[f][d][k];
						}
					}
				}
			});
			forEach(extR, intL, [&](auto const &extIdx, auto const &intIdx) {
				Integer const j = extBox.flatten(extIdx);
				Integer const k = extBox.flatten(intIdx);
				for (Integer f = 0; f < fieldCount; f++) {
					u_[f][j] = u_[f][k];
					if constexpr (order > 2) {
						for (Integer d = 0; d < dimCount; d++) {
							du_[f][d][j] = du_[f][d][k];
						}
					}
				}
			});
		}
	}
	void store() {
		u0_ = u_;
		if constexpr (order > 2) {
			du0_ = du_;
		}
	}
	template <typename F>
	void initialize(F const &foo) {
		if constexpr (order > 2) {
			constexpr Integer quadCount = 1_I << dimCount;
			constexpr Real dx = cellWidth * sqrt(3_R) / 6_R;
			constexpr Real wt0 = inv(quadCount);
			constexpr Real wt1 = sqrt(3_R) * wt0;
			forEach(intBox, [&](auto const &idx) {
				Integer const i = extBox.flatten(idx);
				State<Real, dimCount> u{};
				std::array<Real, dimCount> sgn;
				Vector<State<Real, dimCount>, dimCount> du{};
				for (Integer j = 0; j < quadCount; j++) {
					auto x = cellCenters[i];
					for (Integer d = 0; d < dimCount; d++) {
						sgn[d] = Real(2_I * ((j >> d) & 1_I) - 1_I);
						x[d] += sgn[d] * dx;
					}
					auto const v = foo(x);
					u += wt0 * v;
					for (Integer d = 0; d < dimCount; d++) {
						du[d] += sgn[d] * wt1 * v;
					}
				}
				for (Integer f = 0; f < fieldCount; f++) {
					u_[f][i] = u[f];
					for (Integer d = 0; d < dimCount; d++) {
						du_[f][d][i] = du[d][f];
					}
				}
			});
		} else {
			forEach(intBox, [&](auto const &idx) {
				Integer const i = extBox.flatten(idx);
				auto const u = foo(cellCenters[i]);
				for (Integer f = 0; f < fieldCount; f++) {
					u_[f][i] = u[f];
				}
			});
		}
		boundaries();
		store();
	}
	Grid() {
	}
	virtual ~Grid() {
	}

private:
	using StateArray0 = std::vector<std::array<Real, extVolume>>;
	using StateArray1 = std::vector<std::array<std::array<Real, extVolume>, dimCount>>;
	using StateArray2 = std::vector<std::array<std::array<Real, extVolume>, dimCount *(dimCount + 1) / 2>>;
	using FluxArray0 = std::vector<std::array<std::array<Real, extVolume>, dimCount>>;
	using FluxArray1 = std::vector<std::array<std::array<std::array<Real, extVolume>, dimCount>, dimCount>>;
	FluxArray0 f_{fieldCount};
	StateArray0 u_{fieldCount};
	StateArray0 u0_{fieldCount};
	StateArray1 du_{(order > 1_I) ? fieldCount : 0_I};
	FluxArray1 df_{(order > 2_I) ? fieldCount : 0_I};
	StateArray1 du0_{(order > 2_I) ? fieldCount : 0_I};
	StateArray2 d2u_{(order > 2_I) ? fieldCount : 0_I};
	StateArray0 src_{(order > 2_I) ? fieldCount : 0_I};
};

#endif /* INCLUDE_GRID_HPP_ */
