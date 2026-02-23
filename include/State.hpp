/*
 * State.hpp
 *
 *  Created on: Feb 21, 2026
 *      Author: dmarce1
 */

#ifndef STATE_HPP_
#define STATE_HPP_

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#include <string>
#include <tuple>

#include "DoubleValue.hpp"
#include "Real.hpp"
#include "Vector.hpp"

enum class RiemannSolver : int { LF, HLL, HLLC, HLLE, HLLCE };

constexpr Integer HLLC = 0x1;
constexpr Integer HLLE = 0x2;

constexpr auto select(bool f, auto a, auto b) {
	return f ? a : b;
}

template <Integer D, typename T = Real, RiemannSolver R = RiemannSolver::HLLCE>
struct State : Vector<2 + D> {
	using base_type = Vector<2 + D>;
	static constexpr Integer size() {
		return D + 2;
	}
	static constexpr Real Γ = 5_R / 3_R;
	static auto getFieldNames() {
		static auto names = []() {
			static std::array<std::string, 2 + D> strings;
			strings[0] = "density";
			for (Integer d = 0; d < D; d++) {
				strings[d + 1] = std::string(1, 'x' + d) + "_momentum";
			}
			strings[1 + D] = "energy";
			return strings;
		}();
		return names;
	}
	State() :
		density(*(new(&(*this)[0]) T{})), momentum(*(new(&(*this)[1]) Vector<D, T>{})), energy(*(new(&(*this)[1 + D]) T{})) {
	}
	State(State const &other) :
		density(*(new(&(*this)[0]) T{})), momentum(*(new(&(*this)[1]) Vector<D, T>{})), energy(*(new(&(*this)[1 + D]) T{})) {
		*this = other;
	}
	State(base_type const &other) :
		density(*(new(&(*this)[0]) T{})), momentum(*(new(&(*this)[1]) Vector<D, T>{})), energy(*(new(&(*this)[1 + D]) T{})) {
		*this = other;
	}
	virtual ~State() {
		density.~T();
		energy.~T();
		momentum.~Vector<D, T>();
	}
	State &operator=(State const &other) {
		base_type::operator=(other);
		return *this;
	}
	State &operator=(base_type const &other) {
		base_type::operator=(other);
		return *this;
	}
	auto flux(Real u, Real p, Integer n) const {
		State f = *this * u;
		f.momentum[n] += p;
		f.energy += u * p;
		return f;
	};
	friend auto riemannFlux(State const &UL, State const &UR, Integer n) {
		using std::abs;
		using std::max;
		using std::min;
		using std::sqrt;
		constexpr auto zero = 0.0_R;
		constexpr auto half = 0.5_R;
		constexpr auto one = 1.0_R;
		constexpr auto Γm1 = (Γ - one);
		auto const ρR = UR.density;
		auto const iρR = inv(ρR);
		auto const eR = UR.energy;
		auto const uR = iρR * UR.momentum;
		auto const ekR = half * ρR * uR.dot(uR);
		auto const eiR = max(eR - ekR, zero);
		auto const pR = Γm1 * eiR;
		auto const cR = sqrt(Γ * pR * iρR);
		auto const ρL = UL.density;
		auto const iρL = inv(ρL);
		auto const eL = UL.energy;
		auto const uL = iρL * UL.momentum;
		auto const ekL = half * ρL * uL.dot(uL);
		auto const eiL = max(eL - ekL, zero);
		auto const pL = Γm1 * eiL;
		auto const cL = sqrt(Γ * pL * iρL);
		if constexpr (R == RiemannSolver::LF) {
			auto const FR = UR.flux(uR[n], pR, n);
			auto const FL = UL.flux(uL[n], pL, n);
			auto const sR = abs(uR[n]) + cR;
			auto const sL = abs(uL[n]) + cL;
			auto const s = max(sR, sL);
			auto const F = half * (FR + FL - s * (UR - UL));
			return std::pair(F, s);
		}
		auto sR = max(uR[n] + cR, zero);
		auto sL = min(uL[n] - cL, zero);
		if constexpr ((R == RiemannSolver::HLLE) || (R == RiemannSolver::HLLCE)) {
			auto const wR = one / (one + sqrt(ρL) / sqrt(ρR));
			auto const wL = one - wR;
			auto const hR = iρR * (eR + pR);
			auto const hL = iρL * (eL + pL);
			auto const u0 = wR * uR + wL * uL;
			auto const h0 = wR * hR + wL * hL;
			auto const c0 = sqrt(h0 - half * u0.dot(u0));
			sR = max(sR, u0[n] + c0);
			sL = min(sL, u0[n] - c0);
		} else {
			sR = max(sR, uL[n] + cL);
			sL = min(sL, uR[n] - cR);
		}
		State F;
		if constexpr ((R == RiemannSolver::HLL) || (R == RiemannSolver::HLLE)) {
			auto const FR = UR.flux(uR[n], pR, n);
			auto const FL = UL.flux(uL[n], pL, n);
			F = (sR * FL - sL * FR + sL * sR * (UR - UL)) / (sR - sL);
		} else {
			auto const wL = ρL * (sL - uL[n]);
			auto const wR = ρR * (sR - uR[n]);
			auto const s0 = ((pR - pL) + (wL * uL[n] - wR * uR[n])) / (wL - wR);
			bool const f = s0 > zero;
			auto const ρK = select(f, ρL, ρR);
			auto const pK = select(f, pL, pR);
			auto const uK = select(f, uL, uR);
			auto const sK = select(f, sL, sR);
			auto const eK = select(f, eL, eR);
			auto const wK = select(f, wL, wR);
			auto const w0 = wK / (sK - s0);
			auto const e0 = w0 * (eK / ρK + (s0 - uK[n]) * (s0 + pK / wK));
			auto u0 = uK;
			u0[n] = s0;
			State U0;
			U0.density = w0;
			U0.momentum = w0 * u0;
			U0.energy = e0;
			auto const FK = UK.flux(uK[n], pK, n);
			auto const UK = select(f, UL, UR);
			F = FK + sK * (U0 - UK);
		}
		auto const s = max(sR, -sL);
		return std::pair(F, s);
	}
	void setDensity(Real d) {
		density = d;
	}
	void setPressure(Real p) {
		energy = p * inv(Γ - 1_R) + 0.5_R * inv(density) * momentum.dot(momentum);
	}

private:
	Real &density;
	Vector<D> &momentum;
	Real &energy;
};

#endif /* STATE_HPP_ */
