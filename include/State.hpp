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

#include <tuple>

#include "Real.hpp"
#include "Vector.hpp"

template <Integer D>
struct State : Vector<2 + D> {
	using base_type = Vector<2 + D>;
	static constexpr Real Γ = 5_R / 3_R;
	static constexpr Integer size() {
		return D + 2;
	}
	State() :
		ρ((*this)[0]), e((*this)[1]), s(reinterpret_cast<Vector<D> &>((*this)[2])) {
	}
	State(base_type const &other) :
		ρ((*this)[0]), e((*this)[1]), s(reinterpret_cast<Vector<D> &>((*this)[2])) {
		*this = other;
	}
	State &operator=(State const &other) {
		base_type::operator=(other);
		return *this;
	}
	auto primitives() const {
		auto const iρ = inv(ρ);
		auto const v = iρ * s;
		auto const ek = 0.5_R * v.dot(s);
		auto const ei = std::max(e - ek, 0_R);
		auto const p = (Γ - 1_R) * ei;
		auto const c = sqrt(Γ * p * iρ);
		return std::tuple(v, p, c);
	};
	auto flux(Integer n) const {
		auto const [v, p, c] = primitives();
		auto const a = std::abs(v[n]) + c;
		State f = *this * v[n];
		f.s[n] += p;
		f.e += v[n] * p;
		return std::pair<State, Real>(f, a);
	};
	friend auto riemannFlux(State const &ul, State const &ur, Integer n) {
		auto const [fl, al] = ul.flux(n);
		auto const [fr, ar] = ur.flux(n);
		auto const a = std::max(al, ar);
		auto const f = 0.5_R * ((fl + fr) - a * (ur - ul));
		return std::pair<State, Real>(f, a);
	}
	Real maxEigenvalue(Integer n) const {
		auto const [v, _, c] = primitives();
		return std::abs(v[n]) + c;
	}
	void setDensity(Real d) {
		ρ = d;
	}
	void setPressure(Real p) {
		e = p * inv(Γ - 1_R) + 0.5_R * inv(ρ) * s.dot(s);
	}

private:
	Real &ρ;
	Real &e;
	Vector<D> &s;
};

#endif /* STATE_HPP_ */
