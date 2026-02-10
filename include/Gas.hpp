/*
 * Gas.hpp
 *
 *  Created on: Feb 9, 2026
 *      Author: dmarce1
 */

#ifndef GAS_HPP_
#define GAS_HPP_

#include <cassert>
#include <cmath>
/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#define STATE_NAME Gas
	
#define STATE_DATA \
	SCALAR(ρ) \
	VECTOR(s) \
	SCALAR(E) \
	SCALAR(τ) 
	
#define STATE_CONSTANTS \
	CONSTANT(Γ, 1.4)
	
#include "StateFactory.hpp"

using namespace config;

using ConservedGasState  = Gas::State;

struct PrimitiveGasState {
	Real ρ;
	Vector<Real, dimCount> v;
	Real p;
};


PrimitiveGasState con2prim(ConservedGasState const &con) {
	using namespace Gas;
	using std::pow;
	PrimitiveGasState prim;
	constexpr Real zero = Real(0);
	constexpr Real one = Real(1);
	constexpr Real half = inv(Real(2));
	constexpr Real δ₁ = inv(Real(1000));
	auto const &[_, s, E, τ] = con;
	auto &[ρ, v, p] = prim;
	ρ = con.ρ;
	assert(ρ > zero);
	assert(τ > zero);
	auto const iρ = inv(ρ);
	v = iρ * s;
	auto const eₖ = half * v.dot(s);
//	auto const eᵢ = E - eₖ;
	auto const eᵢ = ((one - δ₁) * E > eₖ) ? (E - eₖ) : pow(τ, Γ);
	assert(eᵢ > zero);
	p = (Γ - 1) * eᵢ;
	return prim;
}


ConservedGasState prim2con(PrimitiveGasState const &prim) {
	using namespace Gas;
	ConservedGasState con;
	constexpr Real half = inv(Real(2));
	auto const &[_, v, p] = prim;
	auto &[ρ, s, E, τ] = con;
	ρ = prim.ρ;
	s = ρ * v;
	auto const eₖ = half * s.dot(v);
	auto const eᵢ = p / (Γ - 1);
	E = eₖ + eᵢ;
	τ = pow(eᵢ, inv(Γ));
	return con;
}


ConservedGasState prim2flux(PrimitiveGasState const &prim, int n) {
	using namespace Gas;
	ConservedGasState F;
	constexpr Real half = inv(Real(2));
	auto const &[ρ, v, p] = prim;
	auto const vₙ = v[n];
	auto const eᵢ = p / (Γ - 1);
	auto const eₖ = half * ρ * v.dot(v);
	F.ρ = vₙ * ρ;
	F.s = vₙ * ρ * v;
	F.s[n] += p; 
	F.E = vₙ * (eₖ + Γ * eᵢ);
	F.τ = vₙ * pow(eᵢ, inv(Γ));
	return F;
}


Real soundSpeed(PrimitiveGasState const &prim) {
	using namespace Gas;
	ConservedGasState f;
	auto const &[ρ, _, p] = prim;
	auto const cₛ = sqrt(Γ * p / ρ);
	return cₛ;
}


Real maxSignalSpeed(PrimitiveGasState const &prim, int n) {
	using std::abs;
	auto const &[ρ, v, p] = prim;
	auto const vₙ = abs(v[n]);
	auto const cₛ = soundSpeed(prim);
	auto const a = cₛ + vₙ;
	return a;
}


ConservedGasState con2flux(ConservedGasState const &Ul, ConservedGasState const &Ur, int n) {
	using std::max;
	ConservedGasState F;
	constexpr Real half = inv(Real(2));
	auto const Vl = con2prim(Ul);
	auto const Fl = prim2flux(Vl, n);
	auto const al = maxSignalSpeed(Vl, n);
	auto const Vr = con2prim(Ur);
	auto const Fr = prim2flux(Vr, n);
	auto const ar = maxSignalSpeed(Vr, n);
	auto const a = max(al, ar);
	F = half * ((Fl + Fr) - a * (Ur - Ul));
	return F;
}

#endif /* GAS_HPP_ */
