/*
 * Gas.hpp
 *
 *  Created on: Feb 9, 2026
 *      Author: dmarce1
 */

#ifndef GAS_HPP_
#define GAS_HPP_

#include "Vector.hpp"
/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
struct ConservedGasState {
	Type ρ;
	Vector<Type, nDim> s;
	Type E;
	Type τ;
};

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
ConservedGasState<Type, nDim, Γ> operator+(ConservedGasState<Type, nDim, Γ> u, ConservedGasState<Type, nDim, Γ> const &v) {
	ConservedGasState<Type, nDim, Γ> w;
	w.ρ = u.ρ + v.ρ;
	w.s = u.s + v.s;
	w.E = u.E + v.E;
	w.τ = u.τ + v.τ;
	return w;
}

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
ConservedGasState<Type, nDim, Γ> operator-(ConservedGasState<Type, nDim, Γ> u, ConservedGasState<Type, nDim, Γ> const &v) {
	ConservedGasState<Type, nDim, Γ> w;
	w.ρ = u.ρ - v.ρ;
	w.s = u.s - v.s;
	w.E = u.E - v.E;
	w.τ = u.τ - v.τ;
	return w;
}

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
ConservedGasState<Type, nDim, Γ> &operator-=(ConservedGasState<Type, nDim, Γ> &u, ConservedGasState<Type, nDim, Γ> const &v) {
	u.ρ -= v.ρ;
	u.s -= v.s;
	u.E -= v.E;
	u.τ -= v.τ;
	return u;
}

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
ConservedGasState<Type, nDim, Γ> operator*(Scalar auto a, ConservedGasState<Type, nDim, Γ> const &u) {
	ConservedGasState<Type, nDim, Γ> v;
	v.ρ = a * u.ρ;
	v.s = a * u.s;
	v.E = a * u.E;
	v.τ = a * u.τ;
	return v;
}

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
struct PrimitiveGasState {
	Type ρ;
	Vector<Type, nDim> v;
	Type p;
};

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
PrimitiveGasState<Type, nDim, Γ> con2prim(ConservedGasState<Type, nDim, Γ> const &con) {
	PrimitiveGasState<Type, nDim, Γ> prim;
	constexpr Type zero = Type(0);
	constexpr Type one = Type(1);
	constexpr Type half = inv(Type(2));
	constexpr Type δ₁ = inv(Type(1000));
	auto const &[_, s, E, τ] = con;
	auto &[ρ, v, p] = prim;
	ρ = con.ρ;
	assert(ρ > zero);
	assert(τ > zero);
	auto const iρ = inv(ρ);
	v = iρ * s;
	auto const eₖ = half * v.dot(s);
	auto const eᵢ = ((one - δ₁) * E < eₖ) ? (E - eₖ) : pow(τ, Γ);
	assert(eᵢ > zero);
	p = (Γ - 1) * eᵢ;
	return prim;
}

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
ConservedGasState<Type, nDim, Γ> prim2con(PrimitiveGasState<Type, nDim, Γ> const &prim) {
	ConservedGasState<Type, nDim, Γ> con;
	constexpr Type half = inv(Type(2));
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

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
ConservedGasState<Type, nDim, Γ> prim2flux(PrimitiveGasState<Type, nDim, Γ> const &prim, int n) {
	ConservedGasState<Type, nDim, Γ> F;
	constexpr Type half = inv(Type(2));
	auto const &[ρ, v, p] = prim;
	auto const vₙ = v[n];
	auto const eᵢ = p / (Γ - 1);
	auto const eₖ = half * ρ * v.dot(v);
	F.ρ = vₙ * ρ;
	F.s = vₙ * ρ * v;
	F.E = vₙ * (eₖ + Γ * eᵢ);
	F.τ = pow(eᵢ, inv(Γ));
	return F;
}

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
Type soundSpeed(PrimitiveGasState<Type, nDim, Γ> const &prim) {
	ConservedGasState<Type, nDim, Γ> f;
	auto const &[ρ, _, p] = prim;
	auto const cₛ = sqrt(Γ * p / ρ);
	return cₛ;
}

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
Type maxSignalSpeed(PrimitiveGasState<Type, nDim, Γ> const &prim, int n) {
	using std::abs;
	auto const &[ρ, v, p] = prim;
	auto const vₙ = abs(v[n]);
	auto const cₛ = soundSpeed(prim);
	auto const a = cₛ + vₙ;
	return a;
}

template <typename Type = double, int nDim = 1, Scalar auto Γ = Type(5) / Type(3)>
ConservedGasState<Type, nDim, Γ> con2flux(ConservedGasState<Type, nDim, Γ> const &Ul, ConservedGasState<Type, nDim, Γ> const &Ur, int n) {
	using std::max;
	ConservedGasState<Type, nDim, Γ> F;
	constexpr Type half = inv(Type(2));
	auto const Vl = con2prim(Ul);
	auto const Fl = prim2flux(Vl, n);
	auto const aₗ = maxSignalSpeed(Vl, n);
	auto const Vr = con2prim(Ur);
	auto const Fr = prim2flux(Vr, n);
	auto const aᵣ = maxSignalSpeed(Vr, n);
	auto const a = max(aᵣ, aₗ);
	F = half * ((Fl + Fr) - a * (Ur - Ul));
	return F;
}

#endif /* GAS_HPP_ */
