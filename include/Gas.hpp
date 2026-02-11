///*
// * Gas.hpp
// *
// *  Created on: Feb 9, 2026
// *      Author: dmarce1
// */
//
//#ifndef GAS_HPP_
//#define GAS_HPP_
//
//#include <cassert>
//#include <cmath>
//#include <type_traits>
///*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
///*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/
//
//#define STATE_NAME Gas
//
//#define STATE_DATA                                                                                                                         \
//	SCALAR(ρ)                                                                                                                              \
//	VECTOR(s)                                                                                                                              \
//	SCALAR(E)                                                                                                                              \
//	SCALAR(τ)
//
//#define STATE_CONSTANTS                                                                                                                    \
//	CONSTANT(Γ, fluidGamma)                                                                                                                \
//	CONSTANT(μ, meanMolecularWeight)                                                                                                       \
//	CONSTANT(δ₁, dualEnergySwitch1)                                                                                                        \
//	CONSTANT(δ₂, dualEnergySwitch2)
//
//#include "StateFactory.hpp"
//
//using namespace config;
//
//
//template <typename GasState>
//GasState con2flux(GasState const &con, int n) {
//	using namespace Gas;
//	GasState F;
//	auto const &[ρ, s, E, τ] = con;
//	auto const v = s / ρ;
//	auto const eₖ = Real(0.5) * v.dot(s);
////	auto const flag = (Real(1) - δ₁) * E > eₖ;
////	auto const p = (Γ - Real(1)) * (flag * (E - eₖ) + (Real(1) - flag) * pow(τ, Γ));
//	auto const p = (Γ - Real(1)) * (E - eₖ);
//	F.ρ = v[n] * ρ;
//	F.s = v[n] * s;
//	F.E = v[n] * (E + p);
//	F.τ = v[n] * τ;
//	F.s[n] += p;
//	return F;
//}
//
//template <typename GasState>
//auto maxSignalSpeed(GasState const &con, int n) {
//	using namespace Gas;
//	using std::abs;
//	using std::sqrt;
//	using Scalar = std::remove_cv_t<decltype(con.ρ)>;
//	auto const [ρ, s, E, τ] = con;
//	auto const v = s / ρ;
//	auto const eₖ = Real(0.5) * v.dot(s);
////	Scalar const flag = Scalar((Real(1) - δ₁) * E > eₖ);
////	Scalar const p = (Γ - Real(1)) * (flag * (E - eₖ) + (Real(1) - flag) * pow(τ, Γ));
//	auto const p = (Γ - Real(1)) * (E - eₖ);
//	auto const c = sqrt(Γ * p  / ρ);
//	auto const a = c + abs(v[n]);
//	return Scalar(a);
//}
//
//template <typename GasState>
//GasState con2flux(GasState const &Ul, GasState const &Ur, int n) {
//	using std::abs;
//	using std::max;
//	GasState F;
//	auto const Fl = con2flux(Ul, n);
//	auto const Fr = con2flux(Ur, n);
//	auto const al = maxSignalSpeed(Ul, n);
//	auto const ar = maxSignalSpeed(Ur, n);
//	auto const a = Real(0.5) * (abs(al) + abs(ar) + abs(al - ar));
//	F = Real(0.5) * ((Fl + Fr) - a * (Ur - Ul));
//	return F;
//}
//
//template <typename GasState>
//void dualEnergyUpdate(GasState &con) {
//// 	using std::pow;
////	auto &[ρ, s, E, τ] = con;
////	auto const eₖ = Real(0.5) * inv(ρ) * s.dot(s);
////	auto const flag = (Real(1) - δ₂) * E > eₖ;
////	τ = flag * pow(E - eₖ, inv(Γ)) + (Real(1) - flag) * τ;
//}
//
//#endif /* GAS_HPP_ */
