/*
 * Constants.hpp
 *
 *  Created on: Feb 25, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_CONSTANTS_HPP_
#define INCLUDE_CONSTANTS_HPP_

#include "Math.hpp"
#include "Real.hpp"
#include "Units.hpp"
/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#define UNITS_CGS

namespace Constants {
inline constexpr auto c = 2.99792458e10_cm / 1_s;
inline constexpr auto kB = 1.380649e-16_erg / 1_K;
inline constexpr auto G = 6.67430e-8_cm3 / (1_g * 1_s2);
inline constexpr auto mᵤ = 1.66053906892e-24_g;
inline constexpr auto σ = 5.67037441918442945397e-5_erg / (1_s * 1_cm2 * 1_K4);
inline constexpr auto π = 3.14159265358979324_R;
inline constexpr auto Rsun = 6.957e10_cm;
inline constexpr auto Lsun = 3.828e33_erg / 1_s;
inline constexpr auto Msun = 1.9884098707e33_g;

} // namespace Constants

struct CodeUnits {
	constexpr CodeUnits(UnitsType auto const &c0, UnitsType auto const &c1, UnitsType auto const &c2, UnitsType auto const &c3) {
		Matrix<Rational, 4> C;
		Vector<Real, 4> v, x;
		C[0] = c0.powers();
		C[1] = c1.powers();
		C[2] = c2.powers();
		C[3] = c3.powers();
		Matrix<Rational, 4> const iC = inv(C);
		v[0] = c0.get();
		v[1] = c1.get();
		v[2] = c2.get();
		v[3] = c3.get();
		for (Integer i = 0; i < 4; i++) {
			x[i] = 1_R;
			for (Integer k = 0; k < 4; k++) {
				x[i] *= pow(v[k], iC[i][k]);
			}
		}
		code2cm = x[0];
		code2g = x[1];
		code2s = x[2];
		code2K = x[3];
		π = Constants::π;
		c = Constants::c.get() * inv(code2cm / code2s);
		σ = Constants::σ.get() * inv(code2g / (code2s * sqr(code2s) * sqr(sqr(code2K))));
		kB = Constants::kB.get() * inv(code2g * sqr(code2cm) / (sqr(code2s) * code2K));
		mᵤ = Constants::mᵤ.get() * inv(code2g);
		G = Constants::G.get() * inv(code2cm * sqr(code2cm) / (code2g * sqr(code2s)));
	}
	Real π;
	Real c;
	Real kB;
	Real G;
	Real mᵤ;
	Real σ;
	Real code2cm;
	Real code2s;
	Real code2g;
	Real code2K;
};

static auto consteval codeUnits() {
	return CodeUnits(Constants::c, Constants::G, Constants::kB / Constants::mᵤ, Constants::σ);
}
#endif /* INCLUDE_CONSTANTS_HPP_ */
