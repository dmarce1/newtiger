/*
 * Interval.hpp
 *
 *  Created on: Feb 10, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_INTERVAL_HPP_
#define INCLUDE_INTERVAL_HPP_

#include "Vector.hpp"
/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

template <typename T, int D>
struct Interval {
	constexpr Interval() = default;
	constexpr Interval(T const init) {
		xmin = zero<T>;
		xmax = init;
	}
	constexpr auto begin(int i) const {
		return xmin[i];
	}
	constexpr auto end(int i) const {
		return xmax[i];
	}
	constexpr auto begin() const {
		return xmin;
	}
	constexpr auto end() const {
		return xmax;
	}
	constexpr auto span() const {
		return xmax - xmin;
	}
	constexpr auto span(int d) const {
		return xmax[d] - xmin[d];
	}

private:
	Vector<T, D> xmin{};
	Vector<T, D> xmax{};
};

#endif /* INCLUDE_INTERVAL_HPP_ */
