/*
 * StateFactory.hpp
 *
 *  Created on: Feb 10, 2026
 *      Author: dmarce1
 */

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#include "Gas.hpp"
#include "Vector.hpp"
#include <config.hpp>
#include <valarray>

#ifndef STATE_DATA
#error "STATE_DATA not defined"
#endif

#ifndef STATE_NAME
#error "STATE_NAME not defined"
#endif

namespace STATE_NAME {

using namespace config;

#define CONSTANT(name, value) static constexpr Real name = Real(value);
STATE_CONSTANTS
#undef CONSANT

struct State {
	State &operator+=(State const &other) {
#define SCALAR(arg) arg += other.arg;
#define VECTOR(arg) arg += other.arg;
		STATE_DATA
#undef SCALAR
#undef VECTOR
		return *this;
	}
	State &operator-=(State const &other) {
#define SCALAR(arg) arg -= other.arg;
#define VECTOR(arg) arg -= other.arg;
		STATE_DATA
#undef SCALAR
#undef VECTOR
		return *this;
	}
	State &operator*=(Scalar auto scale) {
#define SCALAR(arg) arg *= scale;
#define VECTOR(arg) arg *= scale;
		STATE_DATA
#undef SCALAR
#undef VECTOR
		return *this;
	}
	State operator+(State const &other) const {
		State result;
#define SCALAR(arg) result.arg = arg + other.arg;
#define VECTOR(arg) result.arg = arg + other.arg;
		STATE_DATA
#undef SCALAR
#undef VECTOR
		return result;
	}
	State operator-(State const &other) const {
		State result;
#define SCALAR(arg) result.arg = arg - other.arg;
#define VECTOR(arg) result.arg = arg - other.arg;
		STATE_DATA
#undef SCALAR
#undef VECTOR
		return result;
	}
	State operator*(Scalar auto scale) const {
		State result;
#define SCALAR(arg) result.arg = arg * scale;
#define VECTOR(arg) result.arg = arg * scale;
		STATE_DATA
#undef SCALAR
#undef VECTOR
		return result;
	}
	friend State operator*(Scalar auto scale, State const& state) {
		return state * scale;
	}
#define SCALAR(arg) Real arg;
#define VECTOR(arg) Vector<Real, dimCount> arg;
	STATE_DATA
#undef SCALAR
#undef VECTOR
};

struct SoA {
#define SCALAR(arg) std::valarray<Real> arg;
#define VECTOR(arg) Vector<std::valarray<Real>, dimCount> arg;
	STATE_DATA
#undef SCALAR
#undef VECTOR
	size_t size_;
	SoA() = default;
	SoA(size_t size) {
		resize(size);
	}
	size_t size() const {
		return size_;
	}
	void resize(size_t size) {
		size_ = size;
#define SCALAR(arg) arg.resize(size);
#define VECTOR(arg)       \
	for (auto& x : arg) { \
		x.resize(size);   \
	}
		STATE_DATA
#undef SCALAR
#undef VECTOR
	}
};

template<int N>
using AoS = std::array<State, N>;

SoA aos2soa(AoS const& aos) {
	SoA soa(aos.size());
	size_t const count = aos.size();
#define SCALAR(arg)                     \
	for(size_t i = 0; i < count; i++) { \
		soa. arg [i] = aos[i]. arg ;    \
	}
#define VECTOR(arg)                            \
	for(size_t i = 0; i < count; i++) {        \
		for(size_t d = 0; d < dimCount; d++) { \
			soa. arg [d][i] = aos[i]. arg [d]; \
		}                                      \
	}
	STATE_DATA
#undef SCALAR
#undef VECTOR
	return soa;
}

AoS soa2aos(SoA const& soa) {
	AoS aos(soa.size());
	size_t const count = aos.size();
#define SCALAR(arg)                     \
	for(size_t i = 0; i < count; i++) { \
		aos[i]. arg = soa. arg [i] ;    \
	}
#define VECTOR(arg)                            \
	for(size_t i = 0; i < count; i++) {        \
		for(size_t d = 0; d < dimCount; d++) { \
			aos[i]. arg [d] = soa. arg [d][i]; \
		}                                      \
	}
	STATE_DATA
#undef SCALAR
#undef VECTOR
	return aos;
}



} // namespace STATE_NAME

#undef CONSTANT
#undef STATE_DATA
#undef STATE_NAME