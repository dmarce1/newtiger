/*
 * StateFactory.hpp
 *
 *  Created on: Feb 10, 2026
 *      Author: dmarce1
 */

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#include <config.hpp>

#include "SimdVector.hpp"
#include "Vector.hpp"

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

template<size_t N>
struct SoA {
#define SCALAR(arg) SimdVector<Real, N> arg;
#define VECTOR(arg) Vector<SimdVector<Real, N>, dimCount> arg;
	STATE_DATA
#undef SCALAR
#undef VECTOR
	SoA() = default;
	constexpr size_t size() const {
		return N;
	}
	SoA cshift(int shift) const {
		SoA result;
#define SCALAR(arg) result. arg = arg .cshift(shift);
#define VECTOR(arg)                      \
	for (int d = 0; d < dimCount; d++) { \
		result. arg [d] = arg [d].cshift(shift);  \
	}
		STATE_DATA
		return result;
#undef SCALAR
#undef VECTOR
	}
	SoA operator-(SoA const& other) const {
		SoA result;
#define SCALAR(arg) result. arg = arg - other. arg;
#define VECTOR(arg)                                  \
	for (int d = 0; d < dimCount; d++) {             \
		result. arg [d] = arg [d] - other. arg [d];  \
	}
		STATE_DATA
		return result;
#undef SCALAR
#undef VECTOR
	}
	SoA operator+(SoA const &other) const {
		SoA result;
#define SCALAR(arg) result.arg = arg + other.arg;
#define VECTOR(arg)                                                                                                                        \
		for (int d = 0; d < dimCount; d++) {                                                                                                   \
			result.arg[d] = arg[d] + other.arg[d];                                                                                             \
		}
		STATE_DATA
		return result;
	#undef SCALAR
	#undef VECTOR
	}
	SoA operator*(SimdVector<Real, N> const &other) const {
		SoA result;
#define SCALAR(arg) result.arg = arg * other;
#define VECTOR(arg)                                                                                                                        \
		for (int d = 0; d < dimCount; d++) {                                                                                               \
			result.arg[d] = arg[d] * other;                                                                                         \
		}
		STATE_DATA
		return result;
#undef SCALAR
#undef VECTOR
	}
	SoA operator*(Real scalar) const {
		SoA result;
#define SCALAR(arg) result.arg = arg * scalar;
#define VECTOR(arg)                                                                                                                        \
		for (int d = 0; d < dimCount; d++) {                                                                                                   \
			result.arg[d] = arg[d] * scalar;                                                                                             \
		}
		STATE_DATA
		return result;
#undef SCALAR
#undef VECTOR
	}
	SoA &operator+=(SoA const &other) {
#define SCALAR(arg) arg += other.arg;
#define VECTOR(arg)                                                                                                                        \
		for (int d = 0; d < dimCount; d++) {                                                                                                   \
			arg[d] += other.arg[d];                                                                                                            \
		}
		STATE_DATA
		return *this;
#undef SCALAR
#undef VECTOR
	}
	SoA &operator-=(SoA const &other) {
#define SCALAR(arg) arg -= other.arg;
#define VECTOR(arg)                                                                                                                        \
		for (int d = 0; d < dimCount; d++) {                                                                                                   \
			arg[d] -= other.arg[d];                                                                                                            \
		}
		STATE_DATA
		return *this;
#undef SCALAR
#undef VECTOR
	}
	friend SoA operator*(Real scalar, SoA const &soa) {
		return soa * scalar;
	}
	friend SoA operator*(SimdVector<Real, N> scalar, SoA const &soa) {
		return soa * scalar;
	}
};

template<size_t size>
using AoS = std::array<State, size>;

template<size_t N>
SoA<N> aos2soa(AoS<N> const& aos) {
	SoA<N> soa;
#define SCALAR(arg)                     \
	for(size_t i = 0; i < N; i++) { \
		soa. arg [i] = aos[i]. arg ;    \
	}
#define VECTOR(arg)                            \
	for(size_t i = 0; i < N; i++) {        \
		for(size_t d = 0; d < dimCount; d++) { \
			soa. arg [d][i] = aos[i]. arg [d]; \
		}                                      \
	}
	STATE_DATA
#undef SCALAR
#undef VECTOR
	return soa;
}

template<size_t N>
AoS<N> soa2aos(SoA<N> const& soa) {
	AoS<N> aos;
#define SCALAR(arg)                     \
	for(size_t i = 0; i < N; i++) { \
		aos[i]. arg = soa. arg [i] ;    \
	}
#define VECTOR(arg)                            \
	for(size_t i = 0; i < N; i++) {        \
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