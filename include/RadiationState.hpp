/*
 * RadiationState.hpp
 *
 *  Created on: Feb 25, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_RADIATIONSTATE_HPP_
#define INCLUDE_RADIATIONSTATE_HPP_

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#include <string>
#include <tuple>

#include "Concepts.hpp"
#include "Debug.hpp"
#include "GasState.hpp"
#include "Integer.hpp"
#include "Matrix.hpp"
#include "Real.hpp"
#include "Select.hpp"
#include "Vector.hpp"

#include <numeric>
#include <sstream>

template <typename Type, Integer dimCount>
struct RadiationState : Vector<Type, 1 + dimCount> {
private:
	static constexpr Integer eidx = 0;
	static constexpr Integer fidx = 1;
	static constexpr auto zero = 0.0_R;
	static constexpr auto half = 0.5_R;
	static constexpr auto one = 1.0_R;

public:
	static constexpr Integer size() {
		return dimCount + 1;
	}
	using base_type = Vector<Type, size()>;
	static auto getFieldNames() {
		static auto names = []() {
			static std::array<std::string, size()> strings;
			for (Integer d = 0; d < dimCount; d++) {
				strings[fidx + d] = std::string(1, 'x' + d) + std::string("RadiationFluxDensity");
			}
			strings[eidx] = "radiationEnergyDensity";
			return strings;
		}();
		return names;
	}
	RadiationState() :
		E(*(new(&(*this)[0]) Type{})), F(*(new(&(*this)[1]) Vector<Type, dimCount>{})) {
	}
	RadiationState(RadiationState const &other) :
		E(*(new(&(*this)[0]) Type{})), F(*(new(&(*this)[1]) Vector<Type, dimCount>{})) {
		*this = other;
	}
	RadiationState(base_type const &other) :
		E(*(new(&(*this)[0]) Type{})), F(*(new(&(*this)[1]) Vector<Type, dimCount>{})) {
		*this = other;
	}
	virtual ~RadiationState() {
		E.~Type();
		F.~Vector<Type, dimCount>();
	}
	RadiationState &operator=(RadiationState const &other) {
		base_type::operator=(other);
		return *this;
	}
	RadiationState &operator=(base_type const &other) {
		base_type::operator=(other);
		return *this;
	}
	auto eigenstructure(Integer dim) const {
		if (dim != dimCount - 1) std::swap(F[dimCount - 1], F[dim]);
		ASSERT_POSITIVE(E);
		auto const iE = inv(E);
		auto const f2 = F.dot(F) * sqr(iE);
		ASSERT_RANGE(0_R, f2, 1_R);
		auto const Δ = 4_R - 3_R * f2;
		auto const H = (1_R / 3_R) * (2_R + sqrt(Δ)) * E;
		auto const H2 = sqr(H);
		auto const iH = inv(H);
		auto β = F * iH;
		if constexpr (dimCount > 1) {
			constexpr auto eps2 = sqr(Real(dimCount - 1) * eps_R);
			auto const b2 = std::inner_product(β.begin(), β.begin() + dimCount - 1, β.begin(), 0_R);
			if (b2 < eps2 * H2) {
				for (Integer d = 0; d < dimCount - 1; d++) {
					β[d] = eps_R * H;
				}
			}
		}
		auto const β2 = β.dot(β);
		auto const n = F.normalize();
		auto const &βz = β[dimCount - 1];
		auto const βz2 = sqr(βz);
		auto const λs = 2_R * βz;
		auto const den = 3_R - β2;
		auto const num2 = (1_R - β2) * (den - 2_R * βz2);
		auto const λd = sqrt(num2) * inv(den);
		Vector<Type, size()> λ;
		Matrix<Type, size()> R, L;
		λ.front() = λs - λd;
		λ.back() = λs + λd;
		for (Integer d = 1; d < dimCount; d++) {
			λ[d] = βz;
		}
		auto const λ14 = [&](auto const &λ) {
			Vector<Type, size()> r;
			r.front() = 1_R - 2_R * βz * λ + sqr(λ);
			r.back() = -βz + 2_R * λ - βz * sqr(λ);
			for (Integer d = 0; d < dimCount - 1; d++) {
				r[d + 1] = β[d] * (1_R - sqr(λ));
			}
			return r;
		};
		R[0] = λ14(λ.front());
		R[dimCount] = λ14(λ.back());
		if constexpr (dimCount > 1) {
			R[1].front() = sqrt(β2 - βz2);
			R[1].back() = R[1].front() * βz;
			for (Integer d = 0; d < dimCount - 1; d++) {
				R[1][d + 1] = n[d] * (1_R - βz2);
			}
		}
		if constexpr (dimCount > 2) {
			R[2].front() = R[2].back() = zero;
			R[2][1] = -n[1];
			R[2][2] = +n[0];
		}
		R = transpose(R);
		if (dim != dimCount - 1) {
			for (Integer d = 0; d < dimCount; d++) {
				std::swap(R[dim][d + 1], R[dimCount][d]);
				std::swap(L[d + 1][dim], L[d][dimCount]);
			}
		}
		L = inv(R);
		return std::tuple(λ, R, L);
	}
	friend auto riemannFlux(RadiationState const &UL, RadiationState const &UR, Integer n) {
		auto const &ER = UR.E;
		auto const &FR = UR.F;
		auto const iER = inv(ER);
		auto const F2R = FR.dot(FR);
		auto const f2R = F2R * sqr(iER);
		ASSERT_RANGE(0_R, f2R, 1_R);
		auto const ΔR = 4_R - 3_R * f2R;
		auto const sqrtΔR = sqrt(4_R - 3_R * f2R);
		auto const χR = (1_R / 3_R) * (5_R - 2_R * sqrtΔR);
		auto const iF2R = inv(select(F2R != 0_R, F2R, 1_R));
		auto const iFR = sqrt(iF2R);
		auto const βR = 2_R * inv(3_R - χR) * FR[n] * iER;
		auto const ΠR = 0.5_R * (1_R - χR) * ER[n];
		auto const μR = FR[n] * iFR;
		auto const fR = sqrt(f2R);
		auto const ζR = sqrt((2_R / 3_R) * (ΔR - sqrtΔR) + 2_R * sqr(μR) * (2_R - f2R - sqrtΔR));
		auto const λmR = (μR * fR - ζR) * inv(sqrtΔR);
		auto const λpR = (μR * fR + ζR) * inv(sqrtΔR);

		auto const &EL = UL.E;
		auto const &FL = UL.F;
		auto const iEL = inv(EL);
		auto const F2L = FL.dot(FL);
		auto const f2L = F2L * sqr(iEL);
		ASSERT_RANGE(0_R, f2L, 1_R);
		auto const ΔL = 4_R - 3_R * f2L;
		auto const sqrtΔL = 4_R - 3_R * f2L;
		auto const χL = (1_R / 3_R) * (5_R - 2_R * sqrtΔL);
		auto const iF2L = inv(select(F2L != 0_R, F2L, 1_R));
		auto const iFL = sqrt(iF2L);
		auto const βL = 2_R * inv(3_R - χL) * FL[n] * iEL;
		auto const ΠL = 0.5_R * (1_R - χL) * EL;
		auto const μL = FL[n] * iFL;
		auto const fL = sqrt(f2L);
		auto const ζL = sqrt((2_R / 3_R) * (ΔL - sqrtΔL) + 2_R * sqr(μL) * (2_R - f2L - sqrtΔL));
		auto const λmL = (μL * fL - ζL) * inv(sqrtΔL);
		auto const λpL = (μL * fL + ζL) * inv(sqrtΔL);

		auto const λR = max(0_R, max(λpL, λpR));
		auto const λL = min(0_R, min(λmL, λmR));

		auto const AR = λR * ER - FR[n];
		auto const AL = λL * EL - FL[n];
		auto const BR = (λR - βR) * FR[n] - ΠR;
		auto const BL = (λL - βL) * FL[n] - ΠL;

		auto const a = AR * λL - AL * λR;
		auto const i2a = inv(2_R * a);
		auto const b = AL - AR + BL * λR - BR * λL;
		auto const disc = b * b - 2_R * (BR - BL);
		ASSERT_NONNEGATIVE(disc);
		auto const sqrtDisc = sqrt(disc);
		auto const r1 = -(sqrtDisc + b) * i2a;
		auto const r2 = +(sqrtDisc - b) * i2a;
		auto const λ0 = select(abs(r1) < abs(r2), r1, r2);
		ASSERT_RANGE(-1_R, λ0, 1_R);

		auto const flag = λ0 > 0_R;
		auto const AK = select(flag, AL, AR);
		auto const BK = select(flag, BL, BR);
		auto const λK = select(flag, λL, λR);
		auto const EK = select(flag, EL, ER);
		auto const FK = select(flag, FL, FR);
		auto const βK = select(flag, βL, βR);
		auto const ΠK = select(flag, ΠL, ΠR);

		auto const Π0 = (AK * λ0 - BK) / (1_R - λK * λ0);
		auto const λKmλ0 = λK - λ0;
		ASSERT_NONZERO(λKmλ0);
		auto const iden = inv(λKmλ0);
		auto const E0 = iden * (EK * (λK - βK) + Π0 * λ0 - ΠK * βK);
		auto F0 = iden * (FK * (λK - βK));
		F0[n] += Π0 - ΠK;
		auto const F20 = F0.dot(F0);
		auto const iF20 = inv(select(F20 != 0_R, F20, 1_R));
		auto const β0 = (1_R - 3_R * Π0 / E0) * F0[n] * iF20;
		RadiationState flux;
		flux.E = F0[n];
		flux.F = β0 * F0[n];
		flux.F += Π0;
		return std::pair(flux, max(λR, -λL));
	}
	void setEnergy(Type e) {
		E = e;
	}
	void setFlux(Vector<Type, dimCount> const &f) {
		F = f;
	}
	template <typename, Integer>
	friend struct CoupledState;
	friend std::ostream &operator<<(std::ostream &os, RadiationState const &u) {
		os << "(E=" << u.E << ", ";
		for (Integer i = 0; i < dimCount; i++) {
			os << "F_" << std::string(1, 'x' + i) << u.F[i] << ", ";
		}
		os << ")";
		return os;
	}

	Type &E;
	Vector<Type, dimCount> &F;
};

#endif /* INCLUDE_RADIATIONSTATE_HPP_ */
