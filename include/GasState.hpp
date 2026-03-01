/*
 * GasState.hpp
 *
 *  Created on: Feb 21, 2026
 *      Author: dmarce1
 */

#ifndef GasState_HPP_
#define GasState_HPP_

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#include <string>
#include <tuple>

#include "Concepts.hpp"
#include "Integer.hpp"
#include "Matrix.hpp"
#include "Real.hpp"
#include "Select.hpp"
#include "Vector.hpp"

template <typename, Integer>
struct RadiationState;

template <typename Type, Integer dimCount>
struct GasState : Vector<Type, 3 + dimCount> {
	static constexpr auto Γ = 5_R / 3_R;
private:
	static constexpr Integer ρidx = 0;
	static constexpr Integer midx = 1;
	static constexpr Integer eidx = 1 + dimCount;
	static constexpr Integer sidx = 2 + dimCount;
	static constexpr auto zero = 0.0_R;
	static constexpr auto half = 0.5_R;
	static constexpr auto one = 1.0_R;
	static constexpr auto Γm1 = (Γ - one);
	static constexpr auto iΓ = inv(Γ);

public:
	static constexpr Integer size() {
		return dimCount + 3;
	}
	using base_type = Vector<Type, size()>;
	static auto getFieldNames() {
		static auto names = []() {
			static std::array<std::string, size()> strings;
			strings[ρidx] = "massDensity";
			for (Integer d = 0; d < dimCount; d++) {
				strings[midx + d] = std::string(1, 'x' + d) + "Momentum";
			}
			strings[eidx] = "gasEnergyDensity";
			strings[sidx] = "entropyTracer";
			return strings;
		}();
		return names;
	}
	Type internalEnergy(Type const &ek) const {
		auto const ei = e - ek;
		return select(ei < 1e-3_R * e, ei, pow(τ, Γ));
	}
	GasState updateInternalEnergy() const {
		GasState U = *this;
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const ek = half * ρ * u.dot(u);
		auto const ei = e - ek;
		U.τ = select(ei < 1e-1_R * e, τ, pow(ei, iΓ));
		return U;
	}
	GasState() :
		ρ(*(new(&(*this)[0]) Type{})),
		m(*(new(&(*this)[1]) Vector<Type, dimCount>{})),
		e(*(new(&(*this)[1 + dimCount]) Type{})),
		τ(*(new(&(*this)[2 + dimCount]) Type{})) {
	}
	GasState(GasState const &other) :
		ρ(*(new(&(*this)[0]) Type{})),
		m(*(new(&(*this)[1]) Vector<Type, dimCount>{})),
		e(*(new(&(*this)[1 + dimCount]) Type{})),
		τ(*(new(&(*this)[2 + dimCount]) Type{})) {
		*this = other;
	}
	GasState(base_type const &other) :
		ρ(*(new(&(*this)[0]) Type{})),
		m(*(new(&(*this)[1]) Vector<Type, dimCount>{})),
		e(*(new(&(*this)[1 + dimCount]) Type{})),
		τ(*(new(&(*this)[2 + dimCount]) Type{})) {
		*this = other;
	}
	virtual ~GasState() {
		ρ.~Type();
		e.~Type();
		m.~Vector<Type, dimCount>();
	}
	GasState &operator=(GasState const &other) {
		base_type::operator=(other);
		return *this;
	}
	GasState &operator=(base_type const &other) {
		base_type::operator=(other);
		return *this;
	}
	auto flux(Type u, Type p, Integer n) const {
		GasState f = *this * u;
		f.m[n] += p;
		f.e += u * p;
		return f;
	};
	std::tuple<std::array<Type, size()>, Matrix<Type, size()>, Matrix<Type, size()>> eigenstructure(Integer n) const {
		Matrix<Type, size()> R, L;
		std::array<Type, size()> λ;
		constexpr auto shearCount = dimCount - 1;
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const ek = half * ρ * u.dot(u);
		auto const ei = internalEnergy(ek);
		auto const p = Γm1 * ei;
		auto const a = sqrt(Γ * p * iρ);
		auto const ia = inv(a);
		auto const h = iρ * (e + p);
		auto const au = a * u[n];
		auto const q2 = u.dot(u);
		auto const β = Γm1 / (a * a);
		{
			λ.fill(u[n]);
			λ[0] -= a;
			λ[dimCount] += a;
		}
		{
			Vector<Type, size()> acousticL, acousticR, contact, tracer;
			std::array<Vector<Type, size()>, dimCount - 1> shear;
			contact[0] = acousticL[0] = acousticR[0] = one;
			for (Integer d = 0; d < dimCount; d++) {
				contact[d + 1] = acousticL[d + 1] = acousticR[d + 1] = u[d];
			}
			acousticR[n + 1] += a;
			acousticL[n + 1] -= a;
			contact[dimCount] = half * q2;
			acousticL[dimCount] = acousticR[dimCount] = h;
			acousticR[dimCount] += au;
			acousticL[dimCount] -= au;
			for (Integer s = 0, ds = 0; s < shearCount; s++, ds++) {
				if (s == n) ds++;
				shear[s][0] = zero;
				for (Integer d = 0; d < dimCount; d++) {
					shear[s][d + 1] = select(d == ds, one, zero);
				}
				shear[s][dimCount] = u[ds];
			}
			contact.back() = acousticL.back() = acousticR.back() = 0_R;
			std::fill(tracer.begin(), tracer.end() - 1, 0_R);
			tracer.back() = 1_R;
			R.setCol(0, acousticL);
			R.setCol(size() - 1, acousticR);
			R.setCol(n + 1, contact);
			for (Integer s = 0, d = 0; s < shearCount; s++) {
				if (d == n) d++;
				R.setCol(d + 1, shear[s]);
			}
			R.setCol(dimCount + 1, tracer);
		}
		{
			Vector<Type, size()> &acousticL = L[0];
			Vector<Type, size()> &acousticR = L[dimCount];
			Vector<Type, size()> &contact = L[n + 1];
			Vector<Type, size()> &tracer = L[dimCount + 1];
			auto const hia = half * ia;
			auto const hβ = half * β;
			auto const hia_un = hia * u[n];
			auto const βq2 = β * q2;
			acousticR[0] = acousticL[0] = half * βq2;
			contact[0] = one - βq2;
			for (Integer d = 0; d < dimCount; d++) {
				auto const βu = β * u[d];
				contact[d + 1] = βu;
				acousticR[d + 1] = acousticL[d + 1] = -half * βu;
			}
			acousticR[dimCount] = acousticL[dimCount] = hβ;
			contact[dimCount] = -β;
			acousticR[0] -= hia_un;
			acousticL[0] += hia_un;
			acousticR[n + 1] += hia;
			acousticL[n + 1] -= hia;
			for (Integer s = 0, ds = 0; s < shearCount; s++, ds++) {
				if (s == n) ds++;
				Vector<Type, size()> &shear = L[ds];
				shear[0] = -u[ds];
				for (Integer d = 0; d < dimCount; d++) {
					shear[d + 1] = select(d == ds, one, zero);
				}
				shear[dimCount] = zero;
			}
			std::fill(tracer.begin(), tracer.end() - 1, 0_R);
			tracer.back() = 1_R;
		}

		return std::tuple(λ, R, L);
	};
	friend auto riemannFlux(GasState const &UL, GasState const &UR, Integer n) {
		constexpr auto Γm1 = (Γ - one);
		auto const &ρR = UR.ρ;
		auto const &τR = UR.τ;
		auto const &eR = UR.e;
		auto const &ρL = UL.ρ;
		auto const &τL = UL.τ;
		auto const &eL = UL.e;
		auto const iρR = inv(ρR);
		auto const uR = iρR * UR.m;
		auto const eiR = UR.internalEnergy(half * ρR * uR.dot(uR));
		auto const pR = Γm1 * eiR;
		auto const cR = sqrt(Γ * pR * iρR);
		auto const iρL = inv(ρL);
		auto const uL = iρL * UL.m;
		auto const eiL = UL.internalEnergy(half * ρL * uL.dot(uL));
		auto const pL = Γm1 * eiL;
		auto const cL = sqrt(Γ * pL * iρL);
		auto const sR = max(zero, max(uR[n] + cR, uL[n] + cL));
		auto const sL = min(zero, min(uR[n] - cR, uL[n] - cL));
		auto const wL = ρL * (sL - uL[n]);
		auto const wR = ρR * (sR - uR[n]);
		auto const s0 = ((pR - pL) + (wL * uL[n] - wR * uR[n])) / (wL - wR);
		auto const f = s0 > zero;
		auto const ρK = select(f, ρL, ρR);
		auto const τK = select(f, τL, τR);
		auto const pK = select(f, pL, pR);
		auto const sK = select(f, sL, sR);
		auto const eK = select(f, eL, eR);
		auto const wK = select(f, wL, wR);
		auto const uK = select(f, uL, uR);
		auto const UK = select(f, UL, UR);
		auto const w0 = wK / (sK - s0);
		auto const e0 = w0 * (eK / ρK + (s0 - uK[n]) * (s0 + pK / wK));
		GasState U0;
		U0.ρ = w0;
		U0.m = w0 * uK;
		U0.τ = w0 * τK;
		U0.m[n] = w0 * s0;
		U0.e = e0;
		auto const FK = UK.flux(uK[n], pK, n);
		auto const F = FK + sK * (U0 - UK);
		return std::pair(F, max(sR, -sL));
	}
	void setDensity(Type d) {
		ρ = d;
	}
	void setPressure(Type p) {
		auto const ei = p * inv(Γ - 1_R);
		e = ei + 0.5_R * inv(ρ) * m.dot(m);
		τ = pow(ei, iΓ);
	}
	friend std::ostream &operator<<(std::ostream &os, GasState const &u) {
		os << "(ρ=" << u.ρ << ", ";
		for (Integer i = 0; i < dimCount; i++) {
			os << "m_" << std::string(1, 'x' + i) << u.m[i] << ", ";
		}
		os << "e=" << u.e << ")";
		os << "tau=" << u.τ << ")";
		return os;
	}
	template<typename, Integer>
	friend struct CoupledState;

	Type &ρ;
	Vector<Type, dimCount> &m;
	Type &e;
	Type &τ;
};

#endif /* GasState_HPP_ */
