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
	static constexpr Integer eidx = ρidx + 1;
	static constexpr Integer midx = eidx + 1;
	static constexpr Integer τidx = midx + dimCount;

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
			strings[τidx] = "entropyTracer";
			return strings;
		}();
		return names;
	}
	Type internalEnergy(Type const &ek) const {
		auto const ei = max(e - ek, 0_R);
		// return ei;
		return select(ei >= 1e-3_R * e, ei, pow(τ, Γ));
	}
	GasState updateEntropy() const {
		GasState U = *this;
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const ei = e - 0.5_R * ρ * u.dot(u);
		U.τ = select(ei < 1e-1_R * e, τ, pow(ei, inv(Γ)));
		return U;
	}
	GasState &setEntropy() {
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const ek = 0.5_R * ρ * u.dot(u);
		auto const ei = e - ek;
		τ = pow(ei, inv(Γ));
		return *this;
	}
	GasState() :
		ρ(*(new(&(*this)[ρidx]) Type{})),
		m(*(new(&(*this)[midx]) Vector<Type, dimCount>{})),
		e(*(new(&(*this)[eidx]) Type{})),
		τ(*(new(&(*this)[τidx]) Type{})) {
	}
	GasState(GasState const &other) :
		ρ(*(new(&(*this)[ρidx]) Type{})),
		m(*(new(&(*this)[midx]) Vector<Type, dimCount>{})),
		e(*(new(&(*this)[eidx]) Type{})),
		τ(*(new(&(*this)[τidx]) Type{})) {
		*this = other;
	}
	GasState(base_type const &other) :
		ρ(*(new(&(*this)[ρidx]) Type{})),
		m(*(new(&(*this)[midx]) Vector<Type, dimCount>{})),
		e(*(new(&(*this)[eidx]) Type{})),
		τ(*(new(&(*this)[τidx]) Type{})) {
		*this = other;
	}
	virtual ~GasState() {
		ρ.~Type();
		e.~Type();
		τ.~Type();
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
	auto eigenstructure(Integer n) const {
		// M 0  0  iM    0 0
		// U I 0  -U*iM  I 0
		// 0 0 I   0     0 I
		std::array<Type, size()> λ;
		constexpr auto ap = size() - 1;
		constexpr auto am = 0;
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const u2 = u.dot(u);
		auto const ek = 0.5_R * ρ * u2;
		auto const ei = internalEnergy(ek);
		auto const p = (Γ - 1_R) * ei;
		auto const a = sqrt(Γ * p * iρ);
		auto const h = iρ * (e + p);
		{
			λ.fill(u[n]);
			λ[am] -= a;
			λ[ap] += a;
		}
		static constexpr Integer N1 = 3;
		static constexpr Integer N2 = dimCount - 1;
		Matrix<Type, N1, N1> M;
		Matrix<Type, N2, N1> U;
		auto const nidx = midx + n;
		M[ρidx][0] = M[ρidx][1] = M[ρidx][2] = 1_R;
		M[midx][0] = M[midx][1] = M[midx][2] = u[n];
		M[midx][0] -= a;
		M[midx][2] += a;
		M[eidx][0] = M[eidx][2] = h;
		M[eidx][0] -= u[n] * a;
		M[eidx][1] = 0.5_R * u2;
		M[eidx][2] += u[n] * a;
		for (Integer i = 0, k = 0; i < N2; i++, k++) {
			if (i == n) k++;
			for (Integer j = 0; j < N1; j++) {
				U[i][j] = u[k];
			}
		}
		auto const iM = inv(M);
		auto const UiM = U * iM;
		auto const L = [iM, UiM, n, nidx](GasState const &x) {
			GasState y;
			for (Integer i = 0; i < N1; i++) {
				y[i] = iM[i][ρidx] * x[ρidx] + iM[i][midx] * x[nidx] + iM[i][eidx] * x[eidx];
			}
			for (Integer i = 0, k = 0; i < N2; i++, k++) {
				if (i == n) k++;
				y[i + N1] = x[midx + k] - UiM[i][ρidx] * x[ρidx] - UiM[i][midx] * x[nidx] - UiM[i][eidx] * x[eidx];
			}
			y[τidx] = x[τidx];
			return y;
		};
		auto const R = [M, U, n, nidx](GasState const &y) {
			GasState x;
			x[ρidx] = M[ρidx][0] * y[0] + M[ρidx][1] * y[1] + M[ρidx][2] * y[2];
			x[nidx] = M[midx][0] * y[0] + M[midx][1] * y[1] + M[midx][2] * y[2];
			x[eidx] = M[eidx][0] * y[0] + M[eidx][1] * y[1] + M[eidx][2] * y[2];
			for (Integer i = 0, k = 0; i < N2; i++, k++) {
				if (i == n) k++;
				x[midx + k] = y[i + N1] - U[ρidx][i] * y[i] - U[nidx][i] * y[i] - U[eidx][i] * y[i];
			}
			x[τidx] = y[τidx];
			return x;
		};
		return std::tuple(λ, R, L);
	};
	friend auto riemannFlux(GasState const &UL, GasState const &UR, Integer n) {
		auto const &ρR = UR.ρ;
		auto const &τR = UR.τ;
		auto const &eR = UR.e;
		auto const &ρL = UL.ρ;
		auto const &τL = UL.τ;
		auto const &eL = UL.e;
		auto const iρR = inv(ρR);
		auto const uR = iρR * UR.m;
		auto const eiR = UR.internalEnergy(0.5_R * ρR * uR.dot(uR));
		auto const pR = (Γ - 1_R) * eiR;
		auto const cR = sqrt(Γ * pR * iρR);
		auto const iρL = inv(ρL);
		auto const uL = iρL * UL.m;
		auto const eiL = UL.internalEnergy(0.5_R * ρL * uL.dot(uL));
		auto const pL = (Γ - 1_R) * eiL;
		auto const cL = sqrt(Γ * pL * iρL);
		auto const sR = max(0_R, max(uR[n] + cR, uL[n] + cL));
		auto const sL = min(0_R, min(uR[n] - cR, uL[n] - cL));
		auto const wL = ρL * (sL - uL[n]);
		auto const wR = ρR * (sR - uR[n]);
		auto const s0 = ((pR - pL) + (wL * uL[n] - wR * uR[n])) / (wL - wR);
		auto const f = s0 > 0_R;
		auto const ρK = select(f, ρL, ρR);
		auto const τK = select(f, τL, τR);
		auto const pK = select(f, pL, pR);
		auto const sK = select(f, sL, sR);
		auto const eK = select(f, eL, eR);
		auto const wK = select(f, wL, wR);
		auto const uK = select(f, uL, uR);
		auto const UK = select(f, UL, UR);
		auto const w0 = (sK - uK[n]) / (sK - s0);
		GasState U0;
		U0.ρ = w0 * ρK;
		U0.τ = w0 * τK;
		U0.m = U0.ρ * uK;
		U0.e = U0.ρ * (eK / ρK + (s0 - uK[n]) * (s0 + pK / wK));
		U0.m[n] = U0.ρ * s0;
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
		τ = pow(ei, inv(Γ));
	}
	friend std::ostream &operator<<(std::ostream &os, GasState const &u) {
		os << "(ρ=" << u.ρ << ", ";
		for (Integer i = 0; i < dimCount; i++) {
			os << "m_" << std::string(1, 'x' + i) << u.m[i] << ", ";
		}
		os << "e=" << u.e << ", ";
		os << "tau=" << u.τ << ")";
		return os;
	}
	template <typename, Integer>
	friend struct CoupledState;

	Type &ρ;
	Vector<Type, dimCount> &m;
	Type &e;
	Type &τ;
};

#endif /* GasState_HPP_ */
