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

#include "AutoDiff.hpp"
#include "Constants.hpp"
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
	static constexpr auto μ = 4_R / 3_R;
	static constexpr Integer ρidx = 0;
	static constexpr Integer eidx = 1;
	static constexpr Integer midx = 2;
	static constexpr Integer τidx = 2 + dimCount;
	static constexpr auto units = codeUnits();
	static constexpr auto kB = units.kB;
	static constexpr auto mᵤ = units.mᵤ;

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
	auto meanMolecularWeight() const {
		return μ;
	}
	GasState flux(Type u, Type p, Integer n) const {
		GasState f = *this * u;
		f.m[n] += p;
		f.e += u * p;
		return f;
	};
	auto flux(Integer n) const {
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const p = (Γ - 1_R) * internalEnergy(0.5_R * ρ * u.dot(u));
		return flux(u[n], p, n);
	};
	template <Integer count>
	Type enforcePositivity(Vector<GasState, count> const &u) const {
		constexpr auto δ = 1e-3_R;
		using std::abs;
		using std::min;
		using Auto = AutoDiff<Type, 1>;
		auto const δe = δ * internalEnergy();
		Type θ = 1_R;
		for (Integer j = 0; j < count; j++) {
			GasState const du = u[j] - *this;
			{
				auto const δρ = (δ - 1_R) * ρ;
				auto const θn = δρ * inv(select(du.ρ < 0_R, du.ρ, δρ));
				θ = min(θ, θn);
			}
			{
				auto const δτ = (δ - 1_R) * τ;
				auto const θn = δτ * inv(select(du.τ < 0_R, du.τ, δτ));
				θ = min(θ, θn);
			}
			{
				auto const residual = [this, du, δe](auto θ) {
					return e + θ * du.e - 0.5_R * (m.dot(m) + θ * du.m.dot(m) + sqr(θ) * du.m.dot(du.m)) * inv(ρ + θ * du.ρ) - δe;
				};
				if (any(residual(θ) < 0_R)) {
					Type err = 1_R;
					Auto θn = Auto::genVar(θ);
					Vector<Auto, dimCount> m;
					for (Integer i = 0; any(err > 4_R * eps_R); i++) {
						assert(i < 40);
						Auto const f = residual(θn);
						Type const dθn = -f.get() * inv(f.get(0));
						θn = Auto::genVar(θn.get() + dθn);
						err = abs(f.get() * inv(e));
					}
					θ = min(θ, θn.get());
				}
			}
		}
		return θ;
	}
	std::tuple<Vector<Type, size()>, std::function<GasState(GasState const &)>, std::function<GasState(GasState const &)>>
	eigenstructure(Integer n) const {
		FpeGuard fpeGuard{};
		if (n != 0) {
			auto u = *this;
			std::swap(u.m[0], u.m[n]);
			auto const [λ, R, L] = u.eigenstructure(0);
			return std::tuple(
				λ,
				[R, n](GasState const &dc) {
					auto const du = R(dc);
					std::swap(du.m[0], du.m[n]);
					return du;
				},
				[L, n](GasState du) {
					std::swap(du.m[0], du.m[n]);
					auto const dc = L(du);
					return dc;
				});
		}
		constexpr auto am = ρidx;
		constexpr auto ap = eidx;
		constexpr auto ci = midx + 0;
		constexpr auto τi = τidx;
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const u2 = u.dot(u);
		auto const ei = e - 0.5_R * ρ * u2;
		auto const p = (Γ - 1_R) * ei;
		auto const a = sqrt(Γ * p * iρ);
		auto const &un = u[0];
		auto const ua = a * un;
		auto const h = iρ * (e + p);
		Vector<Type, size()> λ;
		{
			λ = u[0];
			λ[am] -= a;
			λ[ap] += a;
		}
		constexpr Integer wCount = 3;
		Matrix<Type, wCount, wCount> M{};
		M[ρidx][am] = M[ρidx][ap] = M[ρidx][ci] = 1_R;
		M[midx][am] = M[midx][ap] = M[midx][ci] = un;
		M[eidx][am] = M[eidx][ap] = h;
		M[eidx][ci] = 0.5_R * u2;
		M[midx][ap] += a;
		M[midx][am] -= a;
		M[eidx][ap] += ua;
		M[eidx][am] -= ua;
		auto const iM = inv(M);
		return std::tuple(
			λ,
			[M, u](GasState const &dc) {
				GasState du{};
				for (Integer j = 0; j < wCount; j++) {
					for (Integer k = 0; k < wCount; k++) {
						du[j] += M[j][k] * dc[k];
					}
				}
				for (Integer j0 = wCount, j1 = 1; j1 < dimCount; j0++, j1++) {
					du[j0] = dc[j0];
					for (Integer k = 0; k < wCount; k++) {
						du[j0] += u[j1] * dc[k];
					}
				}
				du[τidx] = dc[τi];
				return du;
			},
			[iM, u](GasState const &du) {
				GasState dc;
				for (Integer j = 0; j < wCount; j++) {
					for (Integer k = 0; k < wCount; k++) {
						dc[j] += iM[j][k] * du[k];
					}
				}
				for (Integer j0 = wCount, j1 = 1; j1 < dimCount; j0++, j1++) {
					dc[j0] = du[j0] - u[j1] * du[ρidx];
				}
				dc[τi] = du[τidx];
				return dc;
			});
	};
	friend auto riemannFlux(GasState const &UL, GasState const &UR, Integer n) {
		FpeGuard fpeGuard{};
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
		τ = pow(ei, inv(Γ));
		setEnergy();
	}
	void setVelocity(Vector<Type, dimCount> const &v) {
		std::copy(v.begin(), v.end(), m.begin());
		m *= ρ;
		setEnergy();
	}
	Type internalEnergy(Type const &ek) const {
		using std::max;
		auto const ei = max(e - ek, 0_R);
		return select(ei >= 1e-3_R * e, ei, pow(τ, Γ));
	}
	Type internalEnergy() const {
		return internalEnergy(0.5_R * m.dot(m) * inv(ρ));
	}
	GasState updateEntropy() const {
		GasState U = *this;
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const ei = e - 0.5_R * ρ * u.dot(u);
		U.τ = select(ei < 1e-1_R * e, τ, pow(ei, inv(Γ)));
		return U;
	}
	GasState &setEnergy() {
		auto const iρ = inv(ρ);
		auto const u = iρ * m;
		auto const ek = 0.5_R * ρ * u.dot(u);
		auto const ei = pow(τ, Γ);
		e = ei + ek;
		return *this;
	}
	void setTemperature(Type T) {
		constexpr auto constant = kB / mᵤ;
		setPressure(constant * ρ * T * inv(μ));
		setEnergy();
	}
	Type getTemperature() const {
		constexpr auto constant = (Γ - 1_R) * mᵤ / kB;
		return constant * μ * internalEnergy() / ρ;
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

private:
	Type &ρ;
	Vector<Type, dimCount> &m;
	Type &e;
	Type &τ;
};

#endif /* GasState_HPP_ */
