/*
 * CoupledState.hpp
 *
 *  Created on: Feb 26, 2026
 *      Author: dmarce1
 */
// harold joseph 2253376858
#ifndef INCLUDE_COUPLEDSTATE_HPP_
#define INCLUDE_COUPLEDSTATE_HPP_

#include "AutoDiff.hpp"
#include "Constants.hpp"
#include "GasState.hpp"
#include "Material.hpp"
#include "RadiationState.hpp"

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλΛMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

template <typename Type, Integer dimCount>
struct CoupledState {
	using Gas = GasState<Type, dimCount>;
	using Rad = RadiationState<Type, dimCount>;
	static constexpr auto Γ = Gas::Γ;
	static constexpr auto units = codeUnits();
	static constexpr auto kB = units.kB;
	static constexpr auto mᵤ = units.mᵤ;
	static constexpr auto a = 4_R * units.σ;
	static constexpr Integer size() {
		return Gas::size() + Rad::size();
	}
	static auto getFieldNames() {
		static auto const mnames = Gas::getFieldNames();
		static auto const rnames = Rad::getFieldNames();
		std::array<std::string, size()> names;
		std::copy_n(mnames.begin(), Gas::size(), names.begin());
		std::copy_n(rnames.begin(), Rad::size(), names.begin() + Gas::size());
		return names;
	}
	CoupledState() = default;
	CoupledState(CoupledState const &other) = default;
	CoupledState(Vector<Type, size()> const &other) {
		*this = other;
	}
	virtual ~CoupledState() = default;
	CoupledState &operator=(CoupledState const &other) = default;
	CoupledState &operator=(Vector<Type, size()> const &other) {
		M_ = other.template sub<0, Gas::size()>();
		R_ = other.template sub<Gas::size(), size()>();
		return *this;
	}
	//	operator Vector<Type, size()>() const {
	//		return concatenate(M_, R_);
	//	}
	//	operator Vector<Type, size()>&() const {
	//		return reinterpret_cast<Vector<Type, size()>&>(*this);
	//	}
	Type operator[](Integer n) const {
		return (n < Gas::size()) ? M_[n] : R_[n - Gas::size()];
	}
	Type &operator[](Integer n) {
		return (n < Gas::size()) ? M_[n] : R_[n - Gas::size()];
	}
	CoupledState &operator+=(CoupledState const &other) {
		*this = *this + other;
		return *this;
	}
	CoupledState &operator-=(CoupledState const &other) {
		*this = *this - other;
		return *this;
	}
	CoupledState operator+(CoupledState const &other) const {
		CoupledState result;
		result.M_ = M_ + other.M_;
		result.R_ = R_ + other.R_;
		return result;
	}
	CoupledState operator-(CoupledState const &other) const {
		CoupledState result;
		result.M_ = M_ - other.M_;
		result.R_ = R_ - other.R_;
		return result;
	}
	friend CoupledState operator*(Type const &scalar, CoupledState const &state) {
		CoupledState result;
		result.M_ = scalar * state.M_;
		result.R_ = scalar * state.R_;
		return result;
	}
	friend CoupledState operator*(CoupledState const &state, Type const &scalar) {
		return scalar * state;
	}
	CoupledState updateEntropy() const {
		CoupledState u;
		u.matter() = M_.updateEntropy();
		u.radiation() = R_;
		return u;
	};
	auto flux(Integer n) const {
		return CoupledState(concatenate(M_.flux(n), R_.flux(n)));
	};
	auto eigenstructure(Integer n) const {
		auto const [λm, Rm, Lm] = M_.eigenstructure(n);
		auto const [λr, Rr, Lr] = R_.eigenstructure(n);
		return std::tuple(
			concatenate(λm, λr),
			[Rm, Rr](CoupledState const &dc) {
				CoupledState du;
				du.matter() = Rm(dc.M_);
				du.radiation() = Rr(dc.R_);
				return du;
			},
			[Lm, Lr](CoupledState const &du) {
				CoupledState dc;
				dc.matter() = Lm(du.M_);
				dc.radiation() = Lr(du.R_);
				return dc;
			});
	}
	friend auto riemannFlux(CoupledState const &ul, CoupledState const &ur, Integer n) {
		using std::max;
		auto const [Fm, λm] = riemannFlux(ul.M_, ur.M_, n);
		auto const [Fr, λr] = riemannFlux(ul.R_, ur.R_, n);
		CoupledState F;
		F.matter() = Fm;
		F.radiation() = Fr;
		return std::pair(F, max(λm, λr));
	}
	auto solveImplicit(Real μ, Real κₐ, Real κₛ, Real dt) const {
		FpeGuard fpeGuard{};
		using Vector4 = Vector<Type, dimCount + 1>;
		auto const ρ = M_.ρ;
		auto const κ = ρ * κₐ;
		auto const χ = ρ * (κₐ + κₛ);
		auto const Ur0 = concatenate(R_.E_, R_.F_);
		auto const e0 = pow(M_.τ, Γ);
		auto const s0 = M_.m;
		auto const b = a * pow((Γ - 1_R) * μ * mᵤ * inv(kB * ρ), 4_R);
		auto const iρ = inv(ρ);
		auto const computeResidualAndJacobian = [b, dt, κ, χ, iρ, ρ, Ur0, s0, e0](auto const &dx0) {
			auto dx = cvtVector2AutoVector(dx0);
			auto const Ur = Ur0 + dx;
			auto const ds = -dx.template sub<1, 1 + dimCount>();
			auto const E = Ur[0];
			auto const F = Ur.template sub<1, dimCount + 1>();
			auto const s = s0 + ds;
			auto const e = e0 - (dx[0] + iρ * ds.dot(s));
			ASSERT_NONNEGATIVE(E.get());
			auto const F2 = F.dot(F);
			auto const iF2 = inv(F2 + eps_R);
			auto const n = F * sqrt(iF2);
			auto const f2 = F2 * inv(sqr(E));
			auto const Δ2 = 4_R - 3_R * f2;
			ASSERT_NONNEGATIVE(e.get());
			ASSERT_NONNEGATIVE(Δ2.get());
			auto const ξ = (3_R + 4_R * f2) * inv(5_R + 2_R * sqrt(Δ2));
			auto const Dd = 0.5_R - 0.5_R * ξ;
			auto const Ds = 1.5_R * ξ - 0.5_R;
			auto const β = iρ * s;
			auto const B = b * sqr(sqr(e));
			auto const g0 = κ * (E - 2_R * β.dot(F) - B);
			auto const G0 = χ * (F - E * β * (1_R + Dd) - Ds * E * β.dot(n) * n);
			auto const g = g0 + β.dot(G0);
			auto const G = G0 + β * g0;
			auto const R = concatenate(g, G);
			auto const h = Ur - Ur0 + dt * R;
			Vector<Real, dimCount + 1> h_;
			Matrix<Real, dimCount + 1> dhdx;
			for (Integer n = 0; n <= dimCount; n++) {
				h_[n] = h[n].get();
				for (Integer m = 0; m <= dimCount; m++) {
					dhdx[n][m] = h[n].get(m);
				}
			}
			return std::tuple(h_, dhdx);
		};
		constexpr Real tol = 1e-20_R;
		Real err;
		Vector4 U = Ur0;
		U[0] += M_.e;
		Vector4 x(0_R);
		do {
			auto const [f, dfdx] = computeResidualAndJacobian(x);
			//			std::cout << f << std::endl;
			err = sqrt(f.dot(f)) / abs(U);
			x -= f * inv(dfdx);
			//			std::cout << "err = " << err << std::endl;
			//			std::cout << std::endl;
		} while (err > tol);
	}
	void setTemperature(Type T) {
		M_.setTemperature(T);
		R_.setTemperature(T);
	}
	Gas &matter() {
		return M_;
	}
	Rad &radiation() {
		return R_;
	}
	Gas const &matter() const {
		return M_;
	}
	Rad const &radiation() const {
		return R_;
	}
	template <Integer count>
	auto enforcePositivity(Vector<CoupledState, count> &u) const {
		Vector<Gas, count> um;
		Vector<Rad, count> ur;
		for (Integer j = 0; j < count; j++) {
			um[j] = u[j].matter();
			ur[j] = u[j].radiation();
		}
		return std::pair(matter().enforcePositivity(um), radiation().enforcePositivity(ur));
	}

private:
	Gas M_;
	Rad R_;
};

#endif /* INCLUDE_COUPLEDSTATE_HPP_ */
