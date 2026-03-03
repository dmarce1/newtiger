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
	using Radiation = RadiationState<Type, dimCount>;
	static constexpr auto Γ = Gas::Γ;
	static constexpr auto units = codeUnits();
	static constexpr auto kB = units.kB;
	static constexpr auto mᵤ = units.mᵤ;
	static constexpr auto σ = units.σ;
	static constexpr auto a = 4_R * σ;
	CoupledState(Gas &gas, Radiation &radiation) :
		gas_{gas}, radiation_{radiation} {
	}
	auto solveImplicit(Real μ, Real κ, Real κₛ, Real dt) const {
		FpeGuard fpeGuard{};
		using Vector4 = Vector<Type, dimCount + 1>;
		auto const χ = κ + κₛ;
		auto const ρ = gas_.ρ;
		auto const Ur0 = concatenate(radiation_.E, radiation_.F);
		auto const e0 = pow(gas_.τ, Γ);
		auto const s0 = gas_.m;
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
			auto const d = 4_R * sqr(E) - 3_R * F.dot(F);
			ASSERT_NONNEGATIVE(E.get());
			ASSERT_NONNEGATIVE(e.get());
			ASSERT_NONNEGATIVE(d.get());
			auto const H = (1_R / 3_R) * (2_R * E + sqrt(d));
			auto const βr = F * inv(H);
			auto const βg = iρ * s;
			auto const B = b * sqr(sqr(e));
			auto const g0 = κ * H * (0.75_R + (0.25_R * βr - 2_R * βg).dot(βr)) - κ * B;
			auto const G0 = χ * H * ((1_R - βg.dot(βr)) * βr - βg);
			auto const g = g0 + βg.dot(G0);
			auto const G = G0 + βg * g0;
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
		U[0] += gas_.e;
		Vector4 x(0_R);
		do {
			auto const [f, dfdx] = computeResidualAndJacobian(x);
			std::cout << f << std::endl;
			err = sqrt(f.dot(f)) / abs(U);
			x -= f * inv(dfdx);
			std::cout << "err = " << err << std::endl;
			std::cout << std::endl;
		} while (err > tol);
	}

private:
	Gas &gas_;
	Radiation &radiation_;
};

#endif /* INCLUDE_COUPLEDSTATE_HPP_ */
