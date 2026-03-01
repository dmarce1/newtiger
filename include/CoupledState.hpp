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
	static constexpr auto c = units.c;
	static constexpr auto kB = units.kB;
	static constexpr auto mᵤ = units.mᵤ;
	static constexpr auto σ = units.σ;
	static constexpr auto ic = inv(c);
	static constexpr auto a = 4_R * σ * ic;
	CoupledState(Gas &gas, Radiation &radiation) :
		gas_{gas}, radiation_{radiation} {
	}
	auto solveImplicit(Real μ, Real κₐ, Real κₛ, Real dt) const {
		FpeGuard fpeGuard{};
		using Vector4 = Vector<Type, dimCount + 1>;
		//		P =  (Γ - 1_R) * ei = kB / (μ * mᵤ) * ρ * T;
		auto const iCv = ((Γ - 1_R) * μ * mᵤ) / kB;
		auto const κ = κₐ + κₛ;
		auto const ρ = gas_.ρ;
		auto const iρ = inv(ρ);
		auto const Ur0 = concatenate(radiation_.E, radiation_.F);
		auto const Ug0 = concatenate(gas_.e, c * gas_.m);
		auto const computeResidualAndJacobian = [iCv, dt, κₐ, iρ, ρ, κ, Ur0, Ug0](auto const &dx0) {
			auto dx = cvtVector2AutoVector(dx0);
			auto const Ur = Ur0 + dx;
			auto const Ug = Ug0 - dx;
			auto const e = Ug[0];
			auto const E = Ur[0];
			auto const s = ic * Ug.template sub<1, dimCount + 1>();
			auto const F = Ur.template sub<1, dimCount + 1>();
			auto const iE = inv(E);
			ASSERT_NONNEGATIVE(E.get());
			ASSERT_NONNEGATIVE(e.get());
			auto const F2 = F.dot(F);
			auto const n = F * sqrt(inv(F2.get() != 0_R ? F2 : 1_R));
			auto const f2 = F2 * sqr(iE);
			ASSERT_RANGE(0_R, f2.get(), 1_R);
			auto const χ = (3_R + 4_R * f2) / (5_R + 2_R * sqrt(4_R - 3_R * f2));
			auto const Ds = 1.5_R * χ - 0.5_R;
			auto const Dd = 0.5_R - 0.5_R * χ;
			auto const u = iρ * s;
			auto const ei = max(e - 0.5_R * iρ * u.dot(u), 0_R);
			auto const β = ic * u;
			auto const T = iCv * iρ * ei;
			auto const T4 = sqr(sqr(T));
			auto const g0 = κₐ * (E - 2_R * β.dot(F) - a * T4);
			auto const G0 = κ * (F - E * ((1_R + Dd) * β + Ds * β.dot(n) * n));
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
			return std::tuple(h_, dhdx, sqrt(sqrt(E.get() / a)), T, sqrt(f2.get()) * n, β);
		};
		constexpr Real tol = 1e-9_R;
		Real err;
		Vector4 U = Ur0 + Ug0;
		//		Vector4 Ur = Ur0;
		//		Vector4 Ug = Ug0;
		Vector4 x(0_R);
		do {
			auto const [f, dfdx, Tr, Tg, βr, βg] = computeResidualAndJacobian(x);
			std::cout << f << std::endl;
			err = sqrt(f.dot(f)) / abs(U);
			x -= f * inv(dfdx);
			//			Ur = Ur0 + cvtAutoVector2Vector(x);
			//			Ug = Ug0 - cvtAutoVector2Vector(x);
			std::cout << "Tr = " << Tr << " βr = " << βr << std::endl;
			std::cout << "Tg = " << Tg << " βg = " << βg << std::endl;
			std::cout << "err = " << err << std::endl;
			std::cout << std::endl;
		} while (err > tol);
	}

private:
	Gas &gas_;
	Radiation &radiation_;
};

#endif /* INCLUDE_COUPLEDSTATE_HPP_ */
