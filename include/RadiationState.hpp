/*
 * RadiationState.hpp
 *
 *  Created on: Feb 25, 2026
 *      Author: dmarce1
 */

#ifndef INCLUDE_RADIATIONSTATE_HPP_
#define INCLUDE_RADIATIONSTATE_HPP_

/*AαΑBβΒGγΓDδΔEεΕZζΖHηΗThθΘIιΙKκΚLλλMμΜNνΝXξΞOοΟPπΠRρΡSσΣTτΤYυΥFφΦChχΧPsψΨOωΩ*/
/*₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵨ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁱⁿᵃᵇᶜᵈᵉᶠᵍʰⁱʲᵏˡᵐⁿᵒᵖʳˢᵗᵘᵛʷˣʸᶻᵅᵝᵞᵟᵋᵠᵡ*/

#include <string>
#include <tuple>

#include "AutoDiff.hpp"
#include "Constants.hpp"
#include "Debug.hpp"
#include "GasState.hpp"
#include "Integer.hpp"
#include "Matrix.hpp"
#include "Real.hpp"
#include "Vector.hpp"

#include <numeric>

template <typename Type, Integer dimCount>
struct RadiationState : Vector<Type, 1 + dimCount> {
private:
	static constexpr Integer eidx = 0;
	static constexpr Integer fidx = 1;
	static constexpr auto zero = 0.0_R;
	static constexpr auto half = 0.5_R;
	static constexpr auto one = 1.0_R;
	static constexpr auto σ = codeUnits().σ;

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
		E_(*(new(&(*this)[0]) Type{})), F_(*(new(&(*this)[1]) Vector<Type, dimCount>{})) {
	}
	RadiationState(RadiationState const &other) :
		E_(*(new(&(*this)[0]) Type{})), F_(*(new(&(*this)[1]) Vector<Type, dimCount>{})) {
		*this = other;
	}
	RadiationState(base_type const &other) :
		E_(*(new(&(*this)[0]) Type{})), F_(*(new(&(*this)[1]) Vector<Type, dimCount>{})) {
		*this = other;
	}
	virtual ~RadiationState() {
		E_.~Type();
		F_.~Vector<Type, dimCount>();
	}
	RadiationState &operator=(RadiationState const &other) {
		base_type::operator=(other);
		return *this;
	}
	RadiationState &operator=(base_type const &other) {
		base_type::operator=(other);
		return *this;
	}
	template <Integer count>
	Type enforcePositivity(Vector<RadiationState, count> const &u) const {
		constexpr auto δ = 1e-3_R;
		using std::abs;
		using std::min;
		using Auto = AutoDiff<Type, 1>;
		Type θ = 1_R;
		for (Integer j = 0; j < count; j++) {
			RadiationState const du = u[j] - *this;
			{
				auto const δE = (δ - 1_R) * E_;
				auto const θn = δE * inv(select(du.E_ < 0_R, du.E_, δE));
				θ = min(θ, θn);
			}
			{
				auto const residual = [du, this](auto θn) {
					return sqr(E_ + θn * du.E_) - F_.dot(F_) - 2_R * θn * F_.dot(du.F_) - sqr(θn) * du.F_.dot(du.F_);
				};
				if (any(residual(θ) < 0_R)) {
					Type err = 1_R;
					Auto θn = Auto::genVar(θ);
					for (Integer i = 0; any(err > 4_R * eps_R); i++) {
						assert(i < 40);
						Auto const f = residual(θn);
						Type const dθn = -f.get() * inv(f.get(0));
						θn = Auto::genVar(θn.get() + dθn);
						err = abs(f.get() * sqr(inv(E_)));
						std::cout << err << std::endl;
					}
					θ = min(θ, θn.get());
				}
			}
		}
		if (any(θ < 1_R)) {
			for (Integer j = 0; j < count; j++) {
				RadiationState const du = u[j] - *this;
				u[j] -= (1_R - θ) * (u[j] - *this);
			}
		}
	}
	auto eigenstructure(Integer n) const {
		FpeGuard fpeGuard{};
		ASSERT_POSITIVE(E_);
		using std::abs;
		constexpr Integer am = 0;
		constexpr Integer ap = dimCount;
		auto const F2 = F_.dot(F_);
		auto const E2 = sqr(E_);
		ASSERT_NONNEGATIVE(E2 - F2);
		auto const H = (1_R / 3_R) * (2_R * E_ + sqrt(4_R * E2 - 3_R * F2));
		auto const β = F_ * inv(H);
		auto const β2 = β.dot(β);
		auto const &βz = β[n];
		auto const βz2 = sqr(βz);
		auto const Ξ = sqrt((1_R - β2) * (3_R - β2 - 2_R * βz2));
		Matrix<Type, size()> A, B;
		Vector<Type, size()> λ;
		auto const lambda = [β, βz, n](auto const &λ) {
			Vector<Type, size()> r{};
			auto const λ2 = sqr(λ);
			r[0] = 1_R - 2_R * βz * λ + λ2;
			for (Integer k = 0; k < dimCount; k++) {
				r[k + 1] = β[k] * (1_R - λ2);
			}
			r[n + 1] += 2_R * (λ - βz);
			return r;
		};
		λ.front() = (2_R * βz - Ξ) * inv(3_R - β2);
		λ.back() = (2_R * βz + Ξ) * inv(3_R - β2);
		for (Integer k = 1; k < dimCount; k++) {
			λ[k] = βz;
		}
		B.setCol(am, lambda(λ.front()));
		B.setCol(ap, lambda(λ.back()));
		if constexpr (dimCount > 1) {
			B[0][1] = β2 - βz2;
			for (Integer k = 0; k < dimCount; k++) {
				B[k + 1][1] = select(β[k] != 0_R, β[k] * (1_R - βz2), 1_R);
			}
			B[n + 1][1] -= select(β[n] != 0_R, β[n] * (1_R - β2), 1_R);
			if constexpr (dimCount > 2) {
				auto const s1 = (n == 0) ? 1 : 0;
				auto const s2 = (n == 2) ? 1 : 2;
				auto const nz = sqr(β[s1]) + sqr(β[s2]) != 0_R;
				B[s1 + 1][2] = select(nz, -β[s2], -1_R);
				B[s2 + 1][2] = select(nz, +β[s1], +1_R);
			}
		}
		A[0][0] = 1.5_R;
		for (Integer j = 0, k = 1; j < dimCount; j++, k++) {
			A[0][k] = 0.5_R * β[j];
			A[k][0] = β[j];
			A[k][k] = 1_R;
		}
		auto const R = A * B;
		auto const L = inv(R);
		return std::tuple(
			λ,
			[R](RadiationState const &dc) {
				return R * dc;
			},
			[L](RadiationState const &du) {
				return L * du;
			});
	}
	auto flux(Integer k) const {
		RadiationState f;
		auto const [H, β] = primitives();
		return flux(H, β, k);
	}
	friend auto riemannFlux(RadiationState const &uₗ, RadiationState const &uᵣ, Integer n) {
		FpeGuard fpeGuard{};
		auto const [Hᵣ, βᵣ] = uᵣ.primitives();
		auto const [Hₗ, βₗ] = uₗ.primitives();
		auto const β2ᵣ = βᵣ.dot(βᵣ);
		auto const β2ₗ = βₗ.dot(βₗ);
		auto const fₗ = uₗ.flux(Hₗ, βₗ, n);
		auto const fᵣ = uᵣ.flux(Hᵣ, βᵣ, n);
		auto const Ξᵣ = sqrt((1_R - β2ᵣ) * (3_R - β2ᵣ - 2_R * sqr(βᵣ[n])));
		auto const Ξₗ = sqrt((1_R - β2ₗ) * (3_R - β2ₗ - 2_R * sqr(βₗ[n])));
		auto const λpᵣ = (2_R * βᵣ[n] + Ξᵣ) / (3_R - β2ᵣ);
		auto const λpₗ = (2_R * βₗ[n] + Ξₗ) / (3_R - β2ₗ);
		auto const λmᵣ = (2_R * βᵣ[n] - Ξᵣ) / (3_R - β2ᵣ);
		auto const λmₗ = (2_R * βₗ[n] - Ξₗ) / (3_R - β2ₗ);
		auto const λᵣ = max(0_R, max(λpᵣ, λpₗ));
		auto const λₗ = min(0_R, min(λmᵣ, λmₗ));
		auto const f = (λᵣ * fₗ - λₗ * fᵣ + λᵣ * λₗ * (uᵣ - uₗ)) / (λᵣ - λₗ);
		return std::pair<RadiationState, Type>(f, max(λᵣ, -λₗ));
	}
	auto jacobian(Integer n) const {
		using Auto = AutoDiff<Type, dimCount + 1>;
		using AutoVector = Vector<Auto, dimCount>;
		AutoVector F;
		Auto const E = Auto::genVar(E_, 0);
		for (Integer k = 0; k < dimCount; k++) {
			F[k] = Auto::genVar(F_[k], k + 1);
		}
		Auto const F2 = F.dot(F);
		Auto const E2 = sqr(E);
		Auto const H = (1_R / 3_R) * (2_R * E + sqrt(4_R * E2 - 3_R * F2));
		auto const β = F * inv(H);
		Auto fE = F[n];
		AutoVector fF = β[n] * H * β;
		fF[n] += 0.25_R * (1_R - β.dot(β)) * H;
		Matrix<Type, dimCount + 1> J;
		for (Integer j = 0; j <= dimCount; j++) {
			J[0][j] = fE.get(j);
		}
		for (Integer j = 0; j < dimCount; j++) {
			for (Integer k = 0; k <= dimCount; k++) {
				J[j + 1][k] = fF[j].get(k);
			}
		}
		return J;
	}
	void setTemperature(Type T) {
		setEnergy(4_R * σ * sqr(sqr(T)));
	}
	Type getTemperature() const {
		return root(E_ * inv(4_R * σ), 4_I);
	}
	void setEnergy(Type e) {
		E_ = e;
	}
	void setFlux(Vector<Type, dimCount> const &f) {
		F_ = f;
	}
	template <typename, Integer>
	friend struct CoupledState;
	friend std::ostream &operator<<(std::ostream &os, RadiationState const &u) {
		os << "(E=" << u.E_ << ", ";
		for (Integer i = 0; i < dimCount; i++) {
			os << "F_" << std::string(1, 'x' + i) << u.F_[i] << ", ";
		}
		os << ")";
		return os;
	}
	auto flux(Type const &H, Vector<Type, dimCount> const &β, Integer k) const {
		RadiationState f;
		f.E_ = F_[k];
		f.F_ = β[k] * H * β;
		f.F_[k] += 0.25_R * (1_R - β.dot(β)) * H;
		return f;
	}

private:
	auto primitives() const {
		Type const F2 = F_.dot(F_);
		Type const E2 = sqr(E_);
		ASSERT_NONNEGATIVE(E2 - F2);
		Type const H = (1_R / 3_R) * (2_R * E_ + sqrt(4_R * E2 - 3_R * F2));
		Vector<Type, dimCount> const β = F_ * inv(H);
		return std::pair<Type, Vector<Type, dimCount>>(H, β);
	}
	Type &E_;
	Vector<Type, dimCount> &F_;
};

#endif /* INCLUDE_RADIATIONSTATE_HPP_ */
