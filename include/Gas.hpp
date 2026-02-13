/*
 * Gas.hpp
 *
 *  Created on: Feb 9, 2026
 *      Author: dmarce1
 */

#ifndef GAS_HPP_
#define GAS_HPP_

#include <config.hpp>

#include <array>
#include <span>

#include "Vector.hpp"

namespace Gas {
using namespace Config;

template <int N>
struct State {
	static constexpr int fieldCount() {
		return dimCount + 3;
	}
	static constexpr int size() {
		return N;
	}
	State() {
		size_t index = 0;
		for (int d = 0; d < dimCount; d++) {
			momentumDensity[d] = std::span<Real>(state_.begin() + index, N);
			index += N;
		}
		massDensity = std::span<Real>(state_.begin() + index, N);
		index += N;
		energyDensity = std::span<Real>(state_.begin() + index, N);
		index += N;
		entropyTracer = std::span<Real>(state_.begin() + index, N);
		index += N;
	}
	std::array<Real, N> signalSpeed(int direction) const {
		std::array<Real, N> signalSpeed = Real(0);
		for (int n = 0; n < N; n++) {
			Real const invMassDensity = Real(1) / massDensity[n];
			Real const halfInvMassDensity = Real(0.5) * invMassDensity;
			Real kineticEnergyDensity{};
			Real velocity = invMassDensity * momentumDensity[direction][n];
			for (int d = 0; d < dimCount; d++) {
				kineticEnergyDensity += halfInvMassDensity * momentumDensity[d][n];
			}
			Real internalEnergyDensity = energyDensity[n] - kineticEnergyDensity;
			internalEnergyDensity = (internalEnergyDensity > dualEngeryPressureSwitch * energyDensity[n])
										? internalEnergyDensity
										: pow(entropyTracer[n], fluidGamma);
			Real const pressure = (fluidGamma - 1) * internalEnergyDensity;
			Real const soundSpeed = sqrt(fluidGamma * pressure * invMassDensity);
			signalSpeed = std::max(std::abs(velocity + soundSpeed), std::abs(velocity - soundSpeed));
		}
		return signalSpeed;
	}
	auto flux(int direction) const {
		State flux;
		for (int n = 0; n < N; n++) {
			Real const invMassDensity = Real(1) / massDensity[n];
			Real const halfInvMassDensity = Real(0.5) * invMassDensity;
			Real kineticEnergyDensity{};
			Real velocity = invMassDensity * momentumDensity[direction][n];
			for (int d = 0; d < dimCount; d++) {
				kineticEnergyDensity += halfInvMassDensity * momentumDensity[d][n];
			}
			Real internalEnergyDensity = energyDensity[n] - kineticEnergyDensity;
			internalEnergyDensity = (internalEnergyDensity > dualEngeryPressureSwitch * energyDensity[n])
										? internalEnergyDensity
										: pow(entropyTracer[n], fluidGamma);
			Real const pressure = (fluidGamma - 1) * internalEnergyDensity;
			flux.massDensity[n] = momentumDensity[direction][n];
			for (int d = 0; d < dimCount; d++) {
				flux.momentumDensity[d][n] = velocity * momentumDensity[d][n];
			}
			flux.momentumDensity[direction][n] += pressure;
			flux.energyDensity[n] = velocity * (energyDensity[n] + pressure);
			flux.entropyTracer[n] = velocity * entropyTracer[n];
		}
		return flux;
	}
	friend State solveRiemannProblem(State const &ul, State const &ur, int direction) {
		State flux;
		auto const leftSignalSpeed = ul.signalSpeed(direction);
		auto const rightSignalSpeed = ur.signalSpeed(direction);
		auto const leftFlux = ul.flux(direction);
		auto const rightFlux = ur.flux(direction);
		int i = 0;
		for (int n = 0; n < N; n++) {
			auto const maximumSignalSpeed = std::max(leftSignalSpeed[n], rightSignalSpeed[n]);
			for (int f = 0; f < fieldCount(); f++) {
				flux.state_[i] = Real(0.5) * (leftFlux.state_[i] + rightFlux.state_[i]);
				flux.state_[i] += Real(0.5) * maximumSignalSpeed * (ul.state_[i] - ur.state_[i]);
				i++;
			}
		}
		return flux;
	}

private:
	std::array<Real, fieldCount() * size()> state_;
	std::array<std::span<Real>, dimCount> momentumDensity;
	std::span<Real> massDensity;
	std::span<Real> energyDensity;
	std::span<Real> entropyTracer;
};

} // namespace Gas

#endif /* GAS_HPP_ */
