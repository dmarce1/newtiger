#pragma once

#include "Vector.hpp"

namespace config {
using Real = double;

constexpr int dimCount = 1;

constexpr Real fluidGamma = Real(5) / Real(3);
constexpr Real dualEnergySwitch1 = Real(1e-3);
constexpr Real dualEnergySwitch2 = Real(1e-1);
constexpr Real meanMolecularWeight = Real(4) / Real(3);

} // namespace config
