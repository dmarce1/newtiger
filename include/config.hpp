#pragma once

#include "Vector.hpp"

namespace Config {

using Real = double;

constexpr int interiorWidth = 256;
constexpr int ghostWidth = 1;
constexpr int dimCount = 1;
constexpr int exteriorWidth = interiorWidth + 2 * ghostWidth;

constexpr Real dualEngeryPressureSwitch = 0.001;
constexpr Real dualEngerUpdateSwitch = 0.1;
constexpr Real fluidGamma = 5.0 / 3.0;

} // namespace config
