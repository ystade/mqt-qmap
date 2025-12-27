/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#pragma once

#include "ir/QuantumComputation.hpp"
#include "na/zoned/Architecture.hpp"
#include "na/zoned/Types.hpp"
#include "na/zoned/scheduler/SchedulerBase.hpp"

#include <functional>
#include <utility>
#include <vector>

namespace na::zoned {
/**
 * The class MinFlowScheduler implements the min-cost flow scheduling
 * strategy for the zoned neutral atom compiler.
 */
class MinFlowScheduler : public SchedulerBase {
  /// A reference to the zoned neutral atom architecture
  std::reference_wrapper<const Architecture> architecture_;
  /**
   * This value is calculated based on the architecture and indicates the
   * entanglement zone.
   */
  size_t maxTwoQubitGateNumPerLayer_ = 0;

public:
  /// The configuration of the MinFlowScheduler
  struct Config {
    /// The maximal share of traps that are used in the entanglement zone.
    double maxFillingFactor = 0.9;
    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(Config, maxFillingFactor);
  };

private:
  /// The configuration of the MinFlowScheduler
  Config config_;

public:
  /**
   * Create a new MinFlowScheduler.
   * @param architecture is the architecture of the neutral atom system
   * @param config is the configuration for the scheduler
   */
  MinFlowScheduler(const Architecture& architecture, const Config& config);
  /**
   * This function schedules the operations of a quantum computation.
   * @details todo: docstring
   * @param qc is the quantum computation
   * @return a pair of two vectors. The first vector contains the layers of
   * single-qubit operations. The second vector contains the layers of two-qubit
   * operations. A pair of qubits represents every two-qubit operation.
   */
  [[nodiscard]] auto schedule(const qc::QuantumComputation& qc) const
      -> std::pair<std::vector<SingleQubitGateLayer>,
                   std::vector<TwoQubitGateLayer>>;
};
} // namespace na::zoned
