/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "na/zoned/scheduler/MinFlowScheduler.hpp"

#include "ir/Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/zoned/Architecture.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace na::zoned {
MinFlowScheduler::FlowNetwork::FlowNetwork(const size_t reserveNumEdges) {
  edgeTarget_.reserve(reserveNumEdges);
  reverseEdgeTarget_.reserve(reserveNumEdges);
  edgeCapacity_.reserve(reserveNumEdges);
  edgeUnitCost_.reserve(reserveNumEdges);
}
auto MinFlowScheduler::FlowNetwork::validateVertexIndex(
    const VertexIndex i) const -> void {
  if (i >= getNumVertices()) {
    std::ostringstream ss;
    ss << "Vertex index " << i << " is out of range [0, " << getNumVertices()
       << ").";
    throw std::out_of_range(ss.str());
  }
}
auto MinFlowScheduler::FlowNetwork::validateEdgeIndex(const EdgeIndex i) const
    -> void {
  if (i >= 2 * getNumEdges()) {
    std::ostringstream ss;
    ss << "Edge index " << i << " is out of range [0, 2 * " << getNumEdges()
       << ").";
    throw std::out_of_range(ss.str());
  }
}
auto MinFlowScheduler::FlowNetwork::ensureBuilt() const -> void {
  if (!isBuilt) {
    throw std::logic_error("The flow network has not been built yet.");
  }
}
auto MinFlowScheduler::FlowNetwork::ensureNotBuilt() const -> void {
  if (isBuilt) {
    throw std::logic_error("The flow network has already been built.");
  }
}
auto MinFlowScheduler::FlowNetwork::toFlowQuantityWithOverflowCheck(
    const CapacityValue capacity) -> FlowQuantity {
  constexpr CapacityValue maxCapacity =
      std::numeric_limits<FlowQuantity>::max();
  if (capacity > maxCapacity) {
    std::ostringstream ss;
    ss << "Capacity " << capacity << " exceeds maximum capacity "
       << maxCapacity;
    throw std::invalid_argument(ss.str());
  }
  return static_cast<FlowQuantity>(capacity);
}
auto MinFlowScheduler::FlowNetwork::getForwardOutDegree(
    const VertexIndex v) const -> EdgeIndex {
  ensureBuilt();
  validateVertexIndex(v);
  return vertexFirstOutgoingForwardEdge_[v + 1] -
         vertexFirstOutgoingForwardEdge_[v];
}
auto MinFlowScheduler::FlowNetwork::getOutgoingForwardEdges(
    const VertexIndex v) const -> std::ranges::iota_view<EdgeIndex, EdgeIndex> {
  ensureBuilt();
  validateVertexIndex(v);
  return std::views::iota(vertexFirstOutgoingForwardEdge_[v],
                          vertexFirstOutgoingForwardEdge_[v + 1]);
}
auto MinFlowScheduler::FlowNetwork::getForwardSuccessors(
    const VertexIndex v) const -> std::span<const VertexIndex> {
  ensureBuilt();
  validateVertexIndex(v);
  using diff_t = std::iter_difference_t<decltype(edgeTarget_.begin())>;
  const auto start_it = edgeTarget_.cbegin() +
                        static_cast<diff_t>(vertexFirstOutgoingForwardEdge_[v]);
  const auto end_it =
      edgeTarget_.cbegin() +
      static_cast<diff_t>(vertexFirstOutgoingForwardEdge_[v + 1]);
  return {start_it, end_it};
}
auto MinFlowScheduler::FlowNetwork::isBackwardEdge(const EdgeIndex i) const
    -> bool {
  ensureBuilt();
  validateEdgeIndex(i);
  return i >= getNumEdges();
}
auto MinFlowScheduler::FlowNetwork::getReverseEdge(const EdgeIndex i) const
    -> EdgeIndex {
  ensureBuilt();
  validateEdgeIndex(i);
  return reverseEdge_[i];
}
auto MinFlowScheduler::FlowNetwork::getSource(const EdgeIndex i) const
    -> VertexIndex {
  ensureBuilt();
  validateEdgeIndex(i);
  return getTarget(getReverseEdge(i));
}
auto MinFlowScheduler::FlowNetwork::getTarget(const EdgeIndex i) const
    -> VertexIndex {
  ensureBuilt();
  validateEdgeIndex(i);
  if (isBackwardEdge(i)) {
    return reverseEdgeTarget_[i - getNumEdges()];
  }
  return edgeTarget_[i];
}
auto MinFlowScheduler::FlowNetwork::getFlow(EdgeIndex e) const -> FlowQuantity {
  ensureBuilt();
  validateEdgeIndex(e);
  return isBackwardEdge(e) ? -edgeFlow_[getReverseEdge(e)] : edgeFlow_[e];
}
auto MinFlowScheduler::FlowNetwork::residualEdgeCapacity(const EdgeIndex e)
    -> FlowQuantity {
  if (isBackwardEdge(e)) {
    return residualEdgeCapacityFast(getReverseEdge(e), true);
  }
  return residualEdgeCapacityFast(e, false);
}
auto MinFlowScheduler::FlowNetwork::residualEdgeCapacityFast(
    const EdgeIndex forwardEdge, const bool forBackwardEdge) -> FlowQuantity {
  return forBackwardEdge ? edgeFlow_[forwardEdge]
                         : edgeCapacity_[forwardEdge] - edgeFlow_[forwardEdge];
}
auto MinFlowScheduler::FlowNetwork::initializePreflow(const VertexIndex source)
    -> void {
  vertexExcess_.assign(getNumVertices(), 0);
  vertexPotential_.assign(getNumVertices(), 0);
  vertexPotential_[source] = static_cast<CostValue>(getNumVertices());
  edgeFlow_.assign(getNumEdges(), 0);
  activeNodes_ = {};
  std::ranges::for_each(getOutgoingForwardEdges(source),
                        [&](const EdgeIndex e) -> void {
                          const auto c = edgeCapacity_[e];
                          edgeFlow_[e] = c;
                          const auto v = getTarget(e);
                          vertexExcess_[v] = c;
                          vertexExcess_[source] -= c;
                          activeNodes_.push(v);
                        });
}
auto MinFlowScheduler::FlowNetwork::push(const EdgeIndex e) -> void {
  const auto u = getSource(e);
  const auto v = getTarget(e);
  const bool backwardEdge = isBackwardEdge(e);
  const auto forwardEdge = backwardEdge ? getReverseEdge(e) : e;
  const auto delta = std::min(
      vertexExcess_[u], residualEdgeCapacityFast(forwardEdge, backwardEdge));
  assert(delta > 0);
  if (backwardEdge) {
    edgeFlow_[forwardEdge] -= delta;
  } else {
    edgeFlow_[forwardEdge] += delta;
  }
  vertexExcess_[u] -= delta;
  if (vertexExcess_[v] == 0) {
    activeNodes_.push(v);
  }
  vertexExcess_[v] += delta;
}

auto MinFlowScheduler::FlowNetwork::relabel(const VertexIndex u) -> void {
  auto minPotential = std::numeric_limits<CostValue>::max();
  for (const auto e : getAllOutgoingEdges(u)) {
    if (residualEdgeCapacity(e) > 0) {
      const auto v = getTarget(e);
      minPotential = std::min(minPotential, vertexPotential_[v]);
    }
  }
  assert(minPotential < std::numeric_limits<CostValue>::max());
  vertexPotential_[u] = minPotential + 1;
}
auto MinFlowScheduler::FlowNetwork::discharge(const VertexIndex u) -> void {
  assert(vertexExcess_[u] > 0);
  do {
    for (const auto e : getAllOutgoingEdges(u)) {
      if (const auto v = getTarget(e);
          residualEdgeCapacity(e) > 0 &&
          vertexPotential_[u] == vertexPotential_[v] + 1) {
        push(e);
        if (vertexExcess_[u] == 0) {
          return;
        }
      }
    }
    assert(vertexExcess_[u] > 0);
    relabel(u);
  } while (vertexExcess_[u] > 0);
}
auto MinFlowScheduler::FlowNetwork::solveMaxFlow(const VertexIndex source,
                                                 const VertexIndex sink)
    -> void {
  ensureBuilt();
  validateVertexIndex(source);
  validateVertexIndex(sink);
  if (source == sink) {
    throw std::invalid_argument("Source and sink cannot be the same.");
  }
  maximumFlow_ = 0;
  initializePreflow(source);
  while (!activeNodes_.empty()) {
    const auto u = activeNodes_.front();
    activeNodes_.pop();
    if (u != source && u != sink) {
      discharge(u);
    }
  }
  assert(-vertexExcess_[source] == vertexExcess_[sink]);
  maximumFlow_ = vertexExcess_[sink];
}
auto MinFlowScheduler::FlowNetwork::refine() -> void {
  for (VertexIndex u = 0; u < getNumVertices(); ++u) {
    const auto p = vertexPotential_[u];
    for (const auto e : getAllOutgoingEdges(u)) {
      const auto v = getTarget(e);
      // rc = reduced cost
      if (const auto rc = (isBackwardEdge(e) ? -edgeUnitCost_[getReverseEdge(e)]
                                             : edgeUnitCost_[e]) +
                          p - vertexPotential_[v];
          rc < 0) {
        const auto delta = edgeCapacity_[e] - edgeFlow_[e];
        edgeFlow_[e] = edgeCapacity_[e];
        vertexExcess_[u] -= delta;
        vertexExcess_[v] += delta;
      }
    }
  }
  activeNodes_ = {};
  for (VertexIndex u = 0; u < getNumVertices(); ++u) {
    if (vertexExcess_[u] > 0) {
      activeNodes_.push(u);
    }
  }
  while (!activeNodes_.empty()) {
    const auto u = activeNodes_.front();
    activeNodes_.pop();
    discharge2(u);
  }
}
auto MinFlowScheduler::FlowNetwork::discharge2(const VertexIndex u) -> void {
  assert(vertexExcess_[u] > 0);
  do {
    for (const auto e : getAllOutgoingEdges(u)) {
      const auto v = getTarget(e);
      if (residualEdgeCapacity(e) > 0) {
        // rc = reduced cost
        if (const auto rc =
                (isBackwardEdge(e) ? -edgeUnitCost_[getReverseEdge(e)]
                                   : edgeUnitCost_[e]) +
                vertexPotential_[u] - vertexPotential_[v];
            rc < 0) {
          push(e);
          if (vertexExcess_[u] == 0) {
            return;
          }
        }
      }
    }
    assert(vertexExcess_[u] > 0);
    relabel2(u);
  } while (vertexExcess_[u] > 0);
}
auto MinFlowScheduler::FlowNetwork::relabel2(const VertexIndex u) -> void {
  auto maxPotential = std::numeric_limits<CostValue>::min();
  for (const auto e : getAllOutgoingEdges(u)) {
    if (residualEdgeCapacity(e) > 0) {
      const auto v = getTarget(e);
      maxPotential =
          std::max(maxPotential,
                   vertexPotential_[v] -
                       (isBackwardEdge(e) ? -edgeUnitCost_[getReverseEdge(e)]
                                          : edgeUnitCost_[e]) -
                       epsilon_);
    }
  }
  vertexPotential_[u] = maxPotential;
}
auto MinFlowScheduler::FlowNetwork::solveMinCostMaxFlow(
    const VertexIndex source, const VertexIndex sink) -> void {
  solveMaxFlow(source, sink);
  vertexPotential_.assign(getNumVertices(), 0);
  vertexExcess_[source] = 0;
  vertexExcess_[sink] = 0;
  scaleCosts();
  epsilon_ = maxEdgeCost_ * static_cast<CostValue>(getNumEdges());
  do {
    epsilon_ = std::max(CostValue{1}, epsilon_ / 5); // magic number
    refine();
  } while (epsilon_ > 1);
  unscaleCosts();
}
auto MinFlowScheduler::FlowNetwork::scaleCosts() -> void {
  std::ranges::for_each(edgeUnitCost_,
                        [this](auto& cost) { cost *= getNumVertices(); });
}
auto MinFlowScheduler::FlowNetwork::unscaleCosts() -> void {
  std::ranges::for_each(edgeUnitCost_,
                        [this](auto& cost) { cost /= getNumVertices(); });
}
auto MinFlowScheduler::FlowNetwork::addVertex() -> VertexIndex {
  ensureNotBuilt();
  return numVertices_++;
}
auto MinFlowScheduler::FlowNetwork::addEdgeWithCapacityAndUnitCost(
    const VertexIndex source, const VertexIndex target,
    const CapacityValue capacity, const CostValue unitCost) -> EdgeIndex {
  ensureNotBuilt();
  validateVertexIndex(source);
  validateVertexIndex(target);

  reverseEdgeTarget_.emplace_back(source);
  edgeTarget_.emplace_back(target);
  edgeCapacity_.emplace_back(toFlowQuantityWithOverflowCheck(capacity));
  edgeUnitCost_.emplace_back(unitCost);
  maxEdgeCost_ = std::max(maxEdgeCost_, unitCost);

  return getNumEdges() - 1;
}
auto MinFlowScheduler::FlowNetwork::build(
    IVector<EdgeIndex, EdgeIndex>& permutation) -> void {
  ensureNotBuilt();
  permutation.clear();
  isBuilt = true;
  // if the graph has zero vertices, it also has zero edges, and there is
  // nothing to build; hence, we return early
  if (hasZeroVertices()) {
    // currently does not free unused memory in case the capacity of the
    // vectors is > 0
    return;
  }
  // since the following vectors are considered fixed after the end of this
  // function call, we shrink the size to fit its content to release
  // unused memory
  vertexExcess_.shrink_to_fit();
  edgeTarget_.shrink_to_fit();
  edgeCapacity_.shrink_to_fit();
  edgeUnitCost_.shrink_to_fit();
  reverseEdgeTarget_.shrink_to_fit();

  // Count the edges outgoing from each vertex to initialize the
  // vertexFirstEdge_ vector, which will contain for each vertex the index
  // of the first outgoing edge in the edgeTarget_ vector.
  vertexFirstOutgoingForwardEdge_.assign(getNumVertices() + 1, 0);
  vertexFirstOutgoingForwardEdge_.shrink_to_fit();
  std::ranges::for_each(reverseEdgeTarget_, [&](const VertexIndex i) -> void {
    ++vertexFirstOutgoingForwardEdge_[i];
  });
  // Go through the counts and compute the exclusive prefix sum.
  std::exclusive_scan(vertexFirstOutgoingForwardEdge_.begin(),
                      vertexFirstOutgoingForwardEdge_.end(),
                      vertexFirstOutgoingForwardEdge_.begin(), EdgeIndex{0});
  // Check the sentinel value
  assert(vertexFirstOutgoingForwardEdge_[getNumVertices()] == getNumEdges() &&
         "Sentinel value of vertexFirstOutgoingEdge_ is not correct.");
  // Calculate a permutation of edges to sort them by their source vertex.
  // Note this temporarily alters the vertexFirstOutgoingEdge_ vector, whose
  // state is restored afterward.
  permutation.reserve(getNumEdges());
  std::ranges::for_each(reverseEdgeTarget_, [&](const VertexIndex i) {
    permutation.emplace_back(vertexFirstOutgoingForwardEdge_[i]++);
  });
  // Restore the previous state of the vertexFirstOutgoingEdge_ vector.
  for (VertexIndex i = getNumVertices() - 1; i > 0; --i) {
    vertexFirstOutgoingForwardEdge_[i] = vertexFirstOutgoingForwardEdge_[i - 1];
  }
  vertexFirstOutgoingForwardEdge_[0] = 0;
  // Apply the permutation to the target vertices of the edges.
  applyPermutation(permutation, reverseEdgeTarget_);
  applyPermutation(permutation, edgeTarget_);
  applyPermutation(permutation, edgeCapacity_);
  applyPermutation(permutation, edgeUnitCost_);

  // Count the reverse edges outgoing from each vertex to initialize the
  // vertexFirstReverseEdge_ vector, which will contain for each vertex the
  // index of the first outgoing reverse edge in the reverseEdgeTarget_
  // vector.
  vertexFirstOutgoingBackwardEdge_.assign(getNumVertices() + 1, 0);
  vertexFirstOutgoingBackwardEdge_.shrink_to_fit();
  std::ranges::for_each(edgeTarget_, [&](const VertexIndex i) -> void {
    ++vertexFirstOutgoingBackwardEdge_[i];
  });
  // Go through the counts and compute the exclusive prefix sum.
  std::exclusive_scan(vertexFirstOutgoingBackwardEdge_.begin(),
                      vertexFirstOutgoingBackwardEdge_.end(),
                      vertexFirstOutgoingBackwardEdge_.begin(), EdgeIndex{0});
  // Check the sentinel value
  assert(vertexFirstOutgoingBackwardEdge_[getNumVertices()] == getNumEdges() &&
         "Sentinel value of vertexFirstOutgoingReverseEdge_ is not correct.");
  // Calculate a permutation of reverse edges to sort them by their source
  // vertex. Note this temporarily alters the
  // vertexFirstOutgoingReverseEdge_ vector, whose state is restored
  // afterward.
  IVector<EdgeIndex, EdgeIndex> reversePermutation;
  reversePermutation.reserve(getNumEdges());
  std::ranges::for_each(edgeTarget_, [&](const VertexIndex i) -> void {
    reversePermutation.emplace_back(vertexFirstOutgoingBackwardEdge_[i]++);
  });
  // Restore the previous state of the vertexFirstOutgoingReverseEdge_
  // vector.
  for (VertexIndex i = getNumVertices() - 1; i > 0; --i) {
    vertexFirstOutgoingBackwardEdge_[i] =
        vertexFirstOutgoingBackwardEdge_[i - 1];
  }
  vertexFirstOutgoingBackwardEdge_[0] = 0;
  // Apply the permutation to the target vertices of the reverse edges.
  applyPermutation(reversePermutation, reverseEdgeTarget_);
  // Initialize the vector to translate between edges and their reverse
  // edges
  reverseEdge_.assign(2 * getNumEdges(), EdgeIndex{0});
  for (EdgeIndex i = 0; i < getNumEdges(); ++i) {
    const auto edge = permutation[i];
    const auto reverseEdge = reversePermutation[edge] + getNumEdges();
    reverseEdge_[edge] = reverseEdge;
    reverseEdge_[reverseEdge] = edge;
  }
}
MinFlowScheduler::MinFlowScheduler(const Architecture& architecture,
                                   const Config& config)
    : architecture_(architecture), config_(config) {
  // Validate maxFillingFactor
  if (config_.maxFillingFactor < 0.0 || config_.maxFillingFactor > 1.0) {
    std::ostringstream oss;
    oss << "Invalid maxFillingFactor: " << config_.maxFillingFactor
        << ". Value must be in the range [0.0, 1.0].";
    throw std::invalid_argument(oss.str());
  }
  // calculate the maximum possible number of two-qubit gates per layer
  for (const auto& zone : architecture_.get().entanglementZones) {
    maxTwoQubitGateNumPerLayer_ +=
        std::max(static_cast<size_t>(1U),
                 static_cast<size_t>(config_.maxFillingFactor *
                                     static_cast<double>(zone->front().nRows *
                                                         zone->front().nCols)));
  }
  if (maxTwoQubitGateNumPerLayer_ == 0) {
    throw std::invalid_argument("Architecture must contain at least one site "
                                "in an entanglement zone");
  }
}
auto MinFlowScheduler::schedule(const qc::QuantumComputation& qc) const
    -> std::pair<std::vector<SingleQubitGateLayer>,
                 std::vector<TwoQubitGateLayer>> {
  if (qc.empty()) {
    // early exit if there are no operations to schedule
    return std::pair{std::vector<SingleQubitGateLayer>{},
                     std::vector<TwoQubitGateLayer>{}};
  }
  std::vector<SingleQubitGateLayer> singleQubitGateLayers(1);
  std::vector<TwoQubitGateLayer> twoQubitGateLayers(0);
  // the following vector contains a mapping from qubits to the layer where
  // the next two-qubit gate can be scheduled for that qubit, i.e., the layer
  // after the last layer with a two-qubit gate acting on that qubit
  std::vector<size_t> nextLayerForQubit(qc.getNqubits(), 0);
  for (const auto& op : qc) {
    if (op->getType() == qc::Barrier) {
      if (op->getNqubits() < qc.getNqubits()) {
        throw std::invalid_argument("Only global barriers are allowed.");
      }
      // set the next layer for all qubits to the currently last layer
      assert(twoQubitGateLayers.size() + 1 == singleQubitGateLayers.size());
      const auto newNextLayerForQubit = twoQubitGateLayers.size();
      for (qc::Qubit q = 0; q < qc.getNqubits(); ++q) {
        nextLayerForQubit[q] = newNextLayerForQubit;
      }
    } else if (op->isGlobal(qc.getNqubits()) && !op->isControlled() &&
               qc.getNqubits() > 1) {
      const auto maxNextLayerForQubit =
          *std::ranges::max_element(std::as_const(nextLayerForQubit));
      for (qc::Qubit q = 0; q < qc.getNqubits(); ++q) {
        nextLayerForQubit[q] = maxNextLayerForQubit;
      }
      singleQubitGateLayers[maxNextLayerForQubit].emplace_back(*op);
    } else if (op->isStandardOperation()) {
      const auto& stdOp = dynamic_cast<qc::StandardOperation&>(*op);
      if (stdOp.getNtargets() == 1 && stdOp.getNcontrols() == 0) {
        singleQubitGateLayers[nextLayerForQubit[stdOp.getTargets().front()]]
            .emplace_back(stdOp);
      } else if (stdOp.getType() == qc::Z && stdOp.getNtargets() == 1 &&
                 stdOp.getNcontrols() == 1) {
        const auto qubit1 = stdOp.getTargets().front();
        const auto qubit2 = stdOp.getControls().cbegin()->qubit;
        auto layer =
            std::max(nextLayerForQubit[qubit1], nextLayerForQubit[qubit2]);
        while (layer < twoQubitGateLayers.size() &&
               twoQubitGateLayers[layer].size() >=
                   maxTwoQubitGateNumPerLayer_) {
          ++layer;
        }
        assert(layer <= twoQubitGateLayers.size());
        if (layer == twoQubitGateLayers.size()) {
          // add a new layer
          singleQubitGateLayers.emplace_back();
          twoQubitGateLayers.emplace_back();
        }
        twoQubitGateLayers[layer].emplace_back(
            std::array{std::min(qubit1, qubit2), std::max(qubit1, qubit2)});
        nextLayerForQubit[qubit1] = layer + 1;
        nextLayerForQubit[qubit2] = layer + 1;
      } else {
        std::stringstream ss;
        ss << "Operation type not supported: " << stdOp.getType() << " with "
           << stdOp.getNcontrols() << " controls and " << stdOp.getNtargets()
           << " targets";
        throw std::invalid_argument(ss.str());
      }
    } else {
      std::stringstream ss;
      ss << "Operation type not supported: " << op->getType() << " with "
         << op->getNcontrols() << " controls and " << op->getNtargets()
         << " targets";
      throw std::invalid_argument(ss.str());
    }
  }
  return std::pair{singleQubitGateLayers, twoQubitGateLayers};
}
} // namespace na::zoned
