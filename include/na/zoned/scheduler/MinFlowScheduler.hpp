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
#include <queue>
#include <stack>
#include <utility>
#include <vector>

namespace na::zoned {
// A wrapper to add value_type to a view that lacks it.
template <std::ranges::view R>
struct value_type_wrapper
    : public std::ranges::view_interface<value_type_wrapper<R>> {
  R r;

  // Explicitly define value_type
  using value_type = std::iter_value_t<std::ranges::iterator_t<R>>;

  value_type_wrapper() = default;
  constexpr explicit value_type_wrapper(R r) : r(std::move(r)) {}

  constexpr auto begin() { return std::ranges::begin(r); }
  constexpr auto end() { return std::ranges::end(r); }
  constexpr auto begin() const { return std::ranges::begin(r); }
  constexpr auto end() const { return std::ranges::end(r); }
};

// Deduction guide
template <typename R>
value_type_wrapper(R&&) -> value_type_wrapper<std::views::all_t<R>>;

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
  class FlowNetwork {
  public:
    using VertexIndex = size_t;
    using EdgeIndex = size_t;
    using FlowQuantity = int64_t;
    using CapacityValue = uint64_t;
    using CostValue = int64_t;

    /// Wrapper around a vector that allows to specify the index type.
    template <std::unsigned_integral I, typename T>
    class IVector : std::vector<T> {
    public:
      using typename std::vector<T>::value_type;
      using std::vector<T>::vector;
      using std::vector<T>::empty;
      using std::vector<T>::emplace_back;
      using std::vector<T>::clear;
      using std::vector<T>::shrink_to_fit;
      using std::vector<T>::assign;
      using std::vector<T>::begin;
      using std::vector<T>::cbegin;
      using std::vector<T>::end;
      using std::vector<T>::cend;
      using std::vector<T>::resize;
      [[nodiscard]] auto operator[](I index) const -> const T& {
        return std::vector<T>::operator[](
            static_cast<std::vector<T>::size_type>(index));
      }
      [[nodiscard]] auto operator[](I index) -> T& {
        return std::vector<T>::operator[](
            static_cast<std::vector<T>::size_type>(index));
      }
      auto reserve(I index) -> void {
        return std::vector<T>::reserve(
            static_cast<std::vector<T>::size_type>(index));
      }
      [[nodiscard]] constexpr auto size() const -> I {
        return I(std::vector<T>::size());
      }
    };

  private:
    //===------------------------------------------------------------------===//
    // Vertex data
    //===------------------------------------------------------------------===//
    /**
     * Used to store the excess of vertices of a pre-flow during the computation
     * of the actual flow.
     */
    IVector<VertexIndex, FlowQuantity> vertexExcess_;
    /**
     * During the calculation of the max flow, the potential is a measure of
     * distance to the sink. The sink always has potential 0, and the source's
     * potential is equal to the number of vertices.
     */
    IVector<VertexIndex, CostValue> vertexPotential_;
    /**
     * For fast and cache-efficient access of all successors. It is a vector
     * whose values denote the index in the edgeTarget_ vector where the
     * successors of this vertex start.
     */
    IVector<VertexIndex, EdgeIndex> vertexFirstOutgoingForwardEdge_;
    /**
     * Every forward edge in the graph has a corresponding reverse backward
     * edge. The following data structure is analogue to the one above.
     */
    IVector<VertexIndex, EdgeIndex> vertexFirstOutgoingBackwardEdge_;

    //===------------------------------------------------------------------===//
    // Edge data
    //===------------------------------------------------------------------===//
    /**
     * Stores the target vertex of each forward edge. More precisely, the target
     * vertex of the i-th forward edge is stored as the i-th element in the
     * vector.
     */
    IVector<EdgeIndex, VertexIndex> edgeTarget_;
    /**
     * Stores the target vertex of each backward edge. More precisely, the
     * target vertex of the i-th backward edge is stored as the i-th element in
     * the vector. Note that the actual indices of backward edges are shifted by
     * the number of forward edges, i.e., they start with m, where m is the
     * number of forward edges. Also note that this vector is used to retrieve
     * the source vertex of forward edges.
     */
    IVector<EdgeIndex, VertexIndex> reverseEdgeTarget_;
    /**
     * The capacity of each forward edge denotes the maximum possible flow along
     * this forward edge. The capacity of the i-th edge is stored as the i-th
     * element in this vector. Backward edges do not have a capacity.
     */
    IVector<EdgeIndex, FlowQuantity> edgeCapacity_;
    /**
     * The unit cost is the cost per unit of flow along this forward edge. To
     * retrieve the actual cost, the unit cost is multiplied with the flow along
     * this forward edge. The unit cost of the i-th edge is the i-th element in
     * the vector. Backward edges do not have a cost.
     */
    IVector<EdgeIndex, CostValue> edgeUnitCost_;
    /**
     * Stores the finally calculated flow in the flow network. During the run of
     * the push-relabel algorithm, it contains the pre-flow.
     */
    IVector<EdgeIndex, FlowQuantity> edgeFlow_;
    /**
     * For each edge, it maps to the index of the reverse edge. Note that the
     * indices of backward edges start with m, where m is the number of forward
     * edges.
     */
    IVector<EdgeIndex, EdgeIndex> reverseEdge_;
    /**
     * During the calculation of the max flow, this queue stores all active
     * nodes, i.e., nodes with positive excess.
     */
    std::queue<VertexIndex> activeNodes_;

    //===------------------------------------------------------------------===//
    // Global variables
    //===------------------------------------------------------------------===//
    /// Total number of vertices
    VertexIndex numVertices_ = 0;
    /// Whether the graph is already finalized and the flow can be calculated.
    bool isBuilt = false;
    /// The maximum cost of a forward edge
    CostValue maxEdgeCost_ = 0;
    /// The resulting flow of the min-cost flow.
    FlowQuantity maximumFlow_ = 0;
    // todo(yannick): add docstring
    CostValue epsilon_ = 0;

  public:
    FlowNetwork() = default;
    FlowNetwork(size_t reserveNumEdges);
    [[nodiscard]] constexpr auto hasZeroVertices() const -> bool {
      return numVertices_ == 0;
    }
    [[nodiscard]] constexpr auto getNumVertices() const -> VertexIndex {
      return numVertices_;
    }
    [[nodiscard]] constexpr auto getNumEdges() const -> EdgeIndex {
      return edgeTarget_.size();
    }
    auto validateVertexIndex(VertexIndex i) const -> void;
    auto validateEdgeIndex(EdgeIndex i) const -> void;
    auto ensureBuilt() const -> void;
    auto ensureNotBuilt() const -> void;
    static auto toFlowQuantityWithOverflowCheck(CapacityValue capacity)
        -> FlowQuantity;
    [[nodiscard]] auto getForwardOutDegree(VertexIndex v) const -> EdgeIndex;
    [[nodiscard]] auto getOutgoingForwardEdges(VertexIndex v) const
        -> std::ranges::iota_view<EdgeIndex, EdgeIndex>;
    [[nodiscard]] auto getAllOutgoingEdges(const VertexIndex v) const -> auto {
      ensureBuilt();
      validateVertexIndex(v);
      const auto forwardView = getOutgoingForwardEdges(v);
      const auto backwardView = std::views::iota(
          getNumEdges() + vertexFirstOutgoingBackwardEdge_[v],
          getNumEdges() + vertexFirstOutgoingBackwardEdge_[v + 1]);
      // todo(yannick): replace with a more efficient concatenation view when
      //  available in the standard library or write a custom one
      std::vector buffer(forwardView.begin(), forwardView.end());
      buffer.insert(buffer.end(), backwardView.begin(), backwardView.end());
      return buffer;
    }
    [[nodiscard]] auto getForwardSuccessors(VertexIndex v) const
        -> std::span<const VertexIndex>;
    [[nodiscard]] auto isBackwardEdge(EdgeIndex i) const -> bool;
    [[nodiscard]] auto getReverseEdge(EdgeIndex i) const -> EdgeIndex;
    [[nodiscard]] auto getSource(EdgeIndex i) const -> VertexIndex;
    [[nodiscard]] auto getTarget(EdgeIndex i) const -> VertexIndex;
    [[nodiscard]] auto getFlow(EdgeIndex e) const -> FlowQuantity;
    auto addVertex() -> VertexIndex;
    auto addEdgeWithCapacityAndUnitCost(VertexIndex source, VertexIndex target,
                                        CapacityValue capacity,
                                        CostValue unitCost) -> EdgeIndex;
    /**
     * @brief Finalizes the graph structure for efficient access afterward.
     * @note This function must be called before calling @ref solve. Vice versa,
     * @ref addVertexWithSupply and @ref addEdgeWithCapacityAndUnitCost must be
     * called before calling this function.
     * @details This function reorders the edges such that they are sorted by
     * their source vertex. Furthermore, it initializes a vector that carries
     * the index to the first target vertex for each source vertex. This makes
     * iteration over outgoing edges more efficient. In particular, it makes
     * access to edges and associated data, e.g., capacity more cache efficient.
     */
    auto build(IVector<EdgeIndex, EdgeIndex>& permutation) -> void;
    /**
     * Applies a given permutation to a container. The size of the permutation
     * may be less than or equal to the size of the container. In the former
     * case, only the first elements up to the size of the permutation are
     * permuted.
     * @tparam T is the value type of the @ref IVector.
     * @param permutation is the permutation to apply to the data elements.
     * @param data is the container to be permuted.
     */
    template <typename T>
    static auto
    applyPermutation(const IVector<EdgeIndex, EdgeIndex>& permutation,
                     IVector<EdgeIndex, T>& data) -> void {
      IVector<EdgeIndex, T> data_copy = data;
      for (EdgeIndex i = 0; i < permutation.size(); ++i) {
        assert(permutation[i] < data.size());
        data[permutation[i]] = data_copy[i];
      }
    }
    [[nodiscard]] auto residualEdgeCapacity(EdgeIndex e) -> FlowQuantity;
    [[nodiscard]] auto residualEdgeCapacityFast(EdgeIndex forwardEdge,
                                                bool forBackwardEdge = false)
        -> FlowQuantity;
    auto initializePreflow(VertexIndex source) -> void;
    auto push(EdgeIndex e) -> void;
    auto relabel(VertexIndex u) -> void;
    auto discharge(VertexIndex u) -> void;
    auto solveMaxFlow(VertexIndex source, VertexIndex sink) -> void;
    auto scaleCosts() -> void;
    auto unscaleCosts() -> void;
    auto refine() -> void;
    auto discharge2(VertexIndex u) -> void;
    auto push2(EdgeIndex u) -> void;
    auto relabel2(VertexIndex u) -> void;
    auto solveMinCostMaxFlow(VertexIndex source, VertexIndex sink) -> void;
    [[nodiscard]] auto getMaximumFlow() const -> uint64_t {
      return static_cast<uint64_t>(maximumFlow_);
    }
  };

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
   * @param qc is the quantum computation to be scheduled
   * @return a pair of two vectors. The first vector contains the layers of
   * single-qubit operations. The second vector contains the layers of two-qubit
   * operations. A pair of qubits represents every two-qubit operation.
   */
  [[nodiscard]] auto schedule(const qc::QuantumComputation& qc) const
      -> std::pair<std::vector<SingleQubitGateLayer>,
                   std::vector<TwoQubitGateLayer>>;
};
} // namespace na::zoned
