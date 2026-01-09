/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#include "ir/operations/StandardOperation.hpp"
#include "na/zoned/scheduler/MinFlowScheduler.hpp"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>
#include <utility>
#include <vector>

namespace testing {
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
MATCHER_P(RefEq, value, "") { return arg.get() == value; }
} // namespace testing
namespace na::zoned {
constexpr std::string_view architectureJson = R"({
  "name": "min_flow_scheduler_architecture",
  "storage_zones": [{
    "zone_id": 0,
    "slms": [{"id": 0, "site_separation": [3, 3], "r": 20, "c": 20, "location": [0, 0]}],
    "offset": [0, 0],
    "dimension": [60, 60]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {"id": 1, "site_separation": [12, 10], "r": 4, "c": 4, "location": [5, 70]},
      {"id": 2, "site_separation": [12, 10], "r": 4, "c": 4, "location": [7, 70]}
    ],
    "offset": [5, 70],
    "dimension": [50, 40]
  }],
  "aods":[{"id": 0, "site_separation": 2, "r": 20, "c": 20}],
  "rydberg_range": [[[5, 70], [55, 110]]]
})";
class MinFlowSchedulerScheduleTest : public ::testing::Test {
protected:
  Architecture architecture;
  MinFlowScheduler::Config config{.maxFillingFactor = .8};
  MinFlowScheduler scheduler;
  MinFlowSchedulerScheduleTest()
      : architecture(Architecture::fromJSONString(architectureJson)),
        scheduler(architecture, config) {}
};
TEST_F(MinFlowSchedulerScheduleTest, NoGate) {
  qc::QuantumComputation qc;
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);
  EXPECT_THAT(singleQubitGateLayers, ::testing::IsEmpty());
  EXPECT_THAT(twoQubitGateLayers, ::testing::IsEmpty());
}
TEST_F(MinFlowSchedulerScheduleTest, SingleQubitGate) {
  //    ┌───────┐
  // q: ┤ Rz(π) ├
  //    └───────┘
  qc::QuantumComputation qc(1);
  qc.rz(qc::PI, 0);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);
  EXPECT_THAT(singleQubitGateLayers,
              ::testing::ElementsAre(::testing::ElementsAre(::testing::RefEq(
                  static_cast<qc::StandardOperation&>(*qc.at(0))))));
  EXPECT_THAT(twoQubitGateLayers, ::testing::IsEmpty());
}
TEST_F(MinFlowSchedulerScheduleTest, TwoQubitGate) {
  // q_0: ─■─
  //       │
  // q_1: ─■─
  qc::QuantumComputation qc(2);
  qc.cz(0, 1);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);
  EXPECT_THAT(
      singleQubitGateLayers,
      ::testing::ElementsAre(::testing::IsEmpty(), ::testing::IsEmpty()));
  EXPECT_THAT(twoQubitGateLayers,
              ::testing::UnorderedElementsAre(::testing::ElementsAre(
                  ::testing::UnorderedElementsAre(0U, 1U))));
}
TEST_F(MinFlowSchedulerScheduleTest, SingleQubitSandwich) {
  // q_0: ──────────■──────────
  //      ┌───────┐ │ ┌───────┐
  // q_1: ┤ Rz(π) ├─■─┤ Rz(π) ├
  //      └───────┘   └───────┘
  qc::QuantumComputation qc(2);
  qc.rz(qc::PI, 1);
  qc.cz(0, 1);
  qc.rz(qc::PI, 1);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);
  EXPECT_THAT(singleQubitGateLayers,
              ::testing::ElementsAre(
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(0)))),
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(2))))));
  EXPECT_THAT(twoQubitGateLayers,
              ::testing::ElementsAre(::testing::UnorderedElementsAre(
                  ::testing::UnorderedElementsAre(0U, 1U))));
}
TEST_F(MinFlowSchedulerScheduleTest, TwoQubitSequence) {
  // q_0: ─■───────
  //       │
  // q_1: ─■──■────
  //          │
  // q_2: ────■──■─
  //             │
  // q_3: ───────■─
  qc::QuantumComputation qc(4);
  qc.cz(0, 1);
  qc.cz(1, 2);
  qc.cz(2, 3);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);
  EXPECT_THAT(singleQubitGateLayers, ::testing::SizeIs(4));
  EXPECT_THAT(singleQubitGateLayers, ::testing::Each(::testing::IsEmpty()));
  EXPECT_THAT(
      twoQubitGateLayers,
      ::testing::ElementsAre(::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(0U, 1U)),
                             ::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(1U, 2U)),
                             ::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(2U, 3U))));
}
TEST_F(MinFlowSchedulerScheduleTest, Mixed) {
  //            INPUT ORDER                         SCHEDULED ORDER
  // q_0: ─■─────────────────────────  >>>  ─────────░─■─░─────────░───░─
  //       │ ┌───────┐                 >>>           ░ │ ░┌───────┐░   ░
  // q_1: ─■─┤ Rz(π) ├─────────────■─  >>>  ─────────░─■─░┤ Rz(π) ├░─■─░─
  //         └───────┘┌───────┐    │   >>>  ┌───────┐░   ░└───────┘░ │ ░
  // q_2: ────────────┤ Rz(π) ├─■──■─  >>>  ┤ Rz(π) ├░─■─░─────────░─■─░─
  //                  └───────┘ │      >>>  └───────┘░ │ ░         ░   ░
  // q_3: ──────────────────────■────  >>>  ─────────░─■─░─────────░───░─
  qc::QuantumComputation qc(4);
  qc.cz(0, 1);
  qc.rz(qc::PI, 1);
  qc.rz(qc::PI, 2);
  qc.cz(2, 3);
  qc.cz(1, 2);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);
  EXPECT_THAT(singleQubitGateLayers,
              ::testing::ElementsAre(
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(2)))),
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(1)))),
                  ::testing::IsEmpty()));
  EXPECT_THAT(
      twoQubitGateLayers,
      ::testing::ElementsAre(::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(0U, 1U),
                                 ::testing::UnorderedElementsAre(2U, 3U)),
                             ::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(1U, 2U))));
}
TEST_F(MinFlowSchedulerScheduleTest, Barrier) {
  // q_0: ─■─────────░───
  //       │┌───────┐░
  // q_1: ─■┤ Rz(π) ├░───
  //        └───────┘░
  // q_2: ───────────░─■─
  //                 ░ │
  // q_3: ───────────░─■─
  qc::QuantumComputation qc(4);
  qc.cz(0, 1);
  qc.rz(qc::PI, 1);
  qc.barrier();
  qc.cz(2, 3);
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);
  EXPECT_THAT(singleQubitGateLayers,
              ::testing::ElementsAre(
                  ::testing::IsEmpty(),
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(1)))),
                  ::testing::IsEmpty()));
  EXPECT_THAT(
      twoQubitGateLayers,
      ::testing::ElementsAre(::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(0U, 1U)),
                             ::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(2U, 3U))));
}
TEST_F(MinFlowSchedulerScheduleTest, NonGlobalBarrier) {
  // q_0: ─░─
  //
  // q_1: ───
  qc::QuantumComputation qc(2);
  qc.emplace_back<qc::StandardOperation>(0, qc::Barrier);
  EXPECT_THROW(std::ignore = scheduler.schedule(qc), std::invalid_argument);
}
TEST_F(MinFlowSchedulerScheduleTest, NonGlobalCompound) {
  qc::QuantumComputation qc(2);
  qc::CompoundOperation compoundOp;
  compoundOp.emplace_back<qc::StandardOperation>(0, qc::RY,
                                                 std::vector{qc::PI_2});
  qc.emplace_back<qc::CompoundOperation>(compoundOp);
  EXPECT_THROW(std::ignore = scheduler.schedule(qc), std::invalid_argument);
}
TEST_F(MinFlowSchedulerScheduleTest, UnsupportedCXGate) {
  qc::QuantumComputation qc(2);
  qc.cx(0, 1);
  EXPECT_THROW(std::ignore = scheduler.schedule(qc), std::invalid_argument);
}
TEST_F(MinFlowSchedulerScheduleTest, FullEntanglementZone) {
  qc::QuantumComputation qc(26);
  for (qc::Qubit i = 0; i < 13; ++i) {
    qc.cz(2 * i, 2 * i + 1);
  }
  const auto& [singleQubitGateLayers, twoQubitGateLayers] =
      scheduler.schedule(qc);
  EXPECT_THAT(singleQubitGateLayers, ::testing::SizeIs(3));
  EXPECT_THAT(twoQubitGateLayers, ::testing::SizeIs(2));
}
TEST(MinFlowSchedulerConfigTest, InvalidMaxFillingFactor) {
  const auto architecture = Architecture::fromJSONString(architectureJson);
  constexpr MinFlowScheduler::Config config1{.maxFillingFactor = -0.1};
  EXPECT_THROW(MinFlowScheduler scheduler(architecture, config1),
               std::invalid_argument);
  constexpr MinFlowScheduler::Config config2{.maxFillingFactor = 1.1};
  EXPECT_THROW(MinFlowScheduler scheduler(architecture, config2),
               std::invalid_argument);
}
using Permutation = MinFlowScheduler::FlowNetwork::IVector<
    MinFlowScheduler::FlowNetwork::EdgeIndex,
    MinFlowScheduler::FlowNetwork::EdgeIndex>;
TEST(FlowNetwork, ApplyPermutation) {
  MinFlowScheduler::FlowNetwork::IVector<
      MinFlowScheduler::FlowNetwork::EdgeIndex, size_t>
      data = {10, 20, 30};
  MinFlowScheduler::FlowNetwork::applyPermutation<size_t>({2, 1, 0}, data);
  EXPECT_THAT(data, ::testing::ElementsAre(30, 20, 10));
  data = {10, 20, 30};
  MinFlowScheduler::FlowNetwork::applyPermutation<
      MinFlowScheduler::FlowNetwork::EdgeIndex>({1, 2, 0}, data);
  EXPECT_THAT(data, ::testing::ElementsAre(30, 10, 20));
  data = {10, 20, 30, 40};
  MinFlowScheduler::FlowNetwork::applyPermutation<
      MinFlowScheduler::FlowNetwork::EdgeIndex>({3, 2, 0, 1}, data);
  EXPECT_THAT(data, ::testing::ElementsAre(30, 40, 20, 10));
}
TEST(FlowNetwork, PreBuildExceptions) {
  MinFlowScheduler::FlowNetwork g;
  EXPECT_THROW(std::ignore = g.getTarget(0), std::logic_error);
  EXPECT_THROW(std::ignore = g.getSource(0), std::logic_error);
  EXPECT_THROW(std::ignore = g.getReverseEdge(0), std::logic_error);
  EXPECT_THROW(std::ignore = g.isBackwardEdge(0), std::logic_error);
  EXPECT_THROW(std::ignore = g.getForwardSuccessors(0), std::logic_error);
  EXPECT_THROW(std::ignore = g.getOutgoingForwardEdges(0), std::logic_error);
  EXPECT_THROW(std::ignore = g.getForwardOutDegree(0), std::logic_error);
  EXPECT_THROW(g.solveMaxFlow(0, 1), std::logic_error);
  EXPECT_THROW(g.solveMinCostMaxFlow(0, 1), std::logic_error);
}
TEST(FlowNetwork, PostBuildExceptions) {
  MinFlowScheduler::FlowNetwork g;
  Permutation permutation;
  g.build(permutation);
  EXPECT_THROW(std::ignore = g.getTarget(0), std::out_of_range);
  EXPECT_THROW(std::ignore = g.getSource(0), std::out_of_range);
  EXPECT_THROW(std::ignore = g.getReverseEdge(0), std::out_of_range);
  EXPECT_THROW(std::ignore = g.isBackwardEdge(0), std::out_of_range);
  EXPECT_THROW(std::ignore = g.getForwardSuccessors(0), std::out_of_range);
  EXPECT_THROW(std::ignore = g.getOutgoingForwardEdges(0), std::out_of_range);
  EXPECT_THROW(std::ignore = g.getForwardOutDegree(0), std::out_of_range);
}
TEST(FlowNetwork, Empty) {
  MinFlowScheduler::FlowNetwork g;
  Permutation permutation;
  EXPECT_NO_THROW(g.build(permutation));
  EXPECT_EQ(g.getNumVertices(), 0);
  EXPECT_EQ(g.getNumEdges(), 0);
}
TEST(FlowNetwork, OneVertex) {
  MinFlowScheduler::FlowNetwork g;
  g.addVertex();
  Permutation permutation;
  EXPECT_NO_THROW(g.build(permutation));
  EXPECT_EQ(g.getNumVertices(), 1);
  EXPECT_EQ(g.getNumEdges(), 0);
  EXPECT_EQ(g.getForwardOutDegree(0), 0);
  EXPECT_EQ(g.getForwardSuccessors(0).size(), 0);
}
TEST(FlowNetwork, OneEdge) {
  MinFlowScheduler::FlowNetwork g;
  const auto u = g.addVertex();
  const auto v = g.addVertex();
  const auto e = g.addEdgeWithCapacityAndUnitCost(u, v, 1, 1);
  Permutation permutation;
  EXPECT_NO_THROW(g.build(permutation));
  EXPECT_EQ(g.getNumVertices(), 2);
  EXPECT_EQ(g.getNumEdges(), 1);
  EXPECT_EQ(g.getForwardOutDegree(u), 1);
  EXPECT_EQ(g.getForwardOutDegree(v), 0);
  EXPECT_THAT(g.getForwardSuccessors(u), ::testing::ElementsAre(v));
  EXPECT_THAT(g.getForwardSuccessors(v), ::testing::IsEmpty());
  EXPECT_EQ(g.getSource(e), u);
  EXPECT_EQ(g.getTarget(e), v);
  const auto re = g.getReverseEdge(e);
  EXPECT_EQ(g.getReverseEdge(re), e);
  EXPECT_TRUE(g.isBackwardEdge(re));
  EXPECT_FALSE(g.isBackwardEdge(e));
  EXPECT_EQ(g.getSource(re), v);
  EXPECT_EQ(g.getTarget(re), u);
}
TEST(FlowNetwork, DiamondOrdered) {
  //                 ┌───┐
  // ┌───┐     ┌────>│ v ├─────┐     ┌───┐
  // │ u ├─────┤     ├───┤     ├────>│ x │
  // └───┘     └────>│ w ├─────┘     └───┘
  //                 └───┘
  MinFlowScheduler::FlowNetwork g;
  const auto u = g.addVertex();
  const auto v = g.addVertex();
  const auto w = g.addVertex();
  const auto x = g.addVertex();
  auto c = g.addEdgeWithCapacityAndUnitCost(u, v, 1, 1);
  auto d = g.addEdgeWithCapacityAndUnitCost(u, w, 1, 1);
  auto e = g.addEdgeWithCapacityAndUnitCost(v, x, 1, 1);
  auto f = g.addEdgeWithCapacityAndUnitCost(w, x, 1, 1);
  Permutation permutation;
  EXPECT_NO_THROW(g.build(permutation));
  c = permutation[c];
  d = permutation[d];
  e = permutation[e];
  f = permutation[f];
  EXPECT_EQ(g.getNumVertices(), 4);
  EXPECT_EQ(g.getNumEdges(), 4);
  EXPECT_EQ(g.getForwardOutDegree(u), 2);
  EXPECT_EQ(g.getForwardOutDegree(v), 1);
  EXPECT_EQ(g.getForwardOutDegree(w), 1);
  EXPECT_EQ(g.getForwardOutDegree(x), 0);
  EXPECT_THAT(g.getForwardSuccessors(u), ::testing::UnorderedElementsAre(v, w));
  EXPECT_THAT(g.getForwardSuccessors(v), ::testing::UnorderedElementsAre(x));
  EXPECT_THAT(g.getForwardSuccessors(w), ::testing::UnorderedElementsAre(x));
  EXPECT_THAT(g.getForwardSuccessors(x), ::testing::IsEmpty());
  EXPECT_EQ(g.getSource(c), u);
  EXPECT_EQ(g.getTarget(c), v);
  EXPECT_EQ(g.getSource(d), u);
  EXPECT_EQ(g.getTarget(d), w);
  EXPECT_EQ(g.getSource(e), v);
  EXPECT_EQ(g.getTarget(e), x);
  EXPECT_EQ(g.getSource(f), w);
  EXPECT_EQ(g.getTarget(f), x);
  const auto rc = g.getReverseEdge(c);
  const auto rd = g.getReverseEdge(d);
  const auto re = g.getReverseEdge(e);
  const auto rf = g.getReverseEdge(f);
  EXPECT_EQ(g.getReverseEdge(rc), c);
  EXPECT_EQ(g.getReverseEdge(rd), d);
  EXPECT_EQ(g.getReverseEdge(re), e);
  EXPECT_EQ(g.getReverseEdge(rf), f);
  EXPECT_TRUE(g.isBackwardEdge(rc));
  EXPECT_TRUE(g.isBackwardEdge(rd));
  EXPECT_TRUE(g.isBackwardEdge(re));
  EXPECT_TRUE(g.isBackwardEdge(rf));
  EXPECT_FALSE(g.isBackwardEdge(c));
  EXPECT_FALSE(g.isBackwardEdge(d));
  EXPECT_FALSE(g.isBackwardEdge(e));
  EXPECT_FALSE(g.isBackwardEdge(f));
  EXPECT_EQ(g.getSource(rc), v);
  EXPECT_EQ(g.getSource(rd), w);
  EXPECT_EQ(g.getSource(re), x);
  EXPECT_EQ(g.getSource(rf), x);
  EXPECT_EQ(g.getTarget(rc), u);
  EXPECT_EQ(g.getTarget(rd), u);
  EXPECT_EQ(g.getTarget(re), v);
  EXPECT_EQ(g.getTarget(rf), w);
}
TEST(FlowNetwork, DiamondUnordered) {
  //                 ┌───┐
  // ┌───┐     ┌────>│ v ├─────┐     ┌───┐
  // │ u ├─────┤     ├───┤     ├────>│ x │
  // └───┘     └────>│ w ├─────┘     └───┘
  //                 └───┘
  MinFlowScheduler::FlowNetwork g;
  const auto u = g.addVertex();
  const auto v = g.addVertex();
  const auto w = g.addVertex();
  const auto x = g.addVertex();
  auto f = g.addEdgeWithCapacityAndUnitCost(w, x, 1, 1);
  auto e = g.addEdgeWithCapacityAndUnitCost(v, x, 1, 1);
  auto d = g.addEdgeWithCapacityAndUnitCost(u, w, 1, 1);
  auto c = g.addEdgeWithCapacityAndUnitCost(u, v, 1, 1);
  Permutation permutation;
  EXPECT_NO_THROW(g.build(permutation));
  f = permutation[f];
  e = permutation[e];
  d = permutation[d];
  c = permutation[c];
  EXPECT_EQ(g.getNumVertices(), 4);
  EXPECT_EQ(g.getNumEdges(), 4);
  EXPECT_EQ(g.getForwardOutDegree(u), 2);
  EXPECT_EQ(g.getForwardOutDegree(v), 1);
  EXPECT_EQ(g.getForwardOutDegree(w), 1);
  EXPECT_EQ(g.getForwardOutDegree(x), 0);
  EXPECT_THAT(g.getForwardSuccessors(u), ::testing::UnorderedElementsAre(v, w));
  EXPECT_THAT(g.getForwardSuccessors(v), ::testing::UnorderedElementsAre(x));
  EXPECT_THAT(g.getForwardSuccessors(w), ::testing::UnorderedElementsAre(x));
  EXPECT_THAT(g.getForwardSuccessors(x), ::testing::IsEmpty());
  EXPECT_EQ(g.getSource(c), u);
  EXPECT_EQ(g.getTarget(c), v);
  EXPECT_EQ(g.getSource(d), u);
  EXPECT_EQ(g.getTarget(d), w);
  EXPECT_EQ(g.getSource(e), v);
  EXPECT_EQ(g.getTarget(e), x);
  EXPECT_EQ(g.getSource(f), w);
  EXPECT_EQ(g.getTarget(f), x);
  const auto rc = g.getReverseEdge(c);
  const auto rd = g.getReverseEdge(d);
  const auto re = g.getReverseEdge(e);
  const auto rf = g.getReverseEdge(f);
  EXPECT_EQ(g.getReverseEdge(rc), c);
  EXPECT_EQ(g.getReverseEdge(rd), d);
  EXPECT_EQ(g.getReverseEdge(re), e);
  EXPECT_EQ(g.getReverseEdge(rf), f);
  EXPECT_TRUE(g.isBackwardEdge(rc));
  EXPECT_TRUE(g.isBackwardEdge(rd));
  EXPECT_TRUE(g.isBackwardEdge(re));
  EXPECT_TRUE(g.isBackwardEdge(rf));
  EXPECT_FALSE(g.isBackwardEdge(c));
  EXPECT_FALSE(g.isBackwardEdge(d));
  EXPECT_FALSE(g.isBackwardEdge(e));
  EXPECT_FALSE(g.isBackwardEdge(f));
  EXPECT_EQ(g.getSource(rc), v);
  EXPECT_EQ(g.getSource(rd), w);
  EXPECT_EQ(g.getSource(re), x);
  EXPECT_EQ(g.getSource(rf), x);
  EXPECT_EQ(g.getTarget(rc), u);
  EXPECT_EQ(g.getTarget(rd), u);
  EXPECT_EQ(g.getTarget(re), v);
  EXPECT_EQ(g.getTarget(rf), w);
}
TEST(FlowNetwork, MaxFlowForward) {
  //                 ┌───┐
  // ┌───┐     ┌4/5─>│ v ├─4/4─┐     ┌───┐
  // │ u ├─────┤     ├───┤     ├────>│ x │
  // └───┘     └3/3─>│ w ├─3/4─┘     └───┘
  //                 └───┘
  MinFlowScheduler::FlowNetwork g;
  const auto u = g.addVertex();
  const auto v = g.addVertex();
  const auto w = g.addVertex();
  const auto x = g.addVertex();
  g.addEdgeWithCapacityAndUnitCost(u, v, 5, 1);
  g.addEdgeWithCapacityAndUnitCost(u, w, 3, 1);
  g.addEdgeWithCapacityAndUnitCost(v, x, 4, 1);
  g.addEdgeWithCapacityAndUnitCost(w, x, 4, 1);
  Permutation permutation;
  g.build(permutation);
  g.solveMaxFlow(u, x);
  EXPECT_EQ(g.getMaximumFlow(), 7);
}
TEST(FlowNetwork, MaxFlowCycle) {
  //                ┌────0/3─────┐
  //                v            │
  // ┌───┐        ┌───┐        ┌─┴─┐        ┌───┐
  // │ u ├──4/5──>│ v ├──4/6──>│ w ├──4/4──>│ x │
  // └───┘        └───┘        └───┘        └───┘
  MinFlowScheduler::FlowNetwork g;
  const auto u = g.addVertex();
  const auto v = g.addVertex();
  const auto w = g.addVertex();
  const auto x = g.addVertex();
  g.addEdgeWithCapacityAndUnitCost(u, v, 5, 1);
  g.addEdgeWithCapacityAndUnitCost(v, w, 6, 1);
  g.addEdgeWithCapacityAndUnitCost(w, x, 4, 1);
  g.addEdgeWithCapacityAndUnitCost(w, v, 3, 1);
  Permutation permutation;
  g.build(permutation);
  g.solveMaxFlow(u, x);
  EXPECT_EQ(g.getMaximumFlow(), 4);
}
TEST(FlowNetwork, MinCostMaxFlowForward) {
  //                                       ┌───┐
  // ┌───┐              ┌───┐     ┌2/6(2)─>│ w ├─2/5(0)─┐     ┌───┐
  // │ u ├───5/5(1)────>│ v ├─────┤        ├───┤        ├────>│ y │
  // └───┘              └───┘     └3/3(1)─>│ x ├─3/4(0)─┘     └───┘
  //                                       └───┘
  MinFlowScheduler::FlowNetwork g;
  const auto u = g.addVertex();
  const auto v = g.addVertex();
  const auto w = g.addVertex();
  const auto x = g.addVertex();
  const auto y = g.addVertex();
  auto a = g.addEdgeWithCapacityAndUnitCost(u, v, 5, 1);
  auto b = g.addEdgeWithCapacityAndUnitCost(v, w, 6, 2);
  auto c = g.addEdgeWithCapacityAndUnitCost(v, x, 3, 1);
  auto d = g.addEdgeWithCapacityAndUnitCost(w, y, 5, 0);
  auto e = g.addEdgeWithCapacityAndUnitCost(x, y, 4, 0);
  Permutation permutation;
  g.build(permutation);
  a = permutation[a];
  b = permutation[b];
  c = permutation[c];
  d = permutation[d];
  e = permutation[e];
  g.solveMinCostMaxFlow(u, y);
  EXPECT_EQ(g.getMaximumFlow(), 5);
  EXPECT_EQ(g.getFlow(a), 5);
  EXPECT_EQ(g.getFlow(b), 2);
  EXPECT_EQ(g.getFlow(c), 3);
  EXPECT_EQ(g.getFlow(d), 2);
  EXPECT_EQ(g.getFlow(e), 3);
}
} // namespace na::zoned
