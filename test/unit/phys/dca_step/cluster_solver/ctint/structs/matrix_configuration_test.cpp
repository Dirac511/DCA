// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests ct_int_matrix_configuration.hpp.

#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"

#include <array>
#include <cstdint>
#include <vector>

#include "gtest/gtest.h"

using dca::phys::solver::ctint::ConfigRef;
using dca::phys::solver::ctint::InteractionElement;
using dca::phys::solver::ctint::InteractionVertices;
using dca::phys::solver::ctint::MatrixConfiguration;
using dca::phys::solver::ctint::Vertex;

class MatrixConfigurationWrapper : public MatrixConfiguration {
public:
  using BaseClass = MatrixConfiguration;

  MatrixConfigurationWrapper(const InteractionVertices* H_int, int bands)
      : BaseClass(H_int, bands) {}
  using BaseClass::addVertex;
};

TEST(MatrixConfigurationTest, AddVertex) {
  InteractionVertices vertices;
  vertices.initialize(1., false);
  vertices.insertElement({{0, 0, 0, 0}, {0, 0, 2, 2}, 1.});  // up-down
  vertices.insertElement({{0, 0, 0, 0}, {0, 0, 1, 1}, 1.});  // up-up

  MatrixConfigurationWrapper config(&vertices, 2);

  std::array<std::vector<ConfigRef>, 2> config_refs;
  Vertex v1{false, 0, 1, 0.};
  config.addVertex(v1, 0, config_refs);

  EXPECT_EQ(1, config.size(0));
  EXPECT_EQ(1, config.size(1));
  EXPECT_EQ(0, v1.spins[0]);
  EXPECT_EQ(1, v1.spins[1]);
  EXPECT_EQ(1, config_refs[0].size());
  EXPECT_EQ(1, config_refs[1].size());

  Vertex v2{false, 1, 2, 0.1};
  config.addVertex(v2, 1, config_refs);

  EXPECT_EQ(3, config.size(0));
  EXPECT_EQ(1, config.size(1));
  EXPECT_EQ(0, v2.spins[0]);
  EXPECT_EQ(0, v2.spins[1]);
  EXPECT_EQ(3, config_refs[0].size());
  EXPECT_EQ(1, config_refs[1].size());
}
