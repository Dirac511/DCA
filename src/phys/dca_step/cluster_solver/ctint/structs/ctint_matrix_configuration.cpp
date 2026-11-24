// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This file implements ct_int_matrix_configuration.hpp.

#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"

#include <utility>

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

MatrixConfiguration::MatrixConfiguration(const InteractionVertices* H_int, const int bands)
    : H_int_(H_int), n_bands_(bands), sectors_{Sector(), Sector()} {}

MatrixConfiguration& MatrixConfiguration::operator=(const MatrixConfiguration& rhs) {
  n_bands_ = rhs.n_bands_;
  H_int_ = rhs.H_int_;
  sectors_ = rhs.sectors_;
  n_vertices_ = rhs.n_vertices_;
  return *this;
}

MatrixConfiguration& MatrixConfiguration::operator=(MatrixConfiguration&& rhs) {
  n_bands_ = rhs.n_bands_;
  H_int_ = rhs.H_int_;
  sectors_ = std::move(rhs.sectors_);
  n_vertices_ = rhs.n_vertices_;
  return *this;
}

void MatrixConfiguration::addVertex(Vertex& v, unsigned config_id,
                                    std::array<std::vector<ConfigRef>, 2>& config_refs) {
  auto spin = [=](const int nu) { return nu >= n_bands_; };
  auto band = [=](const int nu) -> unsigned short { return nu - n_bands_ * spin(nu); };

  const auto& nu = (*H_int_)[v.interaction_id].nu;
  const auto& r = (*H_int_)[v.interaction_id].r;

  const bool is_ndd = band(nu[0]) != band(nu[1]) || band(nu[2]) != band(nu[3]);

  auto field_type = [&](const Vertex& v, const int leg) -> short {
    const short sign = v.aux_spin ? 1 : -1;
    const InteractionElement& elem = (*H_int_)[v.interaction_id];
    if (is_ndd)
      return leg == 1 ? -3 * sign : 3 * sign;  // non density-density.
    else if (elem.w > 0)
      return leg == 1 ? -1 * sign : 1 * sign;  // positive dd interaction.
    else
      return 2 * sign;  // negative dd interaction.
  };

  v.n_legs = 0;
  for (unsigned short leg = 0; leg < 2; ++leg) {
    assert(spin(nu[0 + 2 * leg]) == spin(nu[1 + 2 * leg]));
    const short s = spin(nu[0 + 2 * leg]);
    const unsigned short b_left = band(nu[0 + 2 * leg]);
    const unsigned short b_right = band(nu[1 + 2 * leg]);
    const short type = field_type(v, leg);

    // Every vertex occupies a single entry per spin sector, as in ALPS/CT-INT: for a four-leg vertex
    // the entry is the off-diagonal operator pair (b_left, b_right) and its D matrix element is the
    // inter-band Green function G(b_left, b_right) minus the per-leg alpha.  Density-density legs
    // are band diagonal (b_left == b_right).
    {
      // ALPS alpha[af_state][rank]: the auxiliary field state selects the row, not the spin sector.
      const int af_state = v.aux_spin ? 1 : 0;
      Sector& sector = sectors_[s];
      v.spins[v.n_legs] = s;
      v.matrix_config_indices[v.n_legs] = sector.size();
      sector.entries_.emplace_back(details::SectorEntry{b_left, r[0 + 2 * leg], b_right,
                                                        r[1 + 2 * leg], v.tau, type, config_id,
                                                        is_ndd ? v.alpha[af_state][leg] : 0.});
      config_refs[s].emplace_back(config_id, v.n_legs);
      ++v.n_legs;
    }
  }
  ++n_vertices_;
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
