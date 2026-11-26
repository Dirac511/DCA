// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Minimal two-orbital Kanamori lattice used by the CT-INT tests.
//
// The interaction is local and contains the density-density terms
//   U  n_{a,up} n_{a,down}
//   U' n_{a,s}  n_{b,s'}
// and the non-density pair-hopping/spin-flip terms proportional to J.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_KANAMORI_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_KANAMORI_LATTICE_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class KanamoriLattice {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
  typedef domains::no_symmetry<2> DCA_point_group;

  constexpr static int DIMENSION = 2;
  constexpr static int BANDS = 2;

  static double* initialize_r_DCA_basis();
  static double* initialize_r_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initialize_H_interaction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <typename Nu, typename RDmn, typename parameters_type>
  static void initializeNonDensityInteraction(
      func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, RDmn>>& non_density_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn,
            typename KDmn>
  static void initialize_H_0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);
};

template <typename point_group_type>
double* KanamoriLattice<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];
  r_DCA[0] = 1.;
  r_DCA[1] = 0.;
  r_DCA[2] = 0.;
  r_DCA[3] = 1.;
  return r_DCA;
}

template <typename point_group_type>
double* KanamoriLattice<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];
  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;
  return r_LDA;
}

template <typename point_group_type>
std::vector<int> KanamoriLattice<point_group_type>::get_flavors() {
  std::vector<int> flavors(BANDS);
  for (int i = 0; i < BANDS; ++i)
    flavors[i] = i;
  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> KanamoriLattice<point_group_type>::get_a_vectors() {
  return std::vector<std::vector<double>>(BANDS, std::vector<double>(DIMENSION, 0.));
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>
KanamoriLattice<point_group_type>::get_orbital_permutations() {
  return {};
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void KanamoriLattice<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS or SpinDmn::dmn_size() != 2)
    throw std::logic_error("KanamoriLattice expects two bands and two spins.");

  const int origin = RDmn::parameter_type::origin_index();
  const double U = parameters.get_U();
  const double U_prime = parameters.get_U_prime();
  const double J = parameters.get_J();

  H_interaction = 0.;

  // The interaction function is indexed as (band_1, spin_1, band_2, spin_2, delta_r).
  // Intra-orbital opposite spin.
  for (int b = 0; b < BANDS; ++b) {
    H_interaction(b, 0, b, 1, origin) = U;
    H_interaction(b, 1, b, 0, origin) = U;
  }

  // Inter-orbital opposite spin.
  for (int b1 = 0; b1 < BANDS; ++b1)
    for (int b2 = 0; b2 < BANDS; ++b2)
      if (b1 != b2) {
        H_interaction(b1, 0, b2, 1, origin) = U_prime;
        H_interaction(b2, 1, b1, 0, origin) = U_prime;
        H_interaction(b1, 1, b2, 0, origin) = U_prime;
        H_interaction(b2, 0, b1, 1, origin) = U_prime;
      }

  // Inter-orbital same spin: U' - J.
  for (int b1 = 0; b1 < BANDS; ++b1)
    for (int b2 = 0; b2 < BANDS; ++b2)
      if (b1 != b2) {
        H_interaction(b1, 0, b2, 0, origin) = U_prime - J;
        H_interaction(b2, 0, b1, 0, origin) = U_prime - J;
        H_interaction(b1, 1, b2, 1, origin) = U_prime - J;
        H_interaction(b2, 1, b1, 1, origin) = U_prime - J;
      }
}

template <typename point_group_type>
template <typename Nu, typename RDmn, typename parameters_type>
void KanamoriLattice<point_group_type>::initializeNonDensityInteraction(
    func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, RDmn>>& non_density_interaction,
    const parameters_type& parameters) {
  const int origin = RDmn::parameter_type::origin_index();
  const double J = parameters.get_J_non_density();
  auto nu = [](int b, int s) { return b + BANDS * s; };

  // The consumer expects each leg to be a spin conserving hopping, i.e. the two indices of a leg
  // belong to the same spin.  With H = sum H(nu1,nu2,nu3,nu4) c^+_{nu1} c_{nu2} c^+_{nu3} c_{nu4}
  // this makes the pair hopping c^+_{b1 up} c_{b2 up} c^+_{b1 down} c_{b2 down} and the spin flip
  // c^+_{b1 up} c_{b2 up} c^+_{b2 down} c_{b1 down}.  The loop over both orders of (b1, b2)
  // generates each term together with its hermitian conjugate.
  constexpr int up = 0, down = 1;

  non_density_interaction = 0.;
  for (int b1 = 0; b1 < BANDS; ++b1)
    for (int b2 = 0; b2 < BANDS; ++b2) {
      if (b1 == b2)
        continue;

      // Spin flip.
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b2, down), nu(b1, down), origin) = J;

      // Pair hopping.
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b1, down), nu(b2, down), origin) = J;
    }
}

template <typename point_group_type>
template <class domain>
void KanamoriLattice<point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;
  for (int b = 0; b < BANDS; ++b)
    for (int s = 0; s < 2; ++s)
      H_symmetries(b + BANDS * s, b + BANDS * s) = b + BANDS * s;
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn,
          typename KDmn>
void KanamoriLattice<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS or SpinDmn::dmn_size() != 2)
    throw std::logic_error("KanamoriLattice expects two bands and two spins.");

  const auto& k_vecs = KDmn::get_elements();
  const auto t = parameters.get_t();

  H_0 = ScalarType(0);
  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
    const ScalarType dispersion =
        -2. * t * (std::cos(k[0]) + std::cos(k[1]));
    for (int b = 0; b < BANDS; ++b)
      for (int s = 0; s < 2; ++s)
        H_0(b, s, b, s, k_ind) = dispersion;
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_KANAMORI_LATTICE_HPP
