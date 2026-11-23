// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Bilayer two-orbital Kanamori-Hubbard model for La3Ni2O7 (arXiv:2508.04554).
//
// Four bands: 0 = d_x2-y2 (layer 1), 1 = d_z2 (layer 1), 2 = d_x2-y2 (layer 2), 3 = d_z2 (layer 2).
//
//   H_0 = sum_k  eps_x(k) (n_0 + n_2) + eps_z(k) (n_1 + n_3)
//                + V(k) (c^+_0 c_1 + c^+_1 c_0 + c^+_2 c_3 + c^+_3 c_2)
//                + t_perp (c^+_1 c_3 + c^+_3 c_1)
//   eps_a(k) = e_a + 2 t_a (cos kx + cos ky) + 4 t_a' cos kx cos ky
//   V(k)     = V (cos kx - cos ky)
//
// The interaction is the intra-atomic Kanamori one of each layer (orbitals 0-1 and 2-3):
//   U      n_{a,up}  n_{a,down}
//   U'     n_{a,s}   n_{b,s'}          (s != s')
//   U'-J   n_{a,s}   n_{b,s}           (s == s')
//   J      spin flip  c^+_{a,up} c_{b,up} c^+_{b,down} c_{a,down}
//   J      pair hopping c^+_{a,up} c_{b,up} c^+_{a,down} c_{b,down}
// with U' = U - 2J.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_KANAMORI_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_KANAMORI_HPP

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
class BilayerKanamori {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
  typedef domains::no_symmetry<2> DCA_point_group;

  constexpr static int DIMENSION = 2;
  constexpr static int BANDS = 4;

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
double* BilayerKanamori<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];
  r_DCA[0] = 1.;
  r_DCA[1] = 0.;
  r_DCA[2] = 0.;
  r_DCA[3] = 1.;
  return r_DCA;
}

template <typename point_group_type>
double* BilayerKanamori<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];
  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;
  return r_LDA;
}

template <typename point_group_type>
std::vector<int> BilayerKanamori<point_group_type>::get_flavors() {
  std::vector<int> flavors(BANDS);
  for (int i = 0; i < BANDS; ++i)
    flavors[i] = i;
  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> BilayerKanamori<point_group_type>::get_a_vectors() {
  return std::vector<std::vector<double>>(BANDS, std::vector<double>(DIMENSION, 0.));
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>
BilayerKanamori<point_group_type>::get_orbital_permutations() {
  return {};
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void BilayerKanamori<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS or SpinDmn::dmn_size() != 2)
    throw std::logic_error("BilayerKanamori expects four bands and two spins.");

  const int origin = RDmn::parameter_type::origin_index();
  const double U = parameters.get_U();
  const double U_prime = parameters.get_U_prime();
  const double same_spin = U_prime - parameters.get_J();  // Kanamori relation

  H_interaction = 0.;

  // Intra-orbital, opposite spin.
  for (int b = 0; b < BANDS; ++b) {
    H_interaction(b, 0, b, 1, origin) = U;
    H_interaction(b, 1, b, 0, origin) = U;
  }

  // Intra-atomic inter-orbital terms: only between d_x2-y2 and d_z2 of the same layer.
  for (int layer = 0; layer < 2; ++layer) {
    const int a = 2 * layer;
    const int b = 2 * layer + 1;
    for (int s1 = 0; s1 < 2; ++s1)
      for (int s2 = 0; s2 < 2; ++s2) {
        const double value = (s1 == s2) ? same_spin : U_prime;
        H_interaction(a, s1, b, s2, origin) = value;
        H_interaction(b, s1, a, s2, origin) = value;
      }
  }
}

template <typename point_group_type>
template <typename Nu, typename RDmn, typename parameters_type>
void BilayerKanamori<point_group_type>::initializeNonDensityInteraction(
    func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, RDmn>>& non_density_interaction,
    const parameters_type& parameters) {
  const int origin = RDmn::parameter_type::origin_index();
  // Amplitude of the four-leg (exchange / pair-hopping) terms; it can be switched off
  // independently of the Kanamori relation used for the density-density part.
  const double J = parameters.get_J_exchange();
  auto nu = [](const int b, const int s) { return b + BANDS * s; };

  constexpr int up = 0, down = 1;

  non_density_interaction = 0.;
  for (int layer = 0; layer < 2; ++layer) {
    const int a0 = 2 * layer;
    const int a1 = 2 * layer + 1;
    for (int pair = 0; pair < 2; ++pair) {
      const int b1 = pair == 0 ? a0 : a1;
      const int b2 = pair == 0 ? a1 : a0;

      // Spin flip.
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b2, down), nu(b1, down), origin) = J;

      // Pair hopping.
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b1, down), nu(b2, down), origin) = J;
    }
  }
}

template <typename point_group_type>
template <class domain>
void BilayerKanamori<point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;
  for (int b = 0; b < BANDS; ++b)
    for (int s = 0; s < 2; ++s)
      H_symmetries(b + BANDS * s, b + BANDS * s) = b + BANDS * s;
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn,
          typename KDmn>
void BilayerKanamori<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS or SpinDmn::dmn_size() != 2)
    throw std::logic_error("BilayerKanamori expects four bands and two spins.");

  const auto& k_vecs = KDmn::get_elements();

  const double ex = parameters.get_ex();
  const double ez = parameters.get_ez();
  const double eb = parameters.get_eb();
  const double et = parameters.get_et();

  const double t11x = parameters.get_t11x();
  const double t11xy = parameters.get_t11xy();
  const double t11xx = parameters.get_t11xx();
  const double t22x = parameters.get_t22x();
  const double t22xy = parameters.get_t22xy();
  const double t22xx = parameters.get_t22xx();
  const double t12x = parameters.get_t12x();
  const double t12xx = parameters.get_t12xx();

  const double s110 = parameters.get_s110();
  const double s11x = parameters.get_s11x();
  const double s11xy = parameters.get_s11xy();
  const double s11xx = parameters.get_s11xx();
  const double s220 = parameters.get_s220();
  const double s22x = parameters.get_s22x();
  const double s22xy = parameters.get_s22xy();
  const double s22xx = parameters.get_s22xx();
  const double s12x = parameters.get_s12x();
  const double s12xx = parameters.get_s12xx();

  H_0 = ScalarType(0);
  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
    const double c1 = std::cos(k[0]);
    const double c2 = std::cos(k[1]);

    const ScalarType Hx0 = ex + eb + 2. * t11x * (c1 + c2) + 4. * t11xy * c1 * c2 +
                           2. * t11xx * (std::cos(2. * k[0]) + std::cos(2. * k[1]));
    const ScalarType Hz1 = ez + eb + 2. * t22x * (c1 + c2) + 4. * t22xy * c1 * c2 +
                           2. * t22xx * (std::cos(2. * k[0]) + std::cos(2. * k[1]));
    const ScalarType Hx2 = ex + et + 2. * t11x * (c1 + c2) + 4. * t11xy * c1 * c2 +
                           2. * t11xx * (std::cos(2. * k[0]) + std::cos(2. * k[1]));
    const ScalarType Hz3 = ez + et + 2. * t22x * (c1 + c2) + 4. * t22xy * c1 * c2 +
                           2. * t22xx * (std::cos(2. * k[0]) + std::cos(2. * k[1]));
    const ScalarType V = 2. * t12x * (c1 - c2) + 2. * t12xx * (std::cos(2. * k[0]) - std::cos(2. * k[1]));

    const ScalarType Hxp = s110 + 2. * s11x * (c1 + c2) + 4. * s11xy * c1 * c2 +
                           2. * s11xx * (std::cos(2. * k[0]) + std::cos(2. * k[1]));
    const ScalarType Hzp = s220 + 2. * s22x * (c1 + c2) + 4. * s22xy * c1 * c2 +
                           2. * s22xx * (std::cos(2. * k[0]) + std::cos(2. * k[1]));
    const ScalarType Vp = 2. * s12x * (c1 - c2) +
                          2. * s12xx * (std::cos(2. * k[0]) - std::cos(2. * k[1]));

    for (int s = 0; s < 2; ++s) {
      H_0(0, s, 0, s, k_ind) = Hx0;
      H_0(1, s, 1, s, k_ind) = Hz1;
      H_0(0, s, 1, s, k_ind) = V;
      H_0(1, s, 0, s, k_ind) = V;

      H_0(2, s, 2, s, k_ind) = Hx2;
      H_0(3, s, 3, s, k_ind) = Hz3;
      H_0(2, s, 3, s, k_ind) = V;
      H_0(3, s, 2, s, k_ind) = V;

      H_0(0, s, 2, s, k_ind) = Hxp;
      H_0(2, s, 0, s, k_ind) = Hxp;
      H_0(1, s, 3, s, k_ind) = Hzp;
      H_0(3, s, 1, s, k_ind) = Hzp;

      H_0(0, s, 3, s, k_ind) = Vp;
      H_0(3, s, 0, s, k_ind) = Vp;
      H_0(1, s, 2, s, k_ind) = Vp;
      H_0(2, s, 1, s, k_ind) = Vp;
    }
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_KANAMORI_HPP
