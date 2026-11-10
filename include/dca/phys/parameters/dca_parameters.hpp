// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the DCA parameters.

#ifndef DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/phys/parameters/self_energy_mixing_scheme.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class DcaParameters {
public:
  DcaParameters()
      : initial_self_energy_("zero"),
        dca_iterations_(1),
        dca_accuracy_(0.),
        self_energy_mixing_factor_(1.),
        self_energy_mixing_scheme_(SelfEnergyMixingScheme::LINEAR),
        self_energy_mixing_history_size_(8),
        self_energy_mixing_regularization_(1.e-8),
        self_energy_mixing_max_coefficient_(10.),
        self_energy_mixing_residual_weight_cutoff_(20),
        interacting_orbitals_{0},

        do_finite_size_qmc_(false),

        k_mesh_recursion_(0),
        coarsegraining_periods_(0),
        quadrature_rule_(1),
        coarsegraining_threads_(1),
        tail_frequencies_(0),

        do_dca_plus_(false),
        deconvolution_iterations_(16),
        deconvolution_tolerance_(1.e-3),
        hts_approximation_(false),
        hts_threads_(1) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_initial_self_energy() const {
    return initial_self_energy_;
  }
  int get_dca_iterations() const {
    return dca_iterations_;
  }
  void set_dca_iterations(const int iterations) {
    assert(iterations > 0);
    dca_iterations_ = iterations;
  }
  double get_dca_accuracy() const {
    return dca_accuracy_;
  }
  double get_self_energy_mixing_factor() const {
    return self_energy_mixing_factor_;
  }
  SelfEnergyMixingScheme get_self_energy_mixing_scheme() const {
    return self_energy_mixing_scheme_;
  }
  int get_self_energy_mixing_history_size() const {
    return self_energy_mixing_history_size_;
  }
  double get_self_energy_mixing_regularization() const {
    return self_energy_mixing_regularization_;
  }
  double get_self_energy_mixing_max_coefficient() const {
    return self_energy_mixing_max_coefficient_;
  }
  int get_self_energy_mixing_residual_weight_cutoff() const {
    return self_energy_mixing_residual_weight_cutoff_;
  }
  const std::vector<int>& get_interacting_orbitals() const {
    return interacting_orbitals_;
  }
  bool do_finite_size_qmc() const {
    return do_finite_size_qmc_;
  }
  int get_k_mesh_recursion() const {
    return k_mesh_recursion_;
  }
  int get_coarsegraining_periods() const {
    return coarsegraining_periods_;
  }
  int get_quadrature_rule() const {
    return quadrature_rule_;
  }
  int get_coarsegraining_threads() const {
    return coarsegraining_threads_;
  }
  int get_tail_frequencies() const {
    return tail_frequencies_;
  }
  bool do_dca_plus() const {
    return do_dca_plus_;
  }
  int get_deconvolution_iterations() const {
    return deconvolution_iterations_;
  }
  double get_deconvolution_tolerance() const {
    return deconvolution_tolerance_;
  }
  bool hts_approximation() const {
    return hts_approximation_;
  }
  int get_hts_threads() const {
    return hts_threads_;
  }

private:
  std::string initial_self_energy_;
  int dca_iterations_;
  double dca_accuracy_;
  double self_energy_mixing_factor_;
  SelfEnergyMixingScheme self_energy_mixing_scheme_;
  int self_energy_mixing_history_size_;
  // Anderson/Pulay mixer settings. In the earlier prototype these were hard-coded constants: the
  // regularization was 1.e-8, the coefficients had no bound, and the residual norm had uniform
  // weights. They are exposed as parameters so that the old behavior can be reproduced exactly
  // (regularization 1.e-8, max coefficient 1.e30, weight cutoff 0) or tuned.
  double self_energy_mixing_regularization_;
  double self_energy_mixing_max_coefficient_;
  int self_energy_mixing_residual_weight_cutoff_;
  std::vector<int> interacting_orbitals_;

  bool do_finite_size_qmc_;

  // coarse-graining
  int k_mesh_recursion_;
  int coarsegraining_periods_;
  int quadrature_rule_;
  int coarsegraining_threads_;
  int tail_frequencies_;

  // DCA+
  bool do_dca_plus_;
  int deconvolution_iterations_;
  double deconvolution_tolerance_;
  bool hts_approximation_;
  int hts_threads_;
};

template <typename Concurrency>
int DcaParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(initial_self_energy_);
  buffer_size += concurrency.get_buffer_size(dca_iterations_);
  buffer_size += concurrency.get_buffer_size(dca_accuracy_);
  buffer_size += concurrency.get_buffer_size(self_energy_mixing_factor_);
  buffer_size += concurrency.get_buffer_size(self_energy_mixing_scheme_);
  buffer_size += concurrency.get_buffer_size(self_energy_mixing_history_size_);
  buffer_size += concurrency.get_buffer_size(self_energy_mixing_regularization_);
  buffer_size += concurrency.get_buffer_size(self_energy_mixing_max_coefficient_);
  buffer_size += concurrency.get_buffer_size(self_energy_mixing_residual_weight_cutoff_);
  buffer_size += concurrency.get_buffer_size(interacting_orbitals_);
  buffer_size += concurrency.get_buffer_size(do_finite_size_qmc_);
  buffer_size += concurrency.get_buffer_size(k_mesh_recursion_);
  buffer_size += concurrency.get_buffer_size(coarsegraining_periods_);
  buffer_size += concurrency.get_buffer_size(quadrature_rule_);
  buffer_size += concurrency.get_buffer_size(coarsegraining_threads_);
  buffer_size += concurrency.get_buffer_size(tail_frequencies_);
  buffer_size += concurrency.get_buffer_size(do_dca_plus_);
  buffer_size += concurrency.get_buffer_size(deconvolution_iterations_);
  buffer_size += concurrency.get_buffer_size(deconvolution_tolerance_);
  buffer_size += concurrency.get_buffer_size(hts_approximation_);
  buffer_size += concurrency.get_buffer_size(hts_threads_);

  return buffer_size;
}

template <typename Concurrency>
void DcaParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                         int& position) const {
  concurrency.pack(buffer, buffer_size, position, initial_self_energy_);
  concurrency.pack(buffer, buffer_size, position, dca_iterations_);
  concurrency.pack(buffer, buffer_size, position, dca_accuracy_);
  concurrency.pack(buffer, buffer_size, position, self_energy_mixing_factor_);
  concurrency.pack(buffer, buffer_size, position, self_energy_mixing_scheme_);
  concurrency.pack(buffer, buffer_size, position, self_energy_mixing_history_size_);
  concurrency.pack(buffer, buffer_size, position, self_energy_mixing_regularization_);
  concurrency.pack(buffer, buffer_size, position, self_energy_mixing_max_coefficient_);
  concurrency.pack(buffer, buffer_size, position, self_energy_mixing_residual_weight_cutoff_);
  concurrency.pack(buffer, buffer_size, position, interacting_orbitals_);
  concurrency.pack(buffer, buffer_size, position, do_finite_size_qmc_);
  concurrency.pack(buffer, buffer_size, position, k_mesh_recursion_);
  concurrency.pack(buffer, buffer_size, position, coarsegraining_periods_);
  concurrency.pack(buffer, buffer_size, position, quadrature_rule_);
  concurrency.pack(buffer, buffer_size, position, coarsegraining_threads_);
  concurrency.pack(buffer, buffer_size, position, tail_frequencies_);
  concurrency.pack(buffer, buffer_size, position, do_dca_plus_);
  concurrency.pack(buffer, buffer_size, position, deconvolution_iterations_);
  concurrency.pack(buffer, buffer_size, position, deconvolution_tolerance_);
  concurrency.pack(buffer, buffer_size, position, hts_approximation_);
  concurrency.pack(buffer, buffer_size, position, hts_threads_);
}

template <typename Concurrency>
void DcaParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                           int& position) {
  concurrency.unpack(buffer, buffer_size, position, initial_self_energy_);
  concurrency.unpack(buffer, buffer_size, position, dca_iterations_);
  concurrency.unpack(buffer, buffer_size, position, dca_accuracy_);
  concurrency.unpack(buffer, buffer_size, position, self_energy_mixing_factor_);
  concurrency.unpack(buffer, buffer_size, position, self_energy_mixing_scheme_);
  concurrency.unpack(buffer, buffer_size, position, self_energy_mixing_history_size_);
  concurrency.unpack(buffer, buffer_size, position, self_energy_mixing_regularization_);
  concurrency.unpack(buffer, buffer_size, position, self_energy_mixing_max_coefficient_);
  concurrency.unpack(buffer, buffer_size, position, self_energy_mixing_residual_weight_cutoff_);
  concurrency.unpack(buffer, buffer_size, position, interacting_orbitals_);
  concurrency.unpack(buffer, buffer_size, position, do_finite_size_qmc_);
  concurrency.unpack(buffer, buffer_size, position, k_mesh_recursion_);
  concurrency.unpack(buffer, buffer_size, position, coarsegraining_periods_);
  concurrency.unpack(buffer, buffer_size, position, quadrature_rule_);
  concurrency.unpack(buffer, buffer_size, position, coarsegraining_threads_);
  concurrency.unpack(buffer, buffer_size, position, tail_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, do_dca_plus_);
  concurrency.unpack(buffer, buffer_size, position, deconvolution_iterations_);
  concurrency.unpack(buffer, buffer_size, position, deconvolution_tolerance_);
  concurrency.unpack(buffer, buffer_size, position, hts_approximation_);
  concurrency.unpack(buffer, buffer_size, position, hts_threads_);
}

template <typename ReaderOrWriter>
void DcaParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  auto try_to_read = [&](const std::string& name, auto& var) {
    try {
      reader_or_writer.execute(name, var);
    }
    catch (const std::exception& r_e) {
    }
  };

  // The self-energy mixing scheme is converted only after the try block below, so that an unknown
  // name is reported to the user instead of being swallowed by the exception handling.
  bool has_mixing_scheme = false;
  std::string mixing_scheme_str;

  try {
    reader_or_writer.open_group("DCA");

    try_to_read("initial-self-energy", initial_self_energy_);
    try_to_read("iterations", dca_iterations_);
    try_to_read("accuracy", dca_accuracy_);
    try_to_read("self-energy-mixing-factor", self_energy_mixing_factor_);

    // Read the self-energy mixing scheme. A missing entry keeps the default, but an unknown name
    // must not be silently replaced by the default (see the conversion after the try block).
    mixing_scheme_str = toString(self_energy_mixing_scheme_);
    try {
      reader_or_writer.execute("mixing-scheme", mixing_scheme_str);
      has_mixing_scheme = true;
    }
    catch (const std::exception& r_e) {
    }

    // Numerical settings of the Anderson/Pulay mixer. A file written for the earlier prototype
    // does not contain them and keeps the defaults below. Those defaults match the old hard-coded
    // values for the regularization (1.e-8) and for the weights (20 positive Matsubara
    // frequencies), but the coefficient bound (10) is new: the old code had no bound at all.
    // Setting "mixing-max-coefficient" to a very large value and
    // "mixing-residual-weight-cutoff" to 0 reproduces the old behavior exactly.
    try_to_read("mixing-history-size", self_energy_mixing_history_size_);
    try_to_read("mixing-regularization", self_energy_mixing_regularization_);
    try_to_read("mixing-max-coefficient", self_energy_mixing_max_coefficient_);
    try_to_read("mixing-residual-weight-cutoff", self_energy_mixing_residual_weight_cutoff_);

    try_to_read("interacting-orbitals", interacting_orbitals_);

    try_to_read("do-finite-size-QMC", do_finite_size_qmc_);

    try {
      reader_or_writer.open_group("coarse-graining");

      try_to_read("k-mesh-recursion", k_mesh_recursion_);
      try_to_read("periods", coarsegraining_periods_);
      try_to_read("quadrature-rule", quadrature_rule_);
      try_to_read("threads", coarsegraining_threads_);
      try_to_read("tail-frequencies", tail_frequencies_);

      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.open_group("DCA+");

      try_to_read("do-DCA+", do_dca_plus_);
      try_to_read("deconvolution-iterations", deconvolution_iterations_);
      try_to_read("deconvolution-tolerance", deconvolution_tolerance_);
      try_to_read("HTS-approximation", hts_approximation_);
      try_to_read("HTS-threads", hts_threads_);

      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }

  if (has_mixing_scheme)
    self_energy_mixing_scheme_ = stringToSelfEnergyMixingScheme(mixing_scheme_str);

  // Check read parameters for consistency.
  if (reader_or_writer.is_reader()) {
    if (do_finite_size_qmc_ && do_dca_plus_)
      throw std::logic_error("Finite-size QMC and DCA+ are mutually exclusive options.");

    if (self_energy_mixing_history_size_ < 1)
      throw std::logic_error("The self-energy mixing history size must be at least 1.");
    if (self_energy_mixing_regularization_ < 0.)
      throw std::logic_error("The self-energy mixing regularization must not be negative.");
    if (self_energy_mixing_max_coefficient_ < 1.)
      throw std::logic_error("The maximal self-energy mixing coefficient must be at least 1.");
    if (self_energy_mixing_residual_weight_cutoff_ < 0)
      throw std::logic_error(
          "The cutoff of the self-energy mixing residual weights must not be negative.");
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP
