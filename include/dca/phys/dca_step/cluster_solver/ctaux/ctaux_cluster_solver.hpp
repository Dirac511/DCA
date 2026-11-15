// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Cluster Monte Carlo integrator based on a continuous-time auxilary field (CT-AUX) expansion.
//
// TODO: Cleanup the computation of Sigma, error bars, etc. and have the same work flow independent
//       of whether the thread jacket (stdthread qmci) is used.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_CLUSTER_SOLVER_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/statistics/util.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_walker.hpp"
#include "dca/phys/dca_step/cluster_solver/self_energy_mixer.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/phys/parameters/self_energy_mixing_scheme.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
class CtauxClusterSolver {
protected:
  using DataType = Data;
  using ParametersType = Parameters;

  using Rng = typename Parameters::random_number_generator;

  using Profiler = typename Parameters::profiler_type;
  using Concurrency = typename Parameters::concurrency_type;

  using Walker = ctaux::CtauxWalker<device_t, Parameters, Data>;
  using Accumulator = ctaux::CtauxAccumulator<device_t, Parameters, Data>;

  static constexpr linalg::DeviceType device = device_t;

private:
  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  using NuNuKClusterWDmn = func::dmn_variadic<nu, nu, KClusterDmn, w>;
  using NuNuRClusterWDmn = func::dmn_variadic<nu, nu, RClusterDmn, w>;

public:
  CtauxClusterSolver(Parameters& parameters_ref, Data& MOMS_ref);

  template <typename Writer>
  void write(Writer& writer);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

  // Computes and returns the local value of the Green's function G(k, \omega), i.e. without
  // averaging it across processes.
  // For testing purposes.
  // Precondition: The accumulator_ data has not been averaged, i.e. finalize has not been called.
  auto local_G_k_w() const;

protected:
  void warmUp(Walker& walker);

  void measure(Walker& walker);

  void computeErrorBars();

private:
  void symmetrize_measurements();

  // Sums/averages the quantities measured by the individual MPI ranks.
  void collect_measurements();

  void compute_G_k_w_from_M_r_w();

  double compute_S_k_w_from_G_k_w();

  void compute_G_k_w_new(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& M_k_w_new,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new) const;

  void compute_S_k_w_new(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& S_k_w_new);

  void set_non_interacting_bands_to_zero();

  void adjust_self_energy_for_double_counting();

  double mix_self_energy(double alpha);

  // Norms of the fixed-point residual Sigma_QMC - Sigma_cg. They are evaluated before mixing so
  // that the convergence criterion does not depend on the mixing scheme.
  double self_energy_residual_l2_norm() const;
  double self_energy_residual_infty_norm() const;

protected:
  Parameters& parameters_;
  Data& data_;
  Concurrency& concurrency_;

  double total_time_;

  Accumulator accumulator_;

  int dca_iteration_;

private:
  Rng rng_;

  double thermalization_time_;
  double mc_integration_time_;

  func::function<std::complex<double>, NuNuKClusterWDmn> Sigma_old_;
  func::function<std::complex<double>, NuNuKClusterWDmn> Sigma_new_;
  // Mixed self-energy with the smallest residual seen so far. When an iteration makes L2 grow
  // again, this is the state that is written to the output instead of the rejected one.
  func::function<std::complex<double>, NuNuKClusterWDmn> best_Sigma_;
  // Un-mixed QMC quantities of the same iteration as best_Sigma_. finalize() derives
  // Sigma_zero_moment and the error bars from them, so they have to be restored together with
  // best_Sigma_, otherwise those diagnostics describe the rejected iteration.
  func::function<std::complex<double>, NuNuKClusterWDmn> best_Sigma_new_;
  func::function<std::complex<double>, NuNuKClusterWDmn> best_Sigma_stdv_;
  func::function<std::complex<double>, NuNuKClusterWDmn> best_G_k_w_stdv_;
  // Cluster Green function measured in the same iteration as best_Sigma_. It is overwritten by
  // compute_G_k_w_from_M_r_w() in every iteration, so without this snapshot the output would pair
  // the best self-energy with the QMC data of the iteration that was just rejected.
  func::function<std::complex<double>, NuNuKClusterWDmn> best_G_k_w_;
  func::function<std::complex<double>, NuNuRClusterWDmn> best_G_r_w_error_;
  // Raw measurements and the previous-iteration self-energy of the saved iteration. They are
  // overwritten in every iteration too, so restoring Sigma alone would not make everything written
  // to the output belong to the saved iteration.
  func::function<std::complex<double>, NuNuKClusterWDmn> best_Sigma_old_;
  func::function<std::complex<double>, NuNuRClusterWDmn> best_M_r_w_;
  func::function<std::complex<double>, NuNuRClusterWDmn> best_M_r_w_squared_;

  SelfEnergyMixer self_energy_mixer_;

  // L2 residual of the previous iteration, used to detect that the Anderson extrapolation started
  // to increase it again, and the resulting stop request.
  double previous_L2_error_ = 0.;
  bool stop_requested_ = false;
  // DCA iteration in which best_L2_error_ was measured, i.e. the iteration whose result is written
  // to the output when the loop is stopped. -1 means no entry yet.
  int best_dca_iteration_ = -1;
  // Smallest residual seen so far. It is restored, together with the matching self-energy and
  // diagnostics, when the loop is stopped.
  double best_L2_error_ = std::numeric_limits<double>::infinity();
  // Chemical potential that belongs to best_Sigma_, and the chemical potential of the iteration
  // that is currently running. The DCA loop updates mu before every iteration, so restoring Sigma
  // without restoring mu would pair the saved self-energy with the chemical potential of a later
  // iteration and make the density and n(k) recomputed afterwards inconsistent with it.
  double best_chemical_potential_ = 0.;
  double chemical_potential_at_start_ = 0.;
  // True if best_G_r_w_error_ holds a jackknife estimate. The estimate is only available for the
  // last DCA iteration, so most saved iterations have none; restoring an empty snapshot would
  // overwrite a valid error bar with zeros.
  bool best_has_error_bars_ = false;
  // Average sign and average expansion order of the QMC run of this iteration and of the saved
  // one. They are diagnostics of the QMC run, not derived from Sigma, so they need their own
  // snapshots: the accumulator is reset in every iteration and cannot be restored. The saved
  // values are reported in the stop message; the row of a rejected iteration is not written.
  // Average sign of the QMC run of this iteration, reduced over all ranks in
  // collect_measurements(). It is the one value used by the log message and the output file, so
  // both report the same number regardless of which rank writes the file.
  double average_sign_ = 1.;
  double current_sign_ = 1.;
  double current_expansion_order_ = 0.;
  double best_sign_ = 1.;
  double best_expansion_order_ = 0.;

  double accumulated_sign_;
  func::function<std::complex<double>, NuNuRClusterWDmn> M_r_w_;
  func::function<std::complex<double>, NuNuRClusterWDmn> M_r_w_squared_;

  bool averaged_;
  bool compute_jack_knife_;
};

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
CtauxClusterSolver<device_t, Parameters, Data>::CtauxClusterSolver(Parameters& parameters_ref,
                                                                   Data& data_ref)
    : parameters_(parameters_ref),
      data_(data_ref),
      concurrency_(parameters_.get_concurrency()),

      total_time_(0),

      accumulator_(parameters_, data_, 0),

      dca_iteration_(-1),

      rng_(concurrency_.id(), concurrency_.number_of_processors(), parameters_.get_seed()),

      thermalization_time_(0),
      mc_integration_time_(0),

      Sigma_old_("Self-Energy-n-1-iteration"),
      Sigma_new_("Self-Energy-n-0-iteration"),
      best_Sigma_("Self-Energy-best-iteration"),
      best_Sigma_new_("Self-Energy-best-iteration-QMC"),
      best_Sigma_stdv_("Self-Energy-best-iteration-stddev"),
      best_G_k_w_stdv_("G_k_w-best-iteration-stddev"),
      best_G_k_w_("G_k_w-best-iteration"),
      best_G_r_w_error_("G_r_w-best-iteration-error"),
      best_Sigma_old_("Self-Energy-n-1-best-iteration"),
      best_M_r_w_("M_r_w-best-iteration"),
      best_M_r_w_squared_("M_r_w_squared-best-iteration"),

      M_r_w_("M_r_w"),
      M_r_w_squared_("M_r_w_squared"),

      averaged_(false) {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t CT-AUX Integrator is born \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <typename Writer>
void CtauxClusterSolver<device_t, Parameters, Data>::write(Writer& writer) {
  writer.open_group("CT-AUX-SOLVER-functions");

  writer.execute(Sigma_old_);
  writer.execute(Sigma_new_);

  accumulator_.write(writer);

  writer.close_group();
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::initialize(int dca_iteration) {
  dca_iteration_ = dca_iteration;

  Sigma_old_ = data_.Sigma;

  if (dca_iteration == 0) {
    self_energy_mixer_.reset();
    previous_L2_error_ = 0.;
    stop_requested_ = false;
    best_L2_error_ = std::numeric_limits<double>::infinity();
    best_dca_iteration_ = -1;
  }

  // The DCA loop has already updated mu for this iteration, so this is the value the QMC run and
  // the mixing of this iteration use.
  chemical_potential_at_start_ = parameters_.get_chemical_potential();

  SelfEnergyMixer::Config mixing_config;
  mixing_config.alpha = parameters_.get_self_energy_mixing_factor();
  mixing_config.max_history = parameters_.get_self_energy_mixing_history_size();
  mixing_config.regularization = parameters_.get_self_energy_mixing_regularization();
  mixing_config.max_coefficient = parameters_.get_self_energy_mixing_max_coefficient();
  self_energy_mixer_.configure(mixing_config);

  accumulator_.initialize(dca_iteration_);

  averaged_ = false;
  compute_jack_knife_ =
      (dca_iteration == parameters_.get_dca_iterations() - 1) &&
      (parameters_.get_error_computation_type() == ErrorComputationType::JACK_KNIFE);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t CT-AUX Integrator has initialized (DCA-iteration : " << dca_iteration
              << ")\n\n";
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::integrate() {
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "QMC integration has started: " << dca::util::print_time() << std::endl;
  }

  Walker walker(parameters_, data_, rng_, 0);

  walker.initialize();

  {
    dca::profiling::WallTime start_time;

    warmUp(walker);

    dca::profiling::WallTime mid_time;

    measure(walker);

    dca::profiling::WallTime end_time;

    dca::profiling::Duration ther_time(mid_time, start_time);
    dca::profiling::Duration meas_time(end_time, mid_time);

    dca::profiling::Duration tot_time(end_time, start_time);

    thermalization_time_ = ther_time.sec + 1.e-6 * ther_time.usec;
    mc_integration_time_ = meas_time.sec + 1.e-6 * meas_time.usec;
    total_time_ = tot_time.sec + 1.e-6 * tot_time.usec;
  }

  accumulator_.get_error_distribution() += walker.get_error_distribution();

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "On-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << parameters_.get_measurements() << std::endl;

    walker.printSummary();
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
template <typename dca_info_struct_t>
double CtauxClusterSolver<device_t, Parameters, Data>::finalize(dca_info_struct_t& dca_info_struct) {
  collect_measurements();
  symmetrize_measurements();

  // Compute new Sigma.
  compute_G_k_w_from_M_r_w();

  // FT<k_DCA,r_DCA>::execute(data_.G_k_w, data_.G_r_w);
  math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(data_.G_k_w, data_.G_r_w);

  // The average sign and the average expansion order describe the QMC run of this iteration. They
  // are evaluated before compute_S_k_w_from_G_k_w() because that call can stop the loop.
  {
    double total = 1.e-6, integral = 0;

    for (int l = 0; l < accumulator_.get_visited_expansion_order_k().size(); l++) {
      total += accumulator_.get_visited_expansion_order_k()(l);
      integral += accumulator_.get_visited_expansion_order_k()(l) * l;
    }

    current_expansion_order_ = integral / total;
    // Computed once over all ranks in collect_measurements(), which runs before this point.
    current_sign_ = average_sign_;
  }

  const double L2_Sigma_difference = compute_S_k_w_from_G_k_w();

  if (compute_jack_knife_ && parameters_.isAccumulatingG4()) {
    for (std::size_t channel = 0; channel < data_.get_G4_error().size(); ++channel)
      data_.get_G4_error()[channel] = concurrency_.jackknifeError(data_.get_G4()[channel], true);
  }

  // When the loop is stopped, only the saved (best) iteration is kept: the row of the rejected
  // iteration and the rows of the iterations after the saved one are dropped, so that nothing is
  // repeated and the last row of the output is the best result. The DCA loop skips its own
  // per-iteration output as well when the solver returns the stop signal below.
  if (!stop_requested_) {
    dca_info_struct.L2_Sigma_difference(dca_iteration_) = L2_Sigma_difference;

    for (int i = 0; i < b::dmn_size() * s::dmn_size(); i++) {
      for (int j = 0; j < KClusterDmn::dmn_size(); j++) {
        std::vector<double> x;
        for (int l = 0; l < w::dmn_size() / 4; l++)
          x.push_back(real(Sigma_new_(i, i, j, l)));

        dca_info_struct.Sigma_zero_moment(i, j, dca_iteration_) =
            math::statistics::util::mean(x);  // real(data_.Sigma(i,i,j,0));
        dca_info_struct.standard_deviation(i, j, dca_iteration_) =
            math::statistics::util::standard_deviation(x);
      }
    }

    dca_info_struct.average_expansion_order(dca_iteration_) = current_expansion_order_;
    dca_info_struct.sign(dca_iteration_) = current_sign_;

    dca_info_struct.thermalization_per_mpi_task(dca_iteration_) =
        thermalization_time_ / double(concurrency_.number_of_processors());
    dca_info_struct.MC_integration_per_mpi_task(dca_iteration_) =
        mc_integration_time_ / double(concurrency_.number_of_processors());

    dca_info_struct.times_per_mpi_task(dca_iteration_) =
        total_time_ / double(concurrency_.number_of_processors());
    dca_info_struct.Gflop_per_mpi_task(dca_iteration_) =
        accumulator_.get_Gflop() / double(concurrency_.number_of_processors());

    dca_info_struct.Gflops_per_mpi_task(dca_iteration_) =
        dca_info_struct.Gflop_per_mpi_task(dca_iteration_) /
        dca_info_struct.times_per_mpi_task(dca_iteration_);
  }
  else
    dca_info_struct.clear_after(best_dca_iteration_);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t CT-AUX Integrator has finalized \n" << std::endl;

  // A negative return value is the stop signal of the solver: an extrapolated step made the L2
  // residual grow again, so further DCA iterations would only accumulate harmful extrapolations.
  // The output was truncated to the saved iteration, so it ends with the best result.
  if (stop_requested_)
    return -1.;

  return L2_Sigma_difference;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::warmUp(Walker& walker) {
  Profiler profiler("thermalization", "QMCI", __LINE__);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t warm-up has started\n" << std::endl;

  for (int i = 0; i < parameters_.get_warm_up_sweeps(); i++) {
    walker.doSweep();
    walker.updateShell(i, parameters_.get_warm_up_sweeps());
  }

  walker.is_thermalized() = true;

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t warm-up has ended\n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::measure(Walker& walker) {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t measuring has started \n" << std::endl;

  const int n_meas = parallel::util::getWorkload(parameters_.get_measurements(), concurrency_);

  for (int i = 0; i < n_meas; i++) {
    {
      Profiler profiler("updating", "QMCI", __LINE__);
      walker.doSweep();
    }

    {
      Profiler profiler("measurements", "QMCI", __LINE__);
      accumulator_.updateFrom(walker);
      accumulator_.measure();
    }

    walker.updateShell(i, n_meas);
  }

  accumulator_.finalize();

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t measuring has ended \n" << std::endl;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::computeErrorBars() {
  if (!accumulator_.compute_std_deviation())
    return;
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t compute-error-bars on Self-energy\t" << dca::util::print_time() << "\n\n";

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w_new(
      "G_k_w_new");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> M_r_w_new(
      "M_r_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> M_k_w_new(
      "M_k_w_new");

  accumulator_.finalize();

  M_r_w_new = accumulator_.get_sign_times_M_r_w();
  M_r_w_new /= accumulator_.get_accumulated_sign();

  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(M_r_w_new, M_k_w_new);

  compute_G_k_w_new(M_k_w_new, G_k_w_new);
  compute_S_k_w_new(G_k_w_new, Sigma_new_);

  concurrency_.average_and_compute_stddev(Sigma_new_, data_.get_Sigma_stdv());
  concurrency_.average_and_compute_stddev(G_k_w_new, data_.get_G_k_w_stdv());

  // sum G4
  if (dca_iteration_ == parameters_.get_dca_iterations() - 1 && parameters_.isAccumulatingG4()) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t compute-error-bars on G4\t" << dca::util::print_time() << "\n\n";

    // Normalize and reduce the G4 channels one at a time. Copying the whole set of G4 tensors up
    // front would duplicate a multi-gigabyte buffer on every rank.
    // Note: the number and the order of the collective calls must not depend on the rank, which
    // holds here since G4_local.size() is identical on all ranks.
    const auto& G4_local = accumulator_.get_sign_times_G4();
    auto& G4_stddev = data_.get_G4_stdv();

    const double normalization =
        parameters_.get_beta() * parameters_.get_beta() * accumulator_.get_accumulated_sign();

    for (std::size_t channel = 0; channel < G4_local.size(); ++channel) {
      // Copy only the channel that is currently being reduced.
      auto G4_channel = G4_local[channel];
      G4_channel /= normalization;
      concurrency_.average_and_compute_stddev(G4_channel, G4_stddev[channel]);
    }
  }
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::collect_measurements() {
  auto collect_delayed = [&](auto& f) {
    if (compute_jack_knife_)
      concurrency_.leaveOneOutSum(f, true);
    else
      concurrency_.delayedSum(f);
  };

  const double local_time = total_time_;
  const bool accumulate_g4 =
      parameters_.isAccumulatingG4() && dca_iteration_ == parameters_.get_dca_iterations() - 1;

  {
    Profiler profiler("QMC-collectives", "CT-AUX solver", __LINE__);
    concurrency_.delayedSum(total_time_);
    concurrency_.delayedSum(accumulator_.get_Gflop());
    accumulated_sign_ = accumulator_.get_accumulated_sign();
    collect_delayed(accumulated_sign_);

    M_r_w_ = accumulator_.get_sign_times_M_r_w();
    collect_delayed(M_r_w_);

    if (accumulator_.compute_std_deviation()) {
      M_r_w_squared_ = accumulator_.get_sign_times_M_r_w_sqr();
      concurrency_.delayedSum(M_r_w_squared_);
    }

    if (parameters_.additional_time_measurements()) {
      Profiler profiler("Additional time measurements.", "QMC-collectives", __LINE__);
      concurrency_.delayedSum(accumulator_.get_G_r_t());
      concurrency_.delayedSum(accumulator_.get_G_r_t_stddev());
      concurrency_.delayedSum(accumulator_.get_charge_cluster_moment());
      concurrency_.delayedSum(accumulator_.get_magnetic_cluster_moment());
      concurrency_.delayedSum(accumulator_.get_densitymatrix_A77());
      concurrency_.delayedSum(accumulator_.get_densitymatrix_A89());
      concurrency_.delayedSum(accumulator_.get_spin_XX_chi());
      concurrency_.delayedSum(accumulator_.get_spin_ZZ_chi());
      concurrency_.delayedSum(accumulator_.get_charge_chi());
      concurrency_.delayedSum(accumulator_.get_dwave_pp_correlator());
    }

    // sum G4
    if (accumulate_g4) {
      for (int channel = 0; channel < data_.get_G4().size(); ++channel) {
        auto& G4 = data_.get_G4()[channel];
        G4 = accumulator_.get_sign_times_G4()[channel];
        if (compute_jack_knife_)
          concurrency_.leaveOneOutSum(G4);
        else
          concurrency_.localSum(G4, concurrency_.first());
      }
    }

    concurrency_.delayedSum(accumulator_.get_visited_expansion_order_k());
    concurrency_.delayedSum(accumulator_.get_error_distribution());

    concurrency_.resolveSums();

    // Average sign of the QMC run, reduced over all ranks. The numerator and the denominator are
    // plain sums: accumulated_sign_ above is a leave-one-out sum when a jackknife estimate is
    // requested, and the requested measurement count can differ from the number of measurements
    // that were actually done. Reducing the raw counters keeps log and output file identical.
    int accumulated_sign = accumulator_.get_accumulated_sign();
    int number_of_measurements = accumulator_.get_number_of_measurements();
    concurrency_.sum(accumulated_sign);
    concurrency_.sum(number_of_measurements);
    if (number_of_measurements > 0)
      average_sign_ = double(accumulated_sign) / double(number_of_measurements);
    else
      average_sign_ = 1.;
  }

  M_r_w_ /= accumulated_sign_;
  M_r_w_squared_ /= accumulated_sign_;
  if (accumulate_g4) {
    for (auto& G4 : data_.get_G4())
      G4 /= accumulated_sign_ * parameters_.get_beta() * parameters_.get_beta();
  }

  if (parameters_.additional_time_measurements()) {
    accumulator_.get_G_r_t() /= accumulated_sign_;
    data_.G_r_t = accumulator_.get_G_r_t();
    accumulator_.get_G_r_t_stddev() /= accumulated_sign_ * std::sqrt(parameters_.get_measurements());
    accumulator_.get_charge_cluster_moment() /= accumulated_sign_;
    accumulator_.get_magnetic_cluster_moment() /= accumulated_sign_;
    accumulator_.get_densitymatrix_A77() /= accumulated_sign_;
    accumulator_.get_densitymatrix_A89() /= accumulated_sign_;
    accumulator_.get_spin_XX_chi() /= accumulated_sign_;
    accumulator_.get_spin_ZZ_chi() /= accumulated_sign_;
    accumulator_.get_charge_chi() /= accumulated_sign_;
    accumulator_.get_dwave_pp_correlator() /= accumulated_sign_;
  }

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t Collect measurements \t" << dca::util::print_time() << "\n"
              << "\n\t\t\t QMC-local-time : " << local_time << " [sec]"
              << "\n\t\t\t QMC-total-time : " << total_time_ << " [sec]"
              << "\n\t\t\t Gflop   : " << accumulator_.get_Gflop() << " [Gf]"
              << "\n\t\t\t Gflop/s   : " << accumulator_.get_Gflop() / local_time << " [Gf/s]"
              << "\n\t\t\t sign     : " << average_sign_ << " \n";

  averaged_ = true;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::symmetrize_measurements() {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t symmetrize measurements has started \t" << dca::util::print_time() << "\n";

  symmetrize::execute(M_r_w_, data_.H_symmetry);
  symmetrize::execute(M_r_w_squared_, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::compute_G_k_w_from_M_r_w() {
  func::function<std::complex<double>, NuNuKClusterWDmn> M_k_w;
  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(M_r_w_, M_k_w);

  int matrix_size = b::dmn_size() * s::dmn_size() * b::dmn_size() * s::dmn_size();
  int matrix_dim = b::dmn_size() * s::dmn_size();

  std::complex<double>* G_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* G0_cluster_excluded_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* M_matrix = new std::complex<double>[matrix_size];
  std::complex<double>* G0_times_M_matrix = new std::complex<double>[matrix_size];

  // G = G0 - G0*M*G0/beta

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      memset(static_cast<void*>(G_matrix), 0, sizeof(std::complex<double>) * matrix_size);
      memset(static_cast<void*>(G0_times_M_matrix), 0, sizeof(std::complex<double>) * matrix_size);

      memcpy(G0_cluster_excluded_matrix, &data_.G0_k_w_cluster_excluded(0, 0, 0, 0, k_ind, w_ind),
             sizeof(std::complex<double>) * matrix_size);
      memcpy(M_matrix, &M_k_w(0, 0, 0, 0, k_ind, w_ind), sizeof(std::complex<double>) * matrix_size);

      // G0 * M --> G0_times_M_matrix
      dca::linalg::blas::gemm("N", "N", matrix_dim, matrix_dim, matrix_dim, 1.,
                              G0_cluster_excluded_matrix, matrix_dim, M_matrix, matrix_dim, 0.,
                              G0_times_M_matrix, matrix_dim);

      // - G0_times_M_matrix * G0 / beta --> G_matrix
      dca::linalg::blas::gemm("N", "N", matrix_dim, matrix_dim, matrix_dim,
                              -1. / parameters_.get_beta(), G0_times_M_matrix, matrix_dim,
                              G0_cluster_excluded_matrix, matrix_dim, 0., G_matrix, matrix_dim);

      // G_matrix + G0_cluster_excluded_matrix --> G_matrix
      for (int l = 0; l < matrix_size; l++)
        G_matrix[l] = G_matrix[l] + G0_cluster_excluded_matrix[l];

      memcpy(&data_.G_k_w(0, 0, 0, 0, k_ind, w_ind), G_matrix,
             sizeof(std::complex<double>) * matrix_size);
    }
  }

  symmetrize::execute(data_.G_k_w, data_.H_symmetry);

  delete[] G_matrix;
  delete[] G0_cluster_excluded_matrix;
  delete[] M_matrix;
  delete[] G0_times_M_matrix;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
double CtauxClusterSolver<device_t, Parameters, Data>::compute_S_k_w_from_G_k_w() {
  const double alpha = parameters_.get_self_energy_mixing_factor();
  //     double L2_difference_norm = 0;
  //     double L2_Sigma_norm      = 0;

  int matrix_dim = b::dmn_size() * s::dmn_size();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_inverted_matrix(matrix_dim);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_cluster_excluded_inverted_matrix(
      matrix_dim);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> sigma_matrix(matrix_dim);

  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

  // Sigma = 1/G0 - 1/G

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      dca::linalg::matrixop::copyArrayToMatrix(matrix_dim, matrix_dim,
                                               &data_.G_k_w(0, 0, 0, 0, k_ind, w_ind), matrix_dim,
                                               G_inverted_matrix);
      dca::linalg::matrixop::inverse(G_inverted_matrix, ipiv, work);

      dca::linalg::matrixop::copyArrayToMatrix(
          matrix_dim, matrix_dim, &data_.G0_k_w_cluster_excluded(0, 0, 0, 0, k_ind, w_ind),
          matrix_dim, G0_cluster_excluded_inverted_matrix);
      dca::linalg::matrixop::inverse(G0_cluster_excluded_inverted_matrix, ipiv, work);

      for (int j = 0; j < sigma_matrix.nrCols(); ++j)
        for (int i = 0; i < sigma_matrix.nrRows(); ++i)
          sigma_matrix(i, j) = G0_cluster_excluded_inverted_matrix(i, j) - G_inverted_matrix(i, j);

      dca::linalg::matrixop::copyMatrixToArray(sigma_matrix, &data_.Sigma(0, 0, 0, 0, k_ind, w_ind),
                                               matrix_dim);
    }
  }

  // Compute error on G and Self Energy.
  if (compute_jack_knife_) {
    data_.get_G_k_w_error() = concurrency_.jackknifeError(data_.G_k_w, true);
    data_.get_G_r_w_error() = concurrency_.jackknifeError(data_.G_r_w, true);
    data_.get_Sigma_error() = concurrency_.jackknifeError(data_.Sigma, true);
  }

  // set_non_interacting_bands_to_zero();

  symmetrize::execute(data_.Sigma, data_.H_symmetry);

  if (parameters_.adjust_self_energy_for_double_counting())
    adjust_self_energy_for_double_counting();

  // Keep the un-mixed QMC self-energy. finalize() derives Sigma_zero_moment and
  // standard_deviation from it; earlier those diagnostics read data_.Sigma after the mixing and
  // therefore described the extrapolated self-energy. Sigma_new_ used to be a member that was
  // written to the output but never assigned.
  Sigma_new_ = data_.Sigma;

  double L2_norm = mix_self_energy(alpha);

  return L2_norm;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::compute_G_k_w_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& M_k_w_new,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new) const {
  //     if(concurrency_.id()==0)
  //       std::cout << "\n\t\t compute-G_k_w_new\t" << dca::util::print_time() << "\n\n";

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix("G_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix("G0_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> M_matrix("M_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_M_matrix("M_matrix", nu::dmn_size());

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = data_.G0_k_w_cluster_excluded(i, j, k_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          M_matrix(i, j) = M_k_w_new(i, j, k_ind, w_ind);

      dca::linalg::matrixop::gemm(G0_matrix, M_matrix, G0_M_matrix);
      dca::linalg::matrixop::gemm(G0_M_matrix, G0_matrix, G_matrix);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_k_w_new(i, j, k_ind, w_ind) = G0_matrix(i, j) - G_matrix(i, j) / parameters_.get_beta();
    }
  }

  symmetrize::execute(G_k_w_new, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::compute_S_k_w_new(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& G_k_w_new,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>>& S_k_w_new) {
  //     if(concurrency_.id()==0)
  //       std::cout << "\n\t\t start compute-S_k_w\t" << dca::util::print_time() << "\n\n";

  int N = nu::dmn_size();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix(N);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix(N);

  for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++) {
    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = data_.G0_k_w_cluster_excluded(i, j, k_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_matrix(i, j) = G_k_w_new(i, j, k_ind, w_ind);

      dca::linalg::matrixop::inverse(G_matrix);
      dca::linalg::matrixop::inverse(G0_matrix);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          S_k_w_new(i, j, k_ind, w_ind) = G0_matrix(i, j) - G_matrix(i, j);
    }
  }

  if (parameters_.adjust_self_energy_for_double_counting())
    adjust_self_energy_for_double_counting();

  //     if(concurrency_.id()==0)
  //       std::cout << "\n\t\t end compute-S_k_w\t" << dca::util::print_time() << "\n\n";

  symmetrize::execute(S_k_w_new, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::set_non_interacting_bands_to_zero() {
  //  for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++){
  //      for(int l2=0; l2<b::dmn_size(); l2++){
  //        for(int l1=0; l1<b::dmn_size(); l1++){
  //
  //          if( !(parameters_.is_interacting_band()[l1] and
  //                parameters_.is_interacting_band()[l2]))
  //            {
  //              data_.Sigma(l1,0,l2,0,k_ind, w_ind) = 0.;
  //              data_.Sigma(l1,1,l2,1,k_ind, w_ind) = 0.;
  //            }
  //        }
  //      }
  //    }
  //  }
  //
  //  for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //      for(int l2=0; l2<2*b::dmn_size(); l2++)
  //        for(int l1=0; l1<2*b::dmn_size(); l1++)
  //          if( !(l1==l2 and parameters_.is_interacting_band()[l1]) )
  //            data_.Sigma(l1,l2,k_ind, w_ind) = 0.;
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
void CtauxClusterSolver<device_t, Parameters, Data>::adjust_self_energy_for_double_counting() {
  set_non_interacting_bands_to_zero();

  //  func::function<double, nu> d_0;
  //  for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //      for(int w_ind=0; w_ind<32; w_ind++)
  //        d_0(l1) += real(data_.Sigma(l1,l1,k_ind,w_ind));
  //
  //  d_0 /= double(32.*KClusterDmn::dmn_size());
  //
  //  for(int l1=0; l1<b::dmn_size()*s::dmn_size(); l1++)
  //    for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //        data_.Sigma(l1,l1,k_ind,w_ind) -= d_0(l1);
  //
  //  if(parameters_.get_double_counting_method()=="constant")
  //    {
  //      std::vector<int>& interacting_bands = parameters_.get_interacting_orbitals();
  //
  //      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
  //        for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++)
  //          for(int s_ind=0; s_ind<s::dmn_size(); s_ind++)
  //            for(int b_ind=0; b_ind<interacting_bands.size(); b_ind++)
  //              data_.Sigma(interacting_bands[b_ind], s_ind,
  //                         interacting_bands[b_ind], s_ind,
  //                         k_ind                   , w_ind) -=
  //  parameters_.get_double_counting_correction();
  //    }
  //
  //  if(parameters_.get_double_counting_method()=="adaptive")
  //    {
  //      std::vector<int>& interacting_bands = parameters_.get_interacting_orbitals();
  //
  //      for(int b_ind=0; b_ind<interacting_bands.size(); b_ind++)
  //        for(int k_ind=0; k_ind<KClusterDmn::dmn_size(); k_ind++){
  //          for(int s_ind=0; s_ind<s::dmn_size(); s_ind++){
  //
  //            double value = real(data_.Sigma(interacting_bands[b_ind], s_ind,
  //                                           interacting_bands[b_ind], s_ind,
  //                                           k_ind                   , 0));
  //
  //            for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
  //
  //              data_.Sigma(interacting_bands[b_ind], s_ind,
  //                         interacting_bands[b_ind], s_ind,
  //                         k_ind                   , w_ind) -= value;
  //            }
  //          }
  //        }
  //    }

  symmetrize::execute(data_.Sigma, data_.H_symmetry);
}

template <dca::linalg::DeviceType device_t, class Parameters, class Data>
double CtauxClusterSolver<device_t, Parameters, Data>::mix_self_energy(double alpha) {
  symmetrize::execute(data_.Sigma, data_.H_symmetry);
  symmetrize::execute(data_.Sigma_cluster, data_.H_symmetry);

  // The convergence criterion of the fixed-point iteration Sigma_cg -> Sigma_QMC is evaluated on
  // the un-mixed quantities, so that its meaning does not depend on the mixing scheme:
  //   now      L2 = || Sigma_QMC - Sigma_cg ||, i.e. the true residual r of the iteration.
  //   earlier  L2 = || Sigma_mixed - Sigma_cg ||, where Sigma_mixed is the output of the mixer.
  //            For LINEAR this is alpha * || r ||, and for ANDERSON it additionally contains the
  //            extrapolation step sum_i c_i (x_i - x_m). The two schemes were therefore not
  //            comparable, and ANDERSON could be reported as converged while r was still large.
  // Consequently, with linear mixing the value reported now is a factor 1 / alpha larger than
  // before, so the DCA::accuracy of an existing input file has to be re-calibrated.
  const double L2_error = self_energy_residual_l2_norm();
  const double error_infty_norm = self_energy_residual_infty_norm();

  // Keep the self-energy with the smallest residual seen so far. L2_error describes
  // data_.Sigma_cluster, i.e. the output of the previous iteration, so storing it together with
  // this value gives a consistent (Sigma, L2) pair to fall back to.
  //
  // best_Sigma_ is overwritten with the mixed self-energy further down, once the mixer has run:
  // that is the self-energy the DCA loop would have written for this iteration, i.e. the one that
  // the density and the chemical potential of this iteration are evaluated with. The input
  // self-energy is only kept as a fallback for the iteration that stops the loop, because the
  // mixing is skipped there.
  const bool new_best = L2_error < best_L2_error_;
  if (new_best) {
    best_L2_error_ = L2_error;
    best_dca_iteration_ = dca_iteration_;
    best_Sigma_ = data_.Sigma_cluster;
    best_Sigma_new_ = Sigma_new_;
    best_Sigma_stdv_ = data_.get_Sigma_stdv();
    best_G_k_w_stdv_ = data_.get_G_k_w_stdv();
    best_G_k_w_ = data_.G_k_w;
    best_Sigma_old_ = Sigma_old_;
    best_M_r_w_ = M_r_w_;
    best_M_r_w_squared_ = M_r_w_squared_;
    best_chemical_potential_ = chemical_potential_at_start_;
    best_sign_ = current_sign_;
    best_expansion_order_ = current_expansion_order_;
    best_has_error_bars_ = compute_jack_knife_;
    if (best_has_error_bars_)
      best_G_r_w_error_ = data_.get_G_r_w_error();
  }

  // Stop as soon as the residual grows again. This is the only criterion and it applies to both
  // mixing schemes: whatever caused the growth, the previous iteration had the smaller residual,
  // so continuing would only move away from the fixed point. The comparison is made against the
  // previous iteration, which is the only linear reference available: evaluating both variants
  // would require an additional QMC run. The criterion has no tolerance: the first increase stops
  // the loop.
  const bool regression = previous_L2_error_ > 0. && L2_error > previous_L2_error_;

  if (regression) {
    stop_requested_ = true;
    // Discard the step of this iteration: the output must be the best iteration, not the one that
    // made L2 grow. best_Sigma_ is the mixed self-energy of the saved iteration and carries the
    // double-counting correction like every self-energy of the loop, so
    // adjust_impurity_self_energy() removes it again as usual.
    data_.Sigma = best_Sigma_;
    // Restore the chemical potential of that iteration as well. best_dca_iteration_ can be older
    // than the previous iteration, so taking the value from the DCA loop would be wrong.
    parameters_.get_chemical_potential() = best_chemical_potential_;
    // Restore the un-mixed quantities of the same iteration, so that the diagnostics derived in
    // finalize() and the error bars written to the output belong to the saved iteration too.
    Sigma_new_ = best_Sigma_new_;
    data_.get_Sigma_stdv() = best_Sigma_stdv_;
    data_.get_G_k_w_stdv() = best_G_k_w_stdv_;
    // G_k_w and G_r_w are measured again in every iteration. Restore the ones that belong to
    // best_Sigma_ and redo the Fourier transform, so that G, G( tau ) and the density derived from
    // them all describe the saved iteration.
    data_.G_k_w = best_G_k_w_;
    math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(data_.G_k_w, data_.G_r_w);
    Sigma_old_ = best_Sigma_old_;
    M_r_w_ = best_M_r_w_;
    M_r_w_squared_ = best_M_r_w_squared_;
    if (best_has_error_bars_)
      data_.get_G_r_w_error() = best_G_r_w_error_;
    else if (compute_jack_knife_)
      // A jackknife estimate only exists for the last DCA iteration, so the saved iteration usually
      // has none. Leaving the array untouched would write the error bars of the rejected run
      // together with the restored self-energy, so clear them instead.
      data_.get_G_r_w_error() = 0.;
    if (concurrency_.id() == concurrency_.first()) {
      std::cout << "\n\t\t iteration " << (dca_iteration_ + 1) << ": stopping, L2 increased from "
                << previous_L2_error_ << " to " << L2_error << "\n\t\t restoring iteration "
                << (best_dca_iteration_ + 1) << " as the final result (L2 = " << best_L2_error_
                << ", sign = " << best_sign_ << ", <k> = " << best_expansion_order_
                << ") and stopping the DCA loop\n";
    }
    return best_L2_error_;
  }
  previous_L2_error_ = L2_error;

  // alpha == 0 freezes the self-energy: the iteration measures at the current Sigma and does not
  // update it. That is what a fixed-self-energy measurement needs, e.g. a one-iteration four-point
  // run on top of a converged sp file. It cannot converge an sp run, so it is reported explicitly.
  // Both schemes are skipped here; with ANDERSON the mixer would otherwise still extrapolate the
  // history of inputs, which is not the same as "do not update".
  if (alpha > 0.) {
    if (parameters_.get_self_energy_mixing_scheme() == SelfEnergyMixingScheme::LINEAR) {
      for (int l = 0; l < data_.Sigma.size(); l++)
        data_.Sigma(l) = alpha * data_.Sigma(l) + (1. - alpha) * data_.Sigma_cluster(l);
    }
    else {
      // Anderson/Pulay (DIIS) acceleration. data_.Sigma holds Sigma_QMC and data_.Sigma_cluster
      // holds the input self-energy; the mixer is in-place safe.
      //
      // Parallel contract: keep this outside of any rank guard. data_.Sigma and data_.Sigma_cluster
      // are globally replicated (M_r_w_ and the coarsegraining are reduced with MPI_Allreduce), so
      // executing the mixer on every rank keeps the replicated histories in sync. Guarding this
      // call by rank, or skipping it on some ranks, would silently make the ranks diverge.
      self_energy_mixer_.mix(data_.Sigma.values(), data_.Sigma_cluster.values(),
                             data_.Sigma.values(), static_cast<int>(data_.Sigma.size()));

      // Report the DIIS coefficients on the first rank only. The mixer itself runs on every rank
      // (see the parallel contract above), but printing must be rank-guarded.
      if (concurrency_.id() == concurrency_.first()) {
        const auto& c = self_energy_mixer_.lastCoefficients();
        std::cout << "\n\t\t Anderson mixing: m = " << self_energy_mixer_.historySize()
                  << ", |r| = " << self_energy_mixer_.lastResidualNorm()
                  << ", |sum c_i r_i| = " << self_energy_mixer_.lastCombinedResidualNorm();
        // The extrapolation only starts once the history is full, so a linear step with fewer
        // entries is the warm-up of a fresh history, not a failed extrapolation.
        if (self_energy_mixer_.lastFellBack()) {
          const int history_size = static_cast<int>(self_energy_mixer_.historySize());
          std::cout << (history_size < parameters_.get_self_energy_mixing_history_size()
                            ? " (warm-up, linear step)"
                            : " (fallback, linear step)");
        }
        if (self_energy_mixer_.lastCoefficientsClipped())
          std::cout << " (coefficients clipped)";
        std::cout << "\n\t\t c = [";
        for (std::size_t j = 0; j < c.size(); ++j)
          std::cout << (j == 0 ? "" : ", ") << c[j];
        std::cout << "]" << std::endl;
      }
    }
  }
  else {
    // alpha == 0: keep the input self-energy, i.e. do not update it. data_.Sigma holds Sigma_QMC
    // at this point, so it has to be reset explicitly; leaving it alone would mean "take the full
    // QMC step" instead of "do not update".
    data_.Sigma = data_.Sigma_cluster;
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t self-energy mixing factor is 0: the self-energy is not updated\n";
  }

  // Only now is the self-energy of this iteration final: the mixer output is what the DCA loop
  // writes and what the density and the chemical potential of this iteration are evaluated with.
  // Saving the input self-energy instead would put a state into the output whose density differs
  // from the recorded one by the size of a mixing step.
  if (new_best)
    best_Sigma_ = data_.Sigma;

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\n\t\t |Sigma_QMC - Sigma_cg|_infty ~ " << error_infty_norm;
    std::cout << "\n\t\t |Sigma_QMC - Sigma_cg|_2 ~ " << L2_error << "\n\n";
  }
  return L2_error;
}

// L2 norm of the fixed-point residual Sigma_QMC - Sigma_cg, restricted to the intra-band diagonal
// at the lowest Matsubara frequency. This is the quantity the DCA::accuracy criterion is based on.
template <dca::linalg::DeviceType device_t, class Parameters, class Data>
double CtauxClusterSolver<device_t, Parameters, Data>::self_energy_residual_l2_norm() const {
  const int offset = std::min(1, w::dmn_size() / 2);

  double diff_L2_norm = 0.;
  for (int w_ind = w::dmn_size() / 2; w_ind < w::dmn_size() / 2 + offset; w_ind++)
    for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++)
      for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
        diff_L2_norm += std::norm(data_.Sigma(l1, l1, k_ind, w_ind) -
                                  data_.Sigma_cluster(l1, l1, k_ind, w_ind));

  return std::sqrt(diff_L2_norm) / double(KClusterDmn::dmn_size());
}

// Infinity norm of the fixed-point residual Sigma_QMC - Sigma_cg over the ten lowest Matsubara
// frequencies.
template <dca::linalg::DeviceType device_t, class Parameters, class Data>
double CtauxClusterSolver<device_t, Parameters, Data>::self_energy_residual_infty_norm() const {
  const int offset = std::min(10, w::dmn_size() / 2);

  double error_infty_norm = 0.;
  for (int w_ind = w::dmn_size() / 2; w_ind < w::dmn_size() / 2 + offset; w_ind++)
    for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); k_ind++)
      for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
        error_infty_norm = std::max(error_infty_norm, std::abs(data_.Sigma(l1, l1, k_ind, w_ind) -
                                                               data_.Sigma_cluster(l1, l1, k_ind,
                                                                                   w_ind)));
  return error_infty_norm;
}



template <dca::linalg::DeviceType device_t, class Parameters, class Data>
auto CtauxClusterSolver<device_t, Parameters, Data>::local_G_k_w() const {
  if (averaged_)
    throw std::logic_error("The local data was already averaged.");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w_new(
      "G_k_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> M_k_w_new(
      "M_k_w_new");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> M_r_w_new(
      accumulator_.get_sign_times_M_r_w(), "M_r_w_new");

  M_r_w_new /= accumulator_.get_accumulated_sign();

  math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(M_r_w_new, M_k_w_new);

  compute_G_k_w_new(M_k_w_new, G_k_w_new);

  return G_k_w_new;
}

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_CTAUX_CLUSTER_SOLVER_HPP
