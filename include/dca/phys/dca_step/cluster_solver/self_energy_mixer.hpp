// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Anderson/Pulay (DIIS) mixer for the cluster self-energy self-consistency loop.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SELF_ENERGY_MIXER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SELF_ENERGY_MIXER_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <deque>
#include <stdexcept>
#include <utility>
#include <vector>

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

// Anderson/Pulay (DIIS) mixer for the self-energy fixed-point iteration x -> F(x), with residual
// r = F(x) - x:
//   x_new = sum_i c_i (x_i + alpha * r_i),
// where the coefficients c_i minimize |sum_i c_i r_i|^2 subject to sum_i c_i = 1.
// For a single history entry this reduces to plain linear mixing, x + alpha * r.
//
// Inner product: the residual norm can be weighted element by element (setResidualWeights).  The
// self-energy spectrum is dominated by the high-frequency tail, where consecutive iterations are
// indistinguishable and where Sigma = G0^-1 - G^-1 suffers from cancellation.  An unweighted norm
// therefore fixes the DIIS coefficients mostly from the least interesting part of the spectrum.
// The weights restrict the minimization to the low-frequency window that carries the convergence
// information; the mixed self-energy always uses all elements.
//
// Robustness:
//   - the bordered Gram system is Tikhonov regularized and solved with a *relative* pivot threshold,
//     so the criterion does not depend on the absolute magnitude of the residual;
//   - the coefficients are clipped to +-max_coefficient and renormalized to sum_i c_i = 1;
//   - every fallback (single history, ill-conditioned or non-finite solve, clipping that cannot be
//     renormalized) degrades gracefully to the linear-mixing step;
//   - the history is dropped when the residual grows by more than residual_growth_limit, so that
//     one bad iteration (or a parameter change in the outer loop) does not poison the next
//     max_history iterations;
//   - non-finite self-energies are rejected instead of being stored in the history.
//
// Difference to the earlier version of this mixer (the prototype that had no safeguards):
//   norm     earlier: A_ij = r_i . r_j summed over all elements with uniform weight. The
//                     high-frequency tail, whose error grows with w, then dominated the solve.
//            now:     A_ij = sum_k w_k r_i(k) r_j(k), see setResidualWeights().
//   pivot    earlier: |pivot| < 1.e-14, an absolute threshold unrelated to the residual scale.
//            now:     |pivot| < pivot_tolerance * max_diag.
//   coeffs   earlier: used exactly as solved, so O(1.e2) values with alternating signs occurred.
//            now:     clipped to +-max_coefficient and renormalized to sum_i c_i = 1.
//   history  earlier: kept for max_history iterations no matter how the residual evolved.
//            now:     dropped once the residual grows by more than residual_growth_limit.
//   input    earlier: a non-finite Sigma was stored and spoiled the next max_history iterations.
//            now:     rejected with an exception.
//   tuning   earlier: alpha and max_history only; the regularization was a hard-coded 1.e-8.
//            now:     Config exposes regularization and max_coefficient as well.
// The algorithm itself (the constrained least-squares problem and the update x_new above) is
// unchanged, and so is the meaning of alpha and of max_history.
//
// Parallel contract: this is a purely local, stateful object and performs no communication. It must
// be called exactly once per DCA iteration on every MPI rank and with identical inputs, so that the
// replicated histories stay in sync. Putting a call behind a rank guard, or calling it a
// rank-dependent number of times, silently makes the ranks diverge.
class SelfEnergyMixer {
public:
  struct Config {
    // Mixing factor alpha. alpha = 1 is pure Anderson, alpha < 1 damps the step.
    double alpha = 1.;
    // Number of history entries. A value < 1 marks an unconfigured mixer.
    int max_history = 0;
    // Relative Tikhonov regularization of the Gram matrix. The earlier version hard-coded 1.e-8;
    // the default keeps that value, but it is now tunable. Larger values pull the coefficients
    // towards the uniform average, which is more robust but less accelerative.
    double regularization = 1.e-8;
    // Clipping bound for the DIIS coefficients. The earlier version had no bound at all; using a
    // very large value (e.g. 1.e30) reproduces that behavior. Must be >= 1.
    double max_coefficient = 10.;
    // Drop the history when the residual grows by more than this factor. The earlier version never
    // dropped the history. Must be > 1; use a very large value to disable the check.
    double residual_growth_limit = 4.;
    // Pivot threshold of the linear solve, relative to the largest Gram diagonal. The earlier
    // version used the absolute threshold 1.e-14, which fails once the residual becomes small.
    double pivot_tolerance = 1.e-12;
  };

  SelfEnergyMixer() = default;

  void configure(const Config& config) {
    validateConfig(config);
    config_ = config;
  }

  void configure(const double alpha, const int max_history) {
    Config config = config_;
    config.alpha = alpha;
    config.max_history = max_history;
    validateConfig(config);
    config_ = config;
  }

  // Uses the weights of the most recent mix() call for the diagnostics.
  double getAlpha() const {
    return config_.alpha;
  }

  int getMaxHistory() const {
    return config_.max_history;
  }

  // Weights of the inner product in which the DIIS residual is minimized. Must have one entry per
  // element of the self-energy. An empty vector restores the uniform norm.
  void setResidualWeights(std::vector<double> weights) {
    for (std::size_t i = 0; i < weights.size(); ++i)
      if (!(weights[i] >= 0.) || !std::isfinite(weights[i]))
        throw std::logic_error("SelfEnergyMixer: the residual weights must be non-negative.");
    weights_ = std::move(weights);
    if (!weights_.empty() && size_ != 0 && static_cast<int>(weights_.size()) != size_)
      throw std::logic_error(
          "SelfEnergyMixer: the number of residual weights must match the size of the "
          "self-energy.");
  }

  void clearResidualWeights() {
    weights_.clear();
  }

  void reset() {
    history_.clear();
    x_new_.clear();
    size_ = 0;
    last_coefficients_.clear();
    last_fell_back_ = false;
    last_clipped_ = false;
    last_residual_norm_ = 0.;
    last_combined_residual_norm_ = 0.;
    previous_residual_norm_ = 0.;
    restarts_ = 0;
  }

  std::size_t historySize() const {
    return history_.size();
  }

  // Coefficients used by the most recent mix() call: the DIIS weights c_j for a regular Anderson
  // step, or {1} for the single-entry / ill-conditioned linear-mixing fallback.
  const std::vector<double>& lastCoefficients() const {
    return last_coefficients_;
  }

  // True if the most recent mix() call fell back to linear mixing (history too short, the DIIS
  // system was ill-conditioned, or the coefficients could not be normalized).
  bool lastFellBack() const {
    return last_fell_back_;
  }

  // True if the coefficients of the most recent mix() call had to be clipped.
  bool lastCoefficientsClipped() const {
    return last_clipped_;
  }

  // Weighted norm of the residual r = F(x) - x of the most recent mix() call, i.e. the quantity a
  // scheme-independent convergence criterion should be based on.
  double lastResidualNorm() const {
    return last_residual_norm_;
  }

  // Weighted norm of the DIIS residual sum_i c_i r_i, i.e. what the coefficient solve minimized.
  double lastCombinedResidualNorm() const {
    return last_combined_residual_norm_;
  }

  // Number of times the history was dropped because the residual grew too fast.
  int restarts() const {
    return restarts_;
  }

  // In-place safe: sigma_out may alias sigma_qmc and/or sigma_in.
  void mix(const std::complex<double>* const sigma_qmc,
           const std::complex<double>* const sigma_in, std::complex<double>* const sigma_out,
           const int size) {
    checkConfiguration();
    if (size <= 0)
      throw std::logic_error("SelfEnergyMixer: the size of the self-energy must be positive.");
    if (sigma_qmc == nullptr || sigma_in == nullptr || sigma_out == nullptr)
      throw std::logic_error("SelfEnergyMixer: the self-energy pointers must not be null.");
    // Guard against a rank-dependent or iteration-dependent size, which would make the stored
    // history inconsistent with the vector that is being mixed.
    if (size_ == 0)
      size_ = size;
    else if (size_ != size)
      throw std::logic_error("SelfEnergyMixer: the size of the self-energy must not change.");
    if (!weights_.empty() && static_cast<int>(weights_.size()) != size)
      throw std::logic_error(
          "SelfEnergyMixer: the number of residual weights must match the size of the "
          "self-energy.");

    // Residual r = F(x) - x. Reject a broken QMC result instead of storing it for max_history
    // iterations.
    std::vector<std::complex<double>> residual(static_cast<std::size_t>(size));
    for (int i = 0; i < size; ++i) {
      residual[i] = sigma_qmc[i] - sigma_in[i];
      if (!std::isfinite(residual[i].real()) || !std::isfinite(residual[i].imag()))
        throw std::runtime_error(
            "SelfEnergyMixer: the self-energy contains non-finite values. This usually means that "
            "the QMC calculation has failed.");
    }

    // A residual that grew a lot means the fixed-point iteration has moved somewhere else (bad
    // iteration, chemical potential update, double counting update, ...). Extrapolating from the
    // old history would be worse than starting over.
    const double residual_norm = weightedNorm(residual);
    if (!history_.empty() && previous_residual_norm_ > 0. &&
        residual_norm > config_.residual_growth_limit * previous_residual_norm_) {
      history_.clear();
      ++restarts_;
    }

    history_.push_back(Entry{std::vector<std::complex<double>>(sigma_in, sigma_in + size),
                             std::move(residual)});
    while (static_cast<int>(history_.size()) > config_.max_history)
      history_.pop_front();

    // The most recent residual, needed for the linear-mixing step.
    const std::vector<std::complex<double>>& latest_residual = history_.back().r;

    bool use_diis = false;
    std::vector<double> coefficients;
    last_clipped_ = false;
    if (history_.size() > 1 && solveCoefficients(coefficients))
      use_diis = clipCoefficients(coefficients);

    x_new_.resize(static_cast<std::size_t>(size));

    if (use_diis) {
      std::fill(x_new_.begin(), x_new_.end(), std::complex<double>(0., 0.));
      const int m = static_cast<int>(history_.size());
      for (int j = 0; j < m; ++j)
        for (int i = 0; i < size; ++i)
          x_new_[i] += coefficients[j] * (history_[j].x[i] + config_.alpha * history_[j].r[i]);
      last_coefficients_ = coefficients;
      last_fell_back_ = false;
      last_combined_residual_norm_ = combinedResidualNorm(coefficients);
    }
    else {
      for (int i = 0; i < size; ++i)
        x_new_[i] = sigma_in[i] + config_.alpha * latest_residual[i];
      last_coefficients_ = std::vector<double>{1.};
      last_fell_back_ = true;
      last_combined_residual_norm_ = config_.alpha * residual_norm;
    }

    std::copy(x_new_.begin(), x_new_.end(), sigma_out);

    last_residual_norm_ = residual_norm;
    previous_residual_norm_ = residual_norm;
  }

private:
  struct Entry {
    std::vector<std::complex<double>> x;
    std::vector<std::complex<double>> r;
  };

  static void validateConfig(const Config& config) {
    if (config.max_history < 1)
      throw std::logic_error("SelfEnergyMixer: max_history must be >= 1.");
    if (!(config.alpha > 0.) || !std::isfinite(config.alpha))
      throw std::logic_error("SelfEnergyMixer: alpha must be positive.");
    if (!(config.regularization >= 0.) || !std::isfinite(config.regularization))
      throw std::logic_error("SelfEnergyMixer: the regularization must be non-negative.");
    if (!(config.max_coefficient >= 1.) || !std::isfinite(config.max_coefficient))
      throw std::logic_error("SelfEnergyMixer: max_coefficient must be >= 1.");
    if (!(config.residual_growth_limit > 1.) || !std::isfinite(config.residual_growth_limit))
      throw std::logic_error("SelfEnergyMixer: residual_growth_limit must be > 1.");
    if (!(config.pivot_tolerance > 0.) || !std::isfinite(config.pivot_tolerance))
      throw std::logic_error("SelfEnergyMixer: pivot_tolerance must be positive.");
  }

  void checkConfiguration() const {
    if (config_.max_history < 1)
      throw std::logic_error("SelfEnergyMixer: configure() must be called before mix().");
  }

  double weight(const std::size_t i) const {
    return weights_.empty() ? 1. : weights_[i];
  }

  // Real inner product of two residuals, treated as real vectors of length 2 * size.
  double weightedInnerProduct(const std::vector<std::complex<double>>& a,
                              const std::vector<std::complex<double>>& b) const {
    const std::size_t size = a.size();
    double dot = 0.;
    for (std::size_t i = 0; i < size; ++i)
      dot += weight(i) * (a[i].real() * b[i].real() + a[i].imag() * b[i].imag());
    return dot;
  }

  double weightedNorm(const std::vector<std::complex<double>>& a) const {
    return std::sqrt(std::max(0., weightedInnerProduct(a, a)));
  }

  double combinedResidualNorm(const std::vector<double>& coefficients) const {
    const std::size_t size = history_.back().r.size();
    double norm_squared = 0.;
    for (std::size_t i = 0; i < size; ++i) {
      std::complex<double> sum(0., 0.);
      for (std::size_t j = 0; j < coefficients.size(); ++j)
        sum += coefficients[j] * history_[j].r[i];
      norm_squared += weight(i) * std::norm(sum);
    }
    return std::sqrt(norm_squared);
  }

  // Solves the DIIS linear system for the coefficients c:
  //   min c^T G c  subject to  sum_i c_i = 1,
  // where G is the (weighted) Gram matrix of the residuals. Returns false if the system is
  // singular or ill-conditioned.
  bool solveCoefficients(std::vector<double>& coefficients) {
    const int m = static_cast<int>(history_.size());
    const std::size_t n = static_cast<std::size_t>(m) + 1;

    std::vector<double> a(n * n, 0.);
    double max_diag = 0.;

    for (int i = 0; i < m; ++i) {
      for (int j = i; j < m; ++j) {
        const double dot = weightedInnerProduct(history_[i].r, history_[j].r);
        a[i * n + j] = dot;
        a[j * n + i] = dot;
      }
      max_diag = std::max(max_diag, std::abs(a[i * n + i]));
    }

    // A vanishing (or non-finite) residual leaves the coefficients undetermined.
    if (!(max_diag > 0.) || !std::isfinite(max_diag))
      return false;

    // Tikhonov regularization, then the Lagrange constraint sum_i c_i = 1.
    const double lambda = config_.regularization * max_diag;
    for (int i = 0; i < m; ++i)
      a[i * n + i] += lambda;
    for (int i = 0; i < m; ++i) {
      a[i * n + m] = 1.;
      a[m * n + i] = 1.;
    }
    a[m * n + m] = 0.;

    std::vector<double> b(n, 0.);
    b[m] = 1.;

    std::vector<double> solution;
    if (!solveLinearSystem(a, b, solution, config_.pivot_tolerance * max_diag))
      return false;

    coefficients.assign(solution.begin(), solution.begin() + m);
    return true;
  }

  // Clips the coefficients to +-max_coefficient and renormalizes them to sum_i c_i = 1.
  // Returns false if the clipped coefficients cannot be renormalized.
  bool clipCoefficients(std::vector<double>& coefficients) {
    double sum = 0.;
    bool clipped = false;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
      double& c = coefficients[i];
      if (!std::isfinite(c))
        return false;
      if (c > config_.max_coefficient) {
        c = config_.max_coefficient;
        clipped = true;
      }
      else if (c < -config_.max_coefficient) {
        c = -config_.max_coefficient;
        clipped = true;
      }
      sum += c;
    }
    if (!(std::abs(sum) > 1.e-12))
      return false;
    for (std::size_t i = 0; i < coefficients.size(); ++i)
      coefficients[i] /= sum;
    last_clipped_ = clipped;
    return true;
  }

  // Dense Gaussian elimination with partial pivoting. Returns false if the system is singular
  // (|pivot| < pivot_threshold) or if the solution is not finite.
  static bool solveLinearSystem(std::vector<double> a, std::vector<double> b,
                                std::vector<double>& x, const double pivot_threshold) {
    const int n = static_cast<int>(b.size());
    for (int col = 0; col < n; ++col) {
      int pivot = col;
      for (int row = col + 1; row < n; ++row)
        if (std::abs(a[row * n + col]) > std::abs(a[pivot * n + col]))
          pivot = row;
      if (!(std::abs(a[pivot * n + col]) > pivot_threshold))
        return false;
      if (pivot != col) {
        for (int k = 0; k < n; ++k)
          std::swap(a[col * n + k], a[pivot * n + k]);
        std::swap(b[col], b[pivot]);
      }
      for (int row = col + 1; row < n; ++row) {
        const double factor = a[row * n + col] / a[col * n + col];
        for (int k = col; k < n; ++k)
          a[row * n + k] -= factor * a[col * n + k];
        b[row] -= factor * b[col];
      }
    }

    x.assign(n, 0.);
    for (int row = n - 1; row >= 0; --row) {
      double sum = b[row];
      for (int k = row + 1; k < n; ++k)
        sum -= a[row * n + k] * x[k];
      x[row] = sum / a[row * n + row];
      if (!std::isfinite(x[row]))
        return false;
    }
    return true;
  }

  Config config_;
  int size_ = 0;
  std::deque<Entry> history_;
  std::vector<double> weights_;
  std::vector<std::complex<double>> x_new_;
  std::vector<double> last_coefficients_;
  bool last_fell_back_ = false;
  bool last_clipped_ = false;
  double last_residual_norm_ = 0.;
  double last_combined_residual_norm_ = 0.;
  double previous_residual_norm_ = 0.;
  int restarts_ = 0;
};

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SELF_ENERGY_MIXER_HPP
