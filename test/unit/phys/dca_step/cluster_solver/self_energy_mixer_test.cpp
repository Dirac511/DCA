// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// This file tests self_energy_mixer.hpp.

#include "dca/phys/dca_step/cluster_solver/self_energy_mixer.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "gtest/gtest.h"

namespace {

using dca::phys::solver::SelfEnergyMixer;
using Complex = std::complex<double>;

// Applies the mixer to the linear contraction x -> fixed + 0.8 * (x - fixed).
void mixLinearContraction(SelfEnergyMixer& mixer, std::vector<Complex>& x, int iterations) {
  const Complex fixed(1., -0.5);
  const int n = static_cast<int>(x.size());
  for (int it = 0; it < iterations; ++it) {
    std::vector<Complex> y(x.size());
    for (int i = 0; i < n; ++i)
      y[i] = fixed + 0.8 * (x[i] - fixed);
    mixer.mix(y.data(), x.data(), x.data(), n);
  }
}

double distanceTo(const std::vector<Complex>& x, const Complex& fixed) {
  double sum = 0.;
  for (std::size_t i = 0; i < x.size(); ++i)
    sum += std::norm(x[i] - fixed);
  return std::sqrt(sum);
}

}  // namespace

TEST(SelfEnergyMixerTest, ThrowsOnInvalidConfiguration) {
  SelfEnergyMixer mixer;
  EXPECT_THROW(mixer.configure(1., 0), std::logic_error);

  std::vector<Complex> v{Complex(0., 0.)};
  EXPECT_THROW(mixer.mix(v.data(), v.data(), v.data(), 1), std::logic_error);

  SelfEnergyMixer::Config config;
  config.alpha = -1.;
  config.max_history = 4;
  EXPECT_THROW(mixer.configure(config), std::logic_error);

  config.alpha = 1.;
  config.max_coefficient = 0.5;
  EXPECT_THROW(mixer.configure(config), std::logic_error);

  config.max_coefficient = 10.;
  config.regularization = -1.;
  EXPECT_THROW(mixer.configure(config), std::logic_error);
}

TEST(SelfEnergyMixerTest, SingleHistoryReducesToLinearMixing) {
  SelfEnergyMixer mixer;
  mixer.configure(0.5, 1);

  const std::vector<Complex> qmc{Complex(2., 1.), Complex(-1., 0.5), Complex(0.25, -3.)};
  const std::vector<Complex> in{Complex(0.5, 0.), Complex(0.25, 0.1), Complex(-1., 2.)};
  std::vector<Complex> out(in.size());

  mixer.mix(qmc.data(), in.data(), out.data(), static_cast<int>(out.size()));

  for (std::size_t i = 0; i < out.size(); ++i) {
    const Complex expected = in[i] + 0.5 * (qmc[i] - in[i]);
    EXPECT_NEAR(0., std::abs(out[i] - expected), 1.e-14);
  }

  ASSERT_EQ(1u, mixer.lastCoefficients().size());
  EXPECT_EQ(1., mixer.lastCoefficients()[0]);
  EXPECT_TRUE(mixer.lastFellBack());
}

TEST(SelfEnergyMixerTest, IsInPlaceSafe) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 1);

  std::vector<Complex> qmc{Complex(2., 1.), Complex(-1., 0.5)};
  const std::vector<Complex> in{Complex(0.5, 0.), Complex(0.25, 0.)};
  const std::vector<Complex> qmc_ref = qmc;

  // With alpha == 1 and a single history entry the output equals the QMC self-energy.
  mixer.mix(qmc.data(), in.data(), qmc.data(), static_cast<int>(qmc.size()));

  for (std::size_t i = 0; i < qmc.size(); ++i)
    EXPECT_NEAR(0., std::abs(qmc[i] - qmc_ref[i]), 1.e-14);
}

// On a linear contraction x -> x* + 0.5 (x - x*) the Anderson/Pulay mixer must converge to the
// fixed point in a couple of iterations (Aitken-like extrapolation), much faster than the plain
// linear mixing (which only halves the error each step).
TEST(SelfEnergyMixerTest, AndersonAcceleratesLinearContraction) {
  const int n = 4;
  const std::vector<Complex> fixed(n, Complex(1., -0.5));

  const auto compute_new = [&](const std::vector<Complex>& x) {
    std::vector<Complex> y(n);
    for (int i = 0; i < n; ++i)
      y[i] = fixed[i] + 0.5 * (x[i] - fixed[i]);
    return y;
  };

  SelfEnergyMixer mixer;
  mixer.configure(1., 4);

  std::vector<Complex> x(n, Complex(0., 0.));
  double residual = 0.;
  for (int it = 0; it < 4; ++it) {
    const std::vector<Complex> y = compute_new(x);
    mixer.mix(y.data(), x.data(), x.data(), n);

    residual = 0.;
    for (int i = 0; i < n; ++i)
      residual += std::norm(x[i] - fixed[i]);
    residual = std::sqrt(residual);
  }

  EXPECT_LT(residual, 1.e-6);
}

// Two history entries are enough to extrapolate a linear problem exactly (up to the Tikhonov
// regularization).
TEST(SelfEnergyMixerTest, TwoHistoriesAreExactForALinearProblem) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 2);

  std::vector<Complex> x(3, Complex(0., 0.));
  mixLinearContraction(mixer, x, 2);

  EXPECT_LT(distanceTo(x, Complex(1., -0.5)), 1.e-6);
}

// The DIIS coefficients are constrained, and clipping plus renormalization keeps that property.
TEST(SelfEnergyMixerTest, CoefficientsSumToOne) {
  SelfEnergyMixer mixer;
  mixer.configure(0.7, 5);

  const int n = 5;
  std::vector<Complex> x(n, Complex(0., 0.));
  for (int it = 0; it < 6; ++it) {
    std::vector<Complex> y(n);
    for (int i = 0; i < n; ++i)
      y[i] = Complex(1.2, -0.3) + 0.6 * (x[i] - Complex(1.2, -0.3)) +
             0.01 * Complex(double((it + 1) * (i + 1) % 3 - 1), 0.);
    mixer.mix(y.data(), x.data(), x.data(), n);

    double sum = 0.;
    for (std::size_t j = 0; j < mixer.lastCoefficients().size(); ++j)
      sum += mixer.lastCoefficients()[j];
    EXPECT_NEAR(1., sum, 1.e-12);
  }
}

// Identical residuals make the DIIS system exactly singular. Without regularization the mixer must
// detect that and fall back to the linear step instead of producing a huge extrapolation.
TEST(SelfEnergyMixerTest, CollinearResidualsFallBackToLinearMixing) {
  SelfEnergyMixer::Config config;
  config.alpha = 1.;
  config.max_history = 4;
  config.regularization = 0.;
  SelfEnergyMixer mixer;
  mixer.configure(config);

  const int n = 3;
  std::vector<Complex> x(n, Complex(0., 0.));
  for (int it = 0; it < 3; ++it) {
    std::vector<Complex> y(n);
    for (int i = 0; i < n; ++i)
      y[i] = x[i] + Complex(1., 0.5);
    mixer.mix(y.data(), x.data(), x.data(), n);
  }

  EXPECT_TRUE(mixer.lastFellBack());
  for (int i = 0; i < n; ++i) {
    EXPECT_TRUE(std::isfinite(x[i].real()));
    EXPECT_LT(std::abs(x[i]), 10.);
  }
}

// With the default regularization the same singular configuration must stay bounded.
TEST(SelfEnergyMixerTest, CollinearResidualsStayBounded) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 4);

  const int n = 3;
  std::vector<Complex> x(n, Complex(0., 0.));
  for (int it = 0; it < 3; ++it) {
    std::vector<Complex> y(n);
    for (int i = 0; i < n; ++i)
      y[i] = x[i] + Complex(1., 0.5);
    mixer.mix(y.data(), x.data(), x.data(), n);
  }

  for (std::size_t j = 0; j < mixer.lastCoefficients().size(); ++j)
    EXPECT_LE(std::abs(mixer.lastCoefficients()[j]), 10.);
  for (int i = 0; i < n; ++i)
    EXPECT_TRUE(std::isfinite(std::abs(x[i])));
}

TEST(SelfEnergyMixerTest, ThrowsOnNonFiniteSelfEnergy) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 4);

  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<Complex> qmc{Complex(nan, 0.), Complex(1., 0.)};
  std::vector<Complex> in{Complex(0., 0.), Complex(1., 0.)};
  std::vector<Complex> out(2);

  EXPECT_THROW(mixer.mix(qmc.data(), in.data(), out.data(), 2), std::runtime_error);
  EXPECT_EQ(0u, mixer.historySize());
}

TEST(SelfEnergyMixerTest, RestartsHistoryWhenResidualGrows) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 8);

  const int n = 2;
  std::vector<Complex> x(n, Complex(0., 0.));
  std::vector<Complex> y{Complex(1., 0.), Complex(1., 0.)};

  mixer.mix(y.data(), x.data(), x.data(), n);
  EXPECT_EQ(0, mixer.restarts());

  // A residual that is ten times larger than the previous one must drop the history.
  for (int i = 0; i < n; ++i)
    y[i] = x[i] + Complex(10., 0.);
  mixer.mix(y.data(), x.data(), x.data(), n);

  EXPECT_EQ(1, mixer.restarts());
  EXPECT_EQ(1u, mixer.historySize());
}

TEST(SelfEnergyMixerTest, ThrowsOnSizeChange) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 4);

  std::vector<Complex> a(2, Complex(0., 0.));
  mixer.mix(a.data(), a.data(), a.data(), 2);

  std::vector<Complex> b(3, Complex(0., 0.));
  EXPECT_THROW(mixer.mix(b.data(), b.data(), b.data(), 3), std::logic_error);
}

TEST(SelfEnergyMixerTest, ThrowsOnInvalidResidualWeights) {
  SelfEnergyMixer mixer;
  EXPECT_THROW(mixer.setResidualWeights(std::vector<double>{1., -1.}), std::logic_error);

  mixer.configure(1., 4);
  mixer.setResidualWeights(std::vector<double>{1., 0., 1.});
  std::vector<Complex> a(2, Complex(0., 0.));
  EXPECT_THROW(mixer.mix(a.data(), a.data(), a.data(), 2), std::logic_error);
}

// The weights select the elements that determine the coefficients.
TEST(SelfEnergyMixerTest, ResidualWeightsSelectElements) {
  const int n = 2;
  SelfEnergyMixer uniform;
  uniform.configure(1., 2);
  SelfEnergyMixer weighted;
  weighted.configure(1., 2);
  weighted.setResidualWeights(std::vector<double>{1., 0.});

  const std::vector<Complex> x0(n, Complex(0., 0.));
  const std::vector<Complex> y0{Complex(1., 0.), Complex(50., 0.)};
  std::vector<Complex> xu(x0);
  std::vector<Complex> xw(x0);
  uniform.mix(y0.data(), x0.data(), xu.data(), n);
  weighted.mix(y0.data(), x0.data(), xw.data(), n);

  // Residuals: r0 = 1, r1 = 2 in the first element.
  const std::vector<Complex> y1{Complex(3., 0.), Complex(120., 0.)};
  uniform.mix(y1.data(), xu.data(), xu.data(), n);
  weighted.mix(y1.data(), xw.data(), xw.data(), n);

  // Analytic one-dimensional solution of the constrained minimization.
  const double r0 = 1.;
  const double r1 = 2.;
  const double delta = r0 - r1;
  const double expected = -r1 * delta / (delta * delta);

  ASSERT_EQ(2u, weighted.lastCoefficients().size());
  EXPECT_NEAR(expected, weighted.lastCoefficients()[0], 1.e-6);
  EXPECT_GT(std::abs(uniform.lastCoefficients()[0] - weighted.lastCoefficients()[0]), 1.e-3);
}

TEST(SelfEnergyMixerTest, ResetClearsTheState) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 4);

  std::vector<Complex> v{Complex(1., 0.), Complex(2., 0.)};
  mixer.mix(v.data(), v.data(), v.data(), 2);
  EXPECT_EQ(1u, mixer.historySize());

  mixer.reset();

  EXPECT_EQ(0u, mixer.historySize());
  EXPECT_TRUE(mixer.lastCoefficients().empty());
  EXPECT_EQ(0, mixer.restarts());
  EXPECT_FALSE(mixer.lastFellBack());
}

// The residual norm reported to the caller must be the un-mixed residual r = F(x) - x.
TEST(SelfEnergyMixerTest, ReportsTheUnmixedResidualNorm) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 4);

  const int n = 2;
  std::vector<Complex> x(n, Complex(0., 0.));
  std::vector<Complex> y{Complex(3., 4.), Complex(0., 0.)};

  mixer.mix(y.data(), x.data(), x.data(), n);

  EXPECT_NEAR(5., mixer.lastResidualNorm(), 1.e-12);
  EXPECT_TRUE(std::isfinite(mixer.lastCombinedResidualNorm()));
}
