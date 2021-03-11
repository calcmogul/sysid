// Copyright (c) FIRST and other WPILib contributors.
// Open Source Software; you can modify and/or share it under the terms of
// the WPILib BSD license file in the root directory of this project.

#include "sysid/analysis/FeedforwardAnalysis.h"

#include <cmath>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <fmt/core.h>
#include <units/math.h>
#include <units/time.h>

#include "sysid/analysis/AnalysisManager.h"
#include "sysid/analysis/FilteringUtils.h"
#include "sysid/analysis/OLS.h"

using namespace sysid;

/**
 * Populates OLS vector for u = Ks + Kv v + Ka a.
 *
 * @param d       List of characterization data.
 * @param type    Type of system being identified.
 * @param olsData Vector of OLS data.
 */
static void PopulateAccelOLSVector(const std::vector<PreparedData>& d,
                                   const AnalysisType& type,
                                   std::vector<double>& olsData) {
  for (const auto& pt : d) {
    // u = Ks sgn(v) + K_v v + Ka a

    // Add the dependent variable (voltage)
    olsData.push_back(pt.voltage);

    // Add the intercept term (for Ks)
    olsData.push_back(std::copysign(1, pt.velocity));

    // Add the velocity term (for Kv)
    olsData.push_back(pt.velocity);

    // Add the acceleration term (for Ka)
    olsData.push_back(pt.acceleration);

    if (type == analysis::kArm) {
      // Add the cosine term (for Kcos)
      olsData.push_back(pt.cos);
    }
  }
}

/**
 * Populates OLS vector for x_k+1 = alpha x_k + beta u_k + gamma sgn(x_k).
 *
 * @param d       List of characterization data.
 * @param type    Type of system being identified.
 * @param olsData Vector of OLS data.
 */
static void PopulateNextVelOLSVector(const std::vector<PreparedData>& d,
                                     const AnalysisType& type,
                                     std::vector<double>& olsData) {
  for (const auto& pt : d) {
    // x_k+1 = alpha x_k + beta u_k + gamma sgn(x_k)

    // Add the dependent variable (next velocity)
    olsData.push_back(pt.nextVelocity);

    // Add the velocity term (for alpha)
    olsData.push_back(pt.velocity);

    // Add the voltage term (for beta)
    olsData.push_back(pt.voltage);

    // Add the intercept term (for gamma)
    olsData.push_back(std::copysign(1, pt.velocity));

    // Add test-specific variables
    if (type == analysis::kElevator) {
      // Add the gravity term (for Kg)
      olsData.push_back(1.0);
    }
  }
}

static std::tuple<std::vector<double>, double> RefineGains(
    const std::vector<double>& data, units::second_t dtMean,
    size_t independentVariables, const std::vector<double>& initialGains) {
  fmt::print(stderr, "RefineGains() initial guess\n");
  fmt::print(stderr, "  {}, {}, {}\n", initialGains[0], initialGains[1],
             initialGains[2]);

  // Perform some quick sanity checks regarding the size of the vector.
  assert(data.size() % (independentVariables + 1) == 0);

  // Get the number of elements.
  size_t n = data.size() / (independentVariables + 1);

  // Create new variables to make things more readable.
  size_t rows = n;
  size_t cols = independentVariables;  // X
  size_t strd = independentVariables + 1;

  // Create y and X matrices.
  Eigen::Map<const Eigen::MatrixXd, 0, Eigen::Stride<1, Eigen::Dynamic>> y(
      data.data() + 0, rows, 1, Eigen::Stride<1, Eigen::Dynamic>(1, strd));

  Eigen::Map<const Eigen::MatrixXd, 0, Eigen::Stride<1, Eigen::Dynamic>> X(
      data.data() + 1, rows, cols, Eigen::Stride<1, Eigen::Dynamic>(1, strd));

  // Implement the Gauss-Newton algorithm for nonlinear least squares
  Eigen::Matrix<double, 3, 1> beta;
  // beta << initialGains[0], initialGains[1] - 0.03e-2, initialGains[2]
  // + 9.9e-3;
  beta << initialGains[0], initialGains[1], initialGains[2];
  Eigen::Matrix<double, 3, 1> nextBeta;

  auto f = [=](double x, double u, double T,
               const Eigen::Matrix<double, 3, 1>& beta) {
    const double& Ks = beta(0);
    const double& Kv = beta(1);
    const double& Ka = beta(2);

    double alpha = std::exp(-Kv / Ka * T);
    return alpha * x + 1.0 / Kv * (1.0 - alpha) * u +
           Ks / Kv * (alpha - 1.0) * wpi::sgn(x);
  };

  double T = dtMean.to<double>();

  Eigen::MatrixXd J{rows, independentVariables};
  bool repeat = true;
  while (repeat) {
    for (int i = 0; i < rows; ++i) {
      const double& x = X(i, 0);
      const double& u = X(i, 1);

      const double& Ks = beta(0);
      const double& Kv = beta(1);
      const double& Ka = beta(2);

      double a = std::exp(-Kv / Ka * T);
      double b = Ks * wpi::sgn(x) - u + Kv * x;
      J(i, 0) = 1.0 / Kv * (a - 1.0) * wpi::sgn(x);
      J(i, 1) = a *
                (Ka * (Ks * wpi::sgn(x) - u) * (std::exp(Kv / Ka * T) - 1.0) -
                 T * Kv * b) /
                (Ka * Kv * Kv);
      J(i, 2) = T * a * b / (Ka * Ka);
    }

    Eigen::MatrixXd b = (J.transpose() * J).llt().solve(J.transpose());

    Eigen::MatrixXd r{rows, 1};
    for (int i = 0; i < rows; ++i) {
      r(i) = y(i, 0) - f(X(i, 0), X(i, 1), T, beta);
    }

    nextBeta = beta + b * r;
    if ((nextBeta - beta).norm() < 1e-5) {
      repeat = false;
    }
    beta = nextBeta;

    fmt::print(stderr, "Step\n");
    fmt::print(stderr, "  {}, {}, {}\n", beta(0), beta(1), beta(2));
  }

  Eigen::MatrixXd r{rows, 1};
  for (int i = 0; i < rows; ++i) {
    r(i) = y(i, 0) - f(X(i, 0), X(i, 1), T, beta);
  }

  // We will now calculate r^2 or the coefficient of determination, which
  // tells us how much of the total variation (variation in y) can be
  // explained by the regression model.

  // We will first calculate the sum of the squares of the error, or the
  // variation in error (SSE).
  double SSE = r.squaredNorm();

  // Now we will calculate the total variation in y, known as SSTO.
  double SSTO = ((y.transpose() * y) - (1 / n) * (y.transpose() * y)).value();

  double rSquared = (SSTO - SSE) / SSTO;
  double adjRSquared = 1 - (1 - rSquared) * ((n - 1.0) / (n - 3));

  std::vector<double> b;
  for (int i = 0; i < 3; ++i) {
    b.emplace_back(beta(i));
  }
  return {b, adjRSquared};
}

std::tuple<std::vector<double>, double> sysid::CalculateFeedforwardGains(
    const Storage& data, const AnalysisType& type) {
  units::second_t dtMean = GetMeanTimeDelta(data);

  // Create a raw vector of doubles with our data in it.
  std::vector<double> olsData;
  olsData.reserve((1 + type.independentVariables) *
                  (data.slow.size() + data.fast.size()));

  // Iterate through the data and add it to our raw vector.
  const auto& [slow, fast] = data;
  if (type == analysis::kArm) {
    PopulateAccelOLSVector(slow, type, olsData);
    PopulateAccelOLSVector(fast, type, olsData);

    // Gains are Ks, Kv, Ka, Kcos
    return sysid::OLS(olsData, type.independentVariables);
  } else {
    // This implements the OLS algorithm defined in
    // https://file.tavsys.net/control/sysid-ols.pdf.

    PopulateNextVelOLSVector(slow, type, olsData);
    PopulateNextVelOLSVector(fast, type, olsData);

    auto ols = sysid::OLS(olsData, type.independentVariables);
    double alpha = std::get<0>(ols)[0];
    double beta = std::get<0>(ols)[1];
    double gamma = std::get<0>(ols)[2];

    // Initialize gains list with Ks, Kv, and Ka
    std::vector<double> gains{
        {-gamma / beta, (1 - alpha) / beta,
         dtMean.to<double>() * (alpha - 1) / (beta * std::log(alpha))}};

    if (type == analysis::kElevator) {
      // Add Kg to gains list
      double delta = std::get<0>(ols)[3];
      gains.emplace_back(-delta / beta);
    }

    return RefineGains(olsData, dtMean, type.independentVariables, gains);
  }
}
