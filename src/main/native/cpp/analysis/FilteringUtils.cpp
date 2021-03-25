// Copyright (c) FIRST and other WPILib contributors.
// Open Source Software; you can modify and/or share it under the terms of
// the WPILib BSD license file in the root directory of this project.

#include "sysid/analysis/FilteringUtils.h"

#include <frc/LinearFilter.h>

namespace sysid {
double GetAccelNoiseFloor(const std::vector<PreparedData>& data, int window) {
  double sum = 0.0;
  size_t step = window / 2;
  frc::LinearFilter<double> averageFilter =
      frc::LinearFilter<double>::MovingAverage(window);
  for (size_t i = 0; i < data.size(); i++) {
    double mean = averageFilter.Calculate(data[i].acceleration);
    if (i >= step) {
      sum += std::pow(data[i - step].acceleration - mean, 2);
    }
  }
  return std::sqrt(sum / (data.size() - step));
}

void TrimStepVoltageData(std::vector<PreparedData>* data,
                         AnalysisManager::Settings& settings,
                         double& minStepTime) {
  double firstTimestamp = data->at(0).timestamp;

  // Trim data before max acceleration
  data->erase(data->begin(),
              std::max_element(
                  data->begin(), data->end(), [](const auto& a, const auto& b) {
                    return std::abs(a.acceleration) < std::abs(b.acceleration);
                  }));

  minStepTime = std::min(data->at(0).timestamp - firstTimestamp, minStepTime);

  // If step duration hasn't been set yet, set calculate a default (find the
  // entry before the acceleration first hits zero)
  if (static_cast<double>(settings.stepTestDuration) <= minStepTime) {
    // Get noise floor
    const double accelNoiseFloor =
        GetAccelNoiseFloor(*data, settings.windowSize);
    // Find latest element with nonzero acceleration
    auto endIt = std::find_if(
        std::reverse_iterator{data->end()},
        std::reverse_iterator{data->begin()}, [&](const PreparedData& entry) {
          return std::abs(entry.acceleration) > accelNoiseFloor;
        });

    // Calculate default duration
    settings.stepTestDuration = static_cast<float>(
        endIt->timestamp - data->front().timestamp + minStepTime);
  }

  // Find first entry greater than the step test duration
  auto maxIt =
      std::find_if(data->begin(), data->end(), [&](PreparedData entry) {
        return entry.timestamp - data->front().timestamp + minStepTime >
               settings.stepTestDuration;
      });

  // Trim data beyond desired step test duration
  if (maxIt != data->end()) {
    data->erase(maxIt, data->end());
  }
}
}  // namespace sysid
