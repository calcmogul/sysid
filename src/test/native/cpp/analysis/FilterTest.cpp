// Copyright (c) FIRST and other WPILib contributors.
// Open Source Software; you can modify and/or share it under the terms of
// the WPILib BSD license file in the root directory of this project.

#include <array>
#include <vector>

#include <wpi/raw_ostream.h>

#include "gtest/gtest.h"
#include "sysid/analysis/AnalysisManager.h"
#include "sysid/analysis/FeedforwardAnalysis.h"
#include "sysid/analysis/FilteringUtils.h"

TEST(FilterTest, MedianFilter) {
  std::vector<std::array<double, 1>> test_data = {{0}, {1},    {10}, {5}, {3},
                                                  {0}, {1000}, {7},  {6}, {5}};
  std::vector<std::array<double, 1>> expected_data = {{1}, {5}, {5}, {3},
                                                      {3}, {7}, {7}, {6}};
  std::vector<std::array<double, 1>> filtered_data =
      sysid::ApplyMedianFilter<1, 0>(test_data, 3);
  EXPECT_EQ(expected_data, filtered_data);
}

TEST(FilterTest, QuasistaticTrim) {
  std::vector<std::array<double, 2>> test_data = {
      {0, 1}, {.5, 2}, {2, 0.1}, {4, 4}, {0, 5}};
  std::vector<std::array<double, 2>> expected_data = {{.5, 2}, {4, 4}};
  sysid::TrimQuasistaticData<2, 0, 1>(&test_data, 0.2);
  EXPECT_EQ(expected_data, test_data);
}

TEST(FilterTests, AccelNoiseFloor) {
  std::vector<sysid::PreparedData> test_data = {
      {0, 1, 2, 3, 0, 0}, {1, 1, 2, 3, 1, 0},    {2, 1, 2, 3, 2, 0},
      {3, 1, 2, 3, 5, 0}, {4, 1, 2, 3, 0.35, 0}, {5, 1, 2, 3, 0.15, 0},
      {6, 1, 2, 3, 0, 0}, {7, 1, 2, 3, 0.02, 0}, {8, 1, 2, 3, 0.01, 0},
      {9, 1, 2, 3, 0, 0}};
  double noiseFloor = GetAccelNoiseFloor(test_data, 2);
  EXPECT_NEAR(0.953, noiseFloor, 0.001);
}

TEST(FilterTest, StepTrim) {
  std::vector<sysid::PreparedData> test_data = {
      {0, 1, 2, 3, 0, 0},    {1, 1, 2, 3, 0.25, 0}, {2, 1, 2, 3, 0.5, 0},
      {3, 1, 2, 3, 0.45, 0}, {4, 1, 2, 3, 0.35, 0}, {5, 1, 2, 3, 0.15, 0},
      {6, 1, 2, 3, 0, 0},    {7, 1, 2, 3, 0.02, 0}, {8, 1, 2, 3, 0.01, 0},
      {9, 1, 2, 3, 0, 0},
  };

  std::vector<sysid::PreparedData> expected_data = {
      {2, 1, 2, 3, 0.5, 0},
      {3, 1, 2, 3, 0.45, 0},
      {4, 1, 2, 3, 0.35, 0},
  };

  double minTime = 100;

  sysid::AnalysisManager::Settings settings;
  sysid::TrimStepVoltageData(&test_data, settings, minTime);

  EXPECT_EQ(expected_data[0].acceleration, test_data[0].acceleration);
  EXPECT_EQ(expected_data.back().acceleration, test_data.back().acceleration);
  EXPECT_EQ(4, settings.stepTestDuration);
  EXPECT_EQ(2, minTime);
}
