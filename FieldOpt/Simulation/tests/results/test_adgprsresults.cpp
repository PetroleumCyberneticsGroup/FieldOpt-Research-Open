/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb-pcg@gmail.com>

This file is part of the FieldOpt project.

FieldOpt is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version
3 of the License, or (at your option) any later version.

FieldOpt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the GNU General Public License for more details.

You should have received a copy of the
GNU General Public License along with FieldOpt.
If not, see <http://www.gnu.org/licenses/>.
***********************************************************/

#include <gtest/gtest.h>
#include "Settings/tests/test_resource_example_file_paths.hpp"
#include "Simulation/results/adgprsresults.h"
#include "Settings/tests/test_resource_settings.hpp"

namespace {

class AdgprsResultsTest : public ::testing::Test,
                          public TestResources::TestResourceSettings {
 protected:
  AdgprsResultsTest() {
    results_ = new Simulation::Results::AdgprsResults(settings_simulator_);
    auto fp = QString::fromStdString(TestResources::ExampleFilePaths::gprs_base_5spot_);
    results_->ReadResults(fp);
  }
  virtual ~AdgprsResultsTest() {}
  virtual void SetUp() {}
  Simulation::Results::AdgprsResults *results_;
};

TEST_F(AdgprsResultsTest, ReadFile) {
  EXPECT_TRUE(true);
}

TEST_F(AdgprsResultsTest, Time) {
  EXPECT_FLOAT_EQ(0.0, results_->GetValue(Simulation::Results::Results::Property::Time, 0));
  EXPECT_FLOAT_EQ(100.0, results_->GetValue(Simulation::Results::Results::Property::Time));
  EXPECT_EQ(8, results_->GetValueVector(Simulation::Results::Results::Property::Time).size());
}

TEST_F(AdgprsResultsTest, FOPT) {
  EXPECT_EQ(8, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeOilProduction).size());
  EXPECT_FLOAT_EQ(0.0, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeOilProduction).front());
  EXPECT_FLOAT_EQ(178150.12, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeOilProduction).back());
}

TEST_F(AdgprsResultsTest, FWPT) {
  EXPECT_EQ(8, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeWaterProduction).size());
  EXPECT_FLOAT_EQ(0.0, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeWaterProduction).front());
  EXPECT_FLOAT_EQ(305866.53, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeWaterProduction).back());
}

TEST_F(AdgprsResultsTest, FGPT) {
  EXPECT_EQ(8, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeGasProduction).size());
  EXPECT_FLOAT_EQ(0.0, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeGasProduction).front());
  EXPECT_FLOAT_EQ(0.0, results_->GetValueVector(Simulation::Results::Results::Property::CumulativeGasProduction).back());
}

}
