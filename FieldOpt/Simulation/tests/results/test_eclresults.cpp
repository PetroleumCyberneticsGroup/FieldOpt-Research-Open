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
#include <QString>
#include <qvector.h>
#include "Simulation/results/eclresults.h"
#include "Settings/tests/test_resource_example_file_paths.hpp"
#include "Settings/tests/test_resource_settings.hpp"

using namespace Simulation::Results;

namespace {

class ECLResultsTest : public ::testing::Test,
                       public TestResources::TestResourceSettings {
 protected:

  ECLResultsTest() {
    results_ = new ECLResults(settings_simulator_);
  }

  virtual ~ECLResultsTest() { }

  virtual void SetUp() { }

  virtual void TearDown() {
    results_->DumpResults();
  }

  Results *results_;
};

TEST_F(ECLResultsTest, ReadSummary) {
  EXPECT_THROW(results_->ReadResults("a"),
               ResultFileNotFoundException);
  auto fp = QString::fromStdString(TestResources::ExampleFilePaths::ecl_base_horzwell);
  EXPECT_NO_THROW(results_->ReadResults(fp));
}

TEST_F(ECLResultsTest, DumpingAndAvailability) {
  // Make results available
  auto fp = QString::fromStdString(TestResources::ExampleFilePaths::ecl_base_horzwell);
  results_->ReadResults(fp);
  EXPECT_TRUE(results_->isAvailable());
  EXPECT_NO_THROW(results_->GetValue(ECLResults::Property::Time));

  // Make results unavailable
  results_->DumpResults();
  EXPECT_FALSE(results_->isAvailable());
  EXPECT_THROW(results_->GetValue(ECLResults::Property::Time),
               RsltsNotAvailExc);

  // Make results available again
  results_->ReadResults(fp);
  EXPECT_TRUE(results_->isAvailable());
  EXPECT_NO_THROW(results_->GetValue(ECLResults::Property::Time));
}

TEST_F(ECLResultsTest, FieldVariables) {
  auto fp = QString::fromStdString(TestResources::ExampleFilePaths::ecl_base_horzwell);
  results_->ReadResults(fp);
  EXPECT_FLOAT_EQ(187866.44, results_->GetValue(ECLResults::Property::FieldOilProdTotal));
  EXPECT_FLOAT_EQ(115870.73, results_->GetValue(ECLResults::Property::FieldOilProdTotal, 10));
  EXPECT_THROW(results_->GetValue(ECLResults::Property::FieldOilProdTotal, 30), std::runtime_error);

  std::vector<double> fopt_vec = results_->GetValueVector(ECLResults::Property::FieldOilProdTotal);
  EXPECT_EQ(21, fopt_vec.size());
  EXPECT_FLOAT_EQ(0, fopt_vec.front());
  EXPECT_FLOAT_EQ(187866.44, fopt_vec.back());
}

TEST_F(ECLResultsTest, MiscVariables) {
  auto fp = QString::fromStdString(TestResources::ExampleFilePaths::ecl_base_horzwell);
  results_->ReadResults(fp);
  EXPECT_FLOAT_EQ(200, results_->GetValue(ECLResults::Property::Time));
  EXPECT_FLOAT_EQ(100, results_->GetValue(ECLResults::Property::Time, 10));
  EXPECT_THROW(results_->GetValue(ECLResults::Property::Time, 30), std::runtime_error);

  std::vector<double> time_vec = results_->GetValueVector(ECLResults::Property::Time);
  EXPECT_EQ(21, time_vec.size());
  EXPECT_EQ(0, time_vec.front());
  EXPECT_EQ(200, time_vec.back());
}

TEST_F(ECLResultsTest, WellVariables) {
  auto fp = QString::fromStdString(TestResources::ExampleFilePaths::ecl_base_horzwell);
  results_->ReadResults(fp);
  EXPECT_FLOAT_EQ(1116.8876, results_->GetValue(ECLResults::Property::WellWatProdTotal, "PROD"));
  EXPECT_FLOAT_EQ(524.5061, results_->GetValue(ECLResults::Property::WellWatProdTotal, "PROD", 10));
}

}
