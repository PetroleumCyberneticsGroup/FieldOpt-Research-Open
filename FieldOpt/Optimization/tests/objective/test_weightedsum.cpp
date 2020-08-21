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
#include "Optimization/objective/weightedsum.h"
#include "Simulation/tests/test_resource_results.h"
#include "Settings/tests/test_resource_settings.hpp"
#include "Model/tests/test_resource_model.h"


using namespace Optimization::Objective;
using namespace Simulation::Results;
using namespace Settings;

namespace {

class WeightedSumTest :
  public ::testing::Test, public TestResources::TestResourceResults,
  public TestResources::TestResourceModel {

 protected:
  WeightedSumTest() { }

  virtual ~WeightedSumTest() { }

  virtual void SetUp() { }

  virtual void TearDown() { }

};

TEST_F(WeightedSumTest, Value) {
  auto *obj = new WeightedSum(settings_optimizer_, results_ecl_horzwell_, model_);
  float wwpt = results_ecl_horzwell_->GetValue(Results::Property::CumulativeWellWaterProduction, "PROD", 10);
  float fopt = results_ecl_horzwell_->GetValue(Results::Property::CumulativeOilProduction);
  EXPECT_FLOAT_EQ(100.0, results_ecl_horzwell_->GetValue(Results::Property::Time, 10));
  EXPECT_FLOAT_EQ(wwpt, 524.5061);
  EXPECT_FLOAT_EQ(fopt, 187866.44);
  EXPECT_FLOAT_EQ(fopt - 0.2*wwpt, obj->value());
}

}
