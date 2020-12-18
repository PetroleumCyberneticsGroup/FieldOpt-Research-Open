/***********************************************************
Copyright (C) 2017
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
#include "constraints/interw_dist.h"
#include "Optimization/tests/test_resource_cases.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Settings/tests/test_resource_settings.hpp"

namespace {

class InterwellDistanceTest : public ::testing::Test,
                              public TestResources::TestResourceCases
{

 public:
  InterwellDistanceTest() {
    iwd_settings_.type = Settings::Optimizer::ConstraintType::WSplineInterwDist;
    iwd_settings_.min = 100;
    iwd_settings_.wells << "TESTW" << "INJE";
    iwd_constraint_ = new Optimization::Constraints::InterwDist(iwd_settings_,
                                                                varcont_two_spline_wells_,
                                                                vp_);
  }

  Settings::Optimizer::Constraint iwd_settings_;
  Optimization::Constraints::InterwDist *iwd_constraint_;
  Optimization::Case *test_case = test_case_two_well_splines_;
  Settings::VerbParams vp_ = {};

  virtual ~InterwellDistanceTest() {}
  virtual void TearDown() {}
  virtual void SetUp() {}
};

TEST_F(InterwellDistanceTest, Penalty) {
  double penalty = iwd_constraint_->Penalty(test_case);
  EXPECT_EQ(160, penalty);

  iwd_settings_.min = 10;
  iwd_constraint_ = new Optimization::Constraints::InterwDist(iwd_settings_,
                                                              varcont_two_spline_wells_,
                                                              vp_);
  penalty = iwd_constraint_->Penalty(test_case);
  EXPECT_EQ(0, penalty);
}

TEST_F(InterwellDistanceTest, Normalization) {
  double penalty = iwd_constraint_->Penalty(test_case);
  iwd_constraint_->InitializeNormalizer(QList<Optimization::Case *>({test_case}));
  long double penalty_norm = iwd_constraint_->PenaltyNormalized(test_case);
  EXPECT_NEAR(1.0, penalty_norm, 0.005);
}

}
