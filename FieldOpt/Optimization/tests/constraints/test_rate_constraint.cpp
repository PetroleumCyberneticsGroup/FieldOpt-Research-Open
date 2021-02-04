/***********************************************************
Copyright (C) 2015-2016
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2021 Mathias Bellout
<chakibbb.pcg@gmail.com>

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
#include <constraints/rate_constraint.h>
#include "test_resource_optimizer.h"

namespace {
class RateConstraintTest : public ::testing::Test,
                           public TestResources::TestResourceOptimizer {

 public:
  RateConstraintTest() : Test() {
    auto *rate_var = new Model::Properties::ContinuousProperty(1200.0);
    rate_var->setName("Rate#INJ#0");
    model_->variables()->AddVariable(rate_var);
    constraint_ = new Optimization::Constraints::RateConstraint(settings_optimizer_->constraints()[2],
                                                                model_->variables(),
                                                                vp);
    c = new Optimization::Case(model_->variables()->GetBinVarValues(),
                               model_->variables()->GetDiscVarValues(),
                               model_->variables()->GetContVarValues());
  }
  Optimization::Constraints::RateConstraint *constraint_;
  Optimization::Case *c;
  Settings::VerbParams vp = settings_full_->global()->verbParams();
};


TEST_F(RateConstraintTest, Constructor) {
  EXPECT_TRUE(true);
}

TEST_F(RateConstraintTest, Initial) {
  EXPECT_TRUE(constraint_->CaseSatisfiesConstraint(c));
}

TEST_F(RateConstraintTest, AfterModification) {
  EXPECT_TRUE(constraint_->CaseSatisfiesConstraint(c));
  auto rate_vars = model_->variables()->GetWellRateVariables("INJ");
  for (auto var : rate_vars) {
    c->set_real_variable_value(var->id(), 900);
  }
  EXPECT_FALSE(constraint_->CaseSatisfiesConstraint(c));
  for (auto var : rate_vars) {
    c->set_real_variable_value(var->id(), 1300);
  }
  EXPECT_TRUE(constraint_->CaseSatisfiesConstraint(c));
  for (auto var : rate_vars) {
    c->set_real_variable_value(var->id(), 1600);
  }
  EXPECT_FALSE(constraint_->CaseSatisfiesConstraint(c));
}


TEST_F(RateConstraintTest, Snapping) {
  EXPECT_TRUE(constraint_->CaseSatisfiesConstraint(c));
  auto rate_vars = model_->variables()->GetWellRateVariables("INJ");

  // Testing snap to minimum
  for (auto var : rate_vars) {
    c->set_real_variable_value(var->id(), 900);
  }

  EXPECT_FALSE(constraint_->CaseSatisfiesConstraint(c));
  constraint_->SnapCaseToConstraints(c);
  EXPECT_TRUE(constraint_->CaseSatisfiesConstraint(c));

  for (auto var : rate_vars) {
    EXPECT_FLOAT_EQ(1200, c->get_real_variable_value(var->id()));
  }

  // Testing snap to maximum
  for (auto var : rate_vars) {
    c->set_real_variable_value(var->id(), 1800);
  }

  EXPECT_FALSE(constraint_->CaseSatisfiesConstraint(c));
  constraint_->SnapCaseToConstraints(c);
  EXPECT_TRUE(constraint_->CaseSatisfiesConstraint(c));

  for (auto var : rate_vars) {
    EXPECT_FLOAT_EQ(1400, c->get_real_variable_value(var->id()));
  }
}
}
