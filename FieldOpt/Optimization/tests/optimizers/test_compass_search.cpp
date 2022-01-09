/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
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
#include <Runner/tests/test_resource_runner.hpp>
#include "Optimization/optimizers/compass_search.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;

namespace {

class CompassSearchTest : public ::testing::Test,
                          public TestResources::TestResourceOptimizer,
                          public TestResources::TestResourceGrids {
 protected:
  CompassSearchTest() {
    base_ = base_case_;
  }
  virtual ~CompassSearchTest() {}
  virtual void SetUp() {}

  Optimization::Case *base_;
};

TEST_F(CompassSearchTest, Constructor) {}

TEST_F(CompassSearchTest, GetNewCases) {
  test_case_1_3i_->set_objf_value(
    Sphere(test_case_1_3i_->GetRealVarVector()));

  Optimization::Optimizer *maximizer =
    new CompassSearch(settings_compass_search_max_unconstr_,
                      test_case_1_3i_,
                      varcont_prod_bhp_,
                      grid_5spot_,
                      logger_);

  // These four cases should change the values of the two first int vars,
  // +50 then -50
  Optimization::Case *new_case_1 = maximizer->GetCaseForEvaluation();
  Optimization::Case *new_case_2 = maximizer->GetCaseForEvaluation();
  Optimization::Case *new_case_3 = maximizer->GetCaseForEvaluation();
  Optimization::Case *new_case_4 = maximizer->GetCaseForEvaluation();

  EXPECT_FALSE(new_case_1->id() == new_case_2->id());
  EXPECT_FALSE(new_case_3->id() == new_case_4->id());

  EXPECT_EQ(test_case_1_3i_->GetIntegerVarVector()[0] + 8,
            new_case_1->GetIntegerVarVector()[0]);

  EXPECT_EQ(test_case_1_3i_->GetIntegerVarVector()[1] + 0,
            new_case_1->GetIntegerVarVector()[1]);

  EXPECT_EQ(test_case_1_3i_->GetIntegerVarVector()[2] + 0,
            new_case_1->GetIntegerVarVector()[2]);

  EXPECT_EQ(test_case_1_3i_->GetIntegerVarVector()[1] + 8,
            new_case_2->GetIntegerVarVector()[1]);

  EXPECT_EQ(test_case_1_3i_->GetIntegerVarVector()[2] + 8,
            new_case_3->GetIntegerVarVector()[2]);

  EXPECT_EQ(test_case_1_3i_->GetIntegerVarVector()[0] - 8,
            new_case_4->GetIntegerVarVector()[0]);
}

TEST_F(CompassSearchTest, TestFunctionSpherical) {
  test_case_2r_->set_objf_value(
    Sphere(test_case_2r_->GetRealVarVector()));

  Optimization::Optimizer *minimizer =
    new CompassSearch(settings_compass_search_min_unconstr_,
                      test_case_2r_,
                      varcont_prod_bhp_,
                      grid_5spot_,
                      logger_
    );

  while (!minimizer->IsFinished()) {
    auto next_case = minimizer->GetCaseForEvaluation();
    next_case->set_objf_value(Sphere(next_case->GetRealVarVector()));
    minimizer->SubmitEvaluatedCase(next_case);
  }
  auto best_case = minimizer->GetTentBestCase();
  EXPECT_NEAR(0.0, best_case->objf_value(), 0.01);
  EXPECT_NEAR(0.0, best_case->GetRealVarVector()[0], 0.01);
  EXPECT_NEAR(0.0, best_case->GetRealVarVector()[1], 0.01);
}

TEST_F(CompassSearchTest, TestFunctionRosenbrock) {

  // First test the Rosenbrock function itself
  Eigen::VectorXd optimum(2); optimum << 1.0, 1.0;
  EXPECT_FLOAT_EQ(0.0, Rosenbrock(optimum));

  test_case_2r_->set_objf_value(
    Rosenbrock(test_case_2r_->GetRealVarVector()));

  Optimization::Optimizer *minimizer =
    new CompassSearch(settings_compass_search_min_unconstr_,
                      test_case_2r_,
                      varcont_prod_bhp_,
                      grid_5spot_,
                      logger_
    );

  while (!minimizer->IsFinished()) {
    auto next_case = minimizer->GetCaseForEvaluation();
    next_case->set_objf_value(Rosenbrock(next_case->GetRealVarVector()));
    minimizer->SubmitEvaluatedCase(next_case);
  }
  auto best_case = minimizer->GetTentBestCase();

  // The Rosenbrock function is hard. We don't expect
  // Compass search to find the optimum exactly.
  EXPECT_NEAR(0.0, best_case->objf_value(), 5);
  EXPECT_NEAR(1.0, best_case->GetRealVarVector()[0], 2.5);
  EXPECT_NEAR(1.0, best_case->GetRealVarVector()[1], 2.5);
}

}

