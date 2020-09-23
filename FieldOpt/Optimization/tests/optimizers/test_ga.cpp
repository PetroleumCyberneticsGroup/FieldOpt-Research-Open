/***********************************************************
Copyright (C) 2019
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
#include <Runner/tests/test_resource_runner.hpp>
#include "optimizers/RGARDD.h"
#include "Optimization/optimizers/GeneticAlgorithm.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;

namespace {

class GeneticAlgorithmTest : public ::testing::Test,
                             public TestResources::TestResourceOptimizer,
                             public TestResources::TestResourceGrids
{
 protected:
  GeneticAlgorithmTest() {
    base_ = base_case_;
  }
  virtual ~GeneticAlgorithmTest() {}
  virtual void SetUp() {}

  Optimization::Case *base_;

};

TEST_F(GeneticAlgorithmTest, Constructor) {
}

TEST_F(GeneticAlgorithmTest, TestFunctionSpherical) {
//    test_case_2r_->set_objective_function_value(abs(Sphere(test_case_2r_->GetRealVarVector())));
//    test_case_ga_spherical_30r_->set_objective_function_value(abs(Sphere(test_case_ga_spherical_30r_->GetRealVarVector())));
  test_case_ga_spherical_6r_->set_objf_value(abs(Sphere(test_case_ga_spherical_6r_->GetRealVarVector())));
  Optimization::Optimizer *minimizer = new RGARDD(settings_ga_min_,
                                                  test_case_ga_spherical_6r_,
                                                  varcont_6r_,
                                                  grid_5spot_,
                                                  logger_
  );

  while (!minimizer->IsFinished()) {
    auto next_case = minimizer->GetCaseForEvaluation();
    next_case->set_objf_value(abs(Sphere(next_case->GetRealVarVector())));
    minimizer->SubmitEvaluatedCase(next_case);
  }
  auto best_case = minimizer->GetTentativeBestCase();
//    EXPECT_NEAR(0.0, best_case->objective_function_value(), 0.1);
//    EXPECT_NEAR(0.0, best_case->GetRealVarVector()[0], 0.1);
//    EXPECT_NEAR(0.0, best_case->GetRealVarVector()[1], 0.1);
}

TEST_F(GeneticAlgorithmTest, TestFunctionRosenbrock) {
  test_case_ga_spherical_6r_->set_objf_value(abs(Rosenbrock(test_case_ga_spherical_6r_->GetRealVarVector())));
  Optimization::Optimizer *minimizer = new RGARDD(settings_ga_min_,
                                                  test_case_ga_spherical_6r_,
                                                  varcont_6r_,
                                                  grid_5spot_,
                                                  logger_
  );

  while (!minimizer->IsFinished()) {
    auto next_case = minimizer->GetCaseForEvaluation();
    next_case->set_objf_value(Rosenbrock(next_case->GetRealVarVector()));
    minimizer->SubmitEvaluatedCase(next_case);
  }
  auto best_case = minimizer->GetTentativeBestCase();

//    EXPECT_NEAR(3.72484, best_case->objective_function_value(), 1);
//    EXPECT_NEAR(1.0, best_case->GetRealVarVector()[0], 2.5);
//    EXPECT_NEAR(1.0, best_case->GetRealVarVector()[1], 2.5);
}

}

