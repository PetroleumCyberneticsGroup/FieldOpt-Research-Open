/***********************************************************
Copyright (C) 2017
Einar J.M. Baumann <einar.baumann@gmail.com>
Created by einar on 11/23/18.

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
#include <Runner/tests/test_resource_runner.hpp>
#include "Optimization/optimizers/PSO.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"

// #include "Optimization/tests/test_resource_cases.h"
#include "Model/tests/test_resource_model.h"

// #include <QtCore/QJsonDocument>

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;

namespace {

class PSOTest : public ::testing::Test,
                public TestResources::TestResourceOptimizer,
                public TestResources::TestResourceGrids
{
 protected:
  PSOTest() {
    base_ = base_case_;
  }
  virtual ~PSOTest() {}
  virtual void SetUp() {}

  Optimization::Case *base_;
  Optimization::Case* test_case_rstrt_;
};

TEST_F(PSOTest, Constructor) {
}

TEST_F(PSOTest, TestFunctionSpherical) {
//    test_case_2r_->set_objective_function_value(abs(Sphere(test_case_2r_->GetRealVarVector())));
//    test_case_ga_spherical_30r_->set_objective_function_value(abs(Sphere(test_case_ga_spherical_30r_->GetRealVarVector())));
  test_case_ga_spherical_6r_->set_objf_value(abs(Sphere(test_case_ga_spherical_6r_->GetRealVarVector())));
  Optimization::Optimizer *minimizer = new PSO(settings_pso_min_,
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
  EXPECT_NEAR(0.0, best_case->objf_value(), 0.12);
  EXPECT_NEAR(0.0, best_case->GetRealVarVector()[0], 0.12);
  EXPECT_NEAR(0.0, best_case->GetRealVarVector()[1], 0.12);
}

TEST_F(PSOTest, TestFunctionRosenbrock) {
  test_case_ga_spherical_6r_->set_objf_value(abs(Rosenbrock(test_case_ga_spherical_6r_->GetRealVarVector())));
  settings_pso_min_->SetRngSeed(5);
  Optimization::Optimizer *minimizer = new PSO(settings_pso_min_,
                                               test_case_ga_spherical_6r_,
                                               varcont_6r_,
                                               grid_5spot_,
                                               logger_);

  while (!minimizer->IsFinished()) {
    auto next_case = minimizer->GetCaseForEvaluation();
    next_case->set_objf_value(Rosenbrock(next_case->GetRealVarVector()));
    minimizer->SubmitEvaluatedCase(next_case);
  }
  auto best_case = minimizer->GetTentativeBestCase();

  EXPECT_NEAR(0.0, best_case->objf_value(), 1);
  EXPECT_NEAR(1.0, best_case->GetRealVarVector()[0], 0.5);
  EXPECT_NEAR(1.0, best_case->GetRealVarVector()[1], 0.5);
}

TEST_F(PSOTest, TestRestart) {

  test_case_rstrt_ = new Optimization::Case(model_rstrt_->variables()->GetBinVarValues(),
                                            model_rstrt_->variables()->GetDiscVarValues(),
                                            model_rstrt_->variables()->GetContVarValues());
  test_case_rstrt_->set_objf_value(0.0);

  test_case_rstrt_->StringRepresentation(model_rstrt_->variables());



  Optimization::Optimizer *optmzr = new PSO(settings_rstrt_opt_,
                                            test_case_rstrt_,
                                            model_rstrt_->variables(),
                                            grid_olympr37_,
                                            logger_);

  for ( auto *c : optmzr->case_handler()->AllCases()) {

  }
}

}
