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
  Optimization::Case* base_case_rstrt_;
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
                                               logger_);

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

  auto vars = model_rstrt_->variables();
  base_case_rstrt_ = new Optimization::Case(vars->GetBinVarValues(),
                                            vars->GetDiscVarValues(),
                                            vars->GetContVarValues());
  base_case_rstrt_->set_objf_value(0.0);

  base_case_rstrt_->StringRepresentation(vars);

  Optimization::Optimizer *optmzr = new PSO(settings_rstrt_opt_,
                                            base_case_rstrt_,
                                            vars,
                                            grid_olympr37_,
                                            logger_);

  map<string, map<string, double>> vmap;
  vmap.insert(
      {"bc", {
          { "Var#BHP#INJD-15#0", 199.96961185494197 },
          { "Var#BHP#INJD-15#1094", 209.54137953870247 },
          { "Var#BHP#INJD-15#1641", 250.0531920881195 },
          { "Var#BHP#INJD-15#2188", 197.0584663055752 },
          { "Var#BHP#INJD-15#2735", 200.6145349062493 },
          { "Var#BHP#INJD-15#3282", 195.43382947558968 },
          { "Var#BHP#INJD-15#3829", 253.89998863950598 },
          { "Var#BHP#INJD-15#4376", 218.0791077686635 },
          { "Var#BHP#INJD-15#4923", 250.38103437442064 },
          { "Var#BHP#INJD-15#547", 244.80148597015165 },
          { "Var#BHP#INJD-16#0", 236.31054225077148 },
          { "Var#BHP#INJD-16#1094", 238.83416476794562 },
          { "Var#BHP#INJD-16#1641", 201.13418297135198 },
          { "Var#BHP#INJD-16#2188", 239.78391825977764 },
          { "Var#BHP#INJD-16#2735", 254.14934644255555 },
          { "Var#BHP#INJD-16#3282", 228.72663503476676 },
          { "Var#BHP#INJD-16#3829", 195.01219735780882 },
          { "Var#BHP#INJD-16#4376", 247.82345714823884 },
          { "Var#BHP#INJD-16#4923", 202.85386597435658 },
          { "Var#BHP#INJD-16#547", 253.05227067702344 },
          { "Var#BHP#PRODX2#0", 167.36580743563388 },
          { "Var#BHP#PRODX2#1094", 105.43188190973554 },
          { "Var#BHP#PRODX2#1641", 117.34227644118684 },
          { "Var#BHP#PRODX2#2188", 157.77387174276436 },
          { "Var#BHP#PRODX2#2735", 173.02961678205642 },
          { "Var#BHP#PRODX2#3282", 106.125474949253 },
          { "Var#BHP#PRODX2#3829", 109.99697898592774 },
          { "Var#BHP#PRODX2#4376", 105.4197999414107 },
          { "Var#BHP#PRODX2#4923", 106.73833844615724 },
          { "Var#BHP#PRODX2#547", 105.00074442714711 },
          { "Var#ICD#PRODX2#0", 0.00653218586329762 },
          { "Var#ICD#PRODX2#1", 0.0001577555527071811 },
          { "Var#ICD#PRODX2#2", 0.00012941330135794304 },
          { "Var#ICD#PRODX2#3", 3.9575624016865985e-06 },
          { "Var#ICD#PRODX2#4", 1.2178067233547564e-05 },
          { "Var#ICD#PRODX2#5", 6.558205015438912e-07 },
          { "Var#ICD#PRODX2#6", 5.045821649893485e-06 },
          { "Var#ICD#PRODX2#7", 1.127772782414261e-05 }, }
      });

  EXPECT_NEAR(0.7571595999999999, optmzr->GetTentativeBestCase()->objf_value(), 1e-12);

  // check base case
  for ( auto v : vmap.at("bc")) {
    auto uuid = vars->GetContVarQUuid(QString::fromStdString(v.first).replace("Var#",""));
    auto val = optmzr->GetTentativeBestCase()->get_real_variable_value(uuid);
    EXPECT_NEAR(val, v.second, 1e-12);
  }


  // for ( auto *c : optmzr->case_handler()->AllCases()) {
  //   for ( auto v : bc_vmap) {
  //     c->GetValues().at("")
  //   }
  // }
}

}
