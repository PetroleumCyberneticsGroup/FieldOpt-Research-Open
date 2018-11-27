//
// Created by bellout on 11/27/18.
//

#include <gtest/gtest.h>

// Define settings for tr-dfo
#include "Optimization/tests/test_resource_optimizer.h"
#include "Optimization/tests/test_resource_cases.h"

#include "Runner/tests/test_resource_runner.hpp"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"

#include "Optimization/optimizers/trust_region/TrustRegionOptimization.h"

#include "Utilities/math.hpp"
#include "Utilities/stringhelpers.hpp"

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;
using namespace std;

namespace {

class TrustRegionTest : public ::testing::Test,
                        public TestResources::TestResourceOptimizer,
                        public TestResources::TestResourceGrids {

 protected:
  TrustRegionTest() {
    base_ = base_case_;
  }

  virtual ~TrustRegionTest() {}
  virtual void SetUp() {}

  Optimization::Case *base_;

};

TEST_F(TrustRegionTest, Constructor) {

  cout << "Computing objective" << endl;
  test_case_tr_dfo_prob1_->set_objective_function_value(
      tr_dfo_prob1(test_case_tr_dfo_prob1_->GetRealVarVector()));

  auto tr_dfo = TrustRegionOptimization(
      settings_tr_opt_max_,
      test_case_tr_dfo_prob1_,
      varcont_tr_dfo_prob1_,
      grid_5spot_,
      logger_);

}

TEST_F(TrustRegionTest, tr_dfo_prob1) {

//  auto values = variables->GetContinousVariableValues();
//  VectorXd test_vec = VectorXd::Zero(values.size());
//  Case * init_case = new Case(base_case);
//  init_case->SetRealVarValues(test_vec);
//  case_handler_->AddNewCase(init_case);

  cout << "Computing objective" << endl;
  test_case_tr_dfo_prob1_->set_objective_function_value(
      tr_dfo_prob1(test_case_tr_dfo_prob1_->GetRealVarVector()));

  cout << "Setting up constructor, computing base case?" << endl;
  Optimization::Optimizer *tr_dfo_ = new TrustRegionOptimization(
      settings_tr_opt_max_,
      test_case_tr_dfo_prob1_,
      varcont_tr_dfo_prob1_,
      grid_5spot_,
      logger_);

  while (!tr_dfo_->IsFinished()) {
    cout << "Get case from optimizer" << endl;
    auto next_case = tr_dfo_->GetCaseForEvaluation();
    next_case->set_objective_function_value(tr_dfo_prob1(next_case->GetRealVarVector()));
    tr_dfo_->SubmitEvaluatedCase(next_case);
  }

  auto best_case = tr_dfo_->GetTentativeBestCase();
  cout << best_case->objective_function_value() << endl;

}

}
