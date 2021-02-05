/***********************************************************
Created by Mathias C. Bellout on 27.11.18
Copyright (C) 2018-2020 Mathias Bellout
<chakibbb.pcg@gmail.com>

Modified 2018-2020 Thiago Lima Silva
<thiagolims@gmail.com>

Modified 2018-2020 Caio Giuliani
<caiogiuliani@gmail.com>

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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <objective/objective.h>

#include "Runner/tests/test_resource_runner.hpp"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Optimization/tests/test_resource_cases.h"
#include "Settings/tests/test_resource_settings.hpp"

#include "Optimization/optimizers/trust_region/TrustRegionOptimization.h"
#include "Runner/runners/ensemble_helper.h"
#include "Simulation/simulator_interfaces/eclsimulator.h"
#include "Optimization/objective/weightedsum.h"

#include "Utilities/math.hpp"
#include "Utilities/printer.hpp"
#include "Utilities/debug.hpp"
#include "Utilities/colors.hpp"
// #include "Utilities/stringhelpers.hpp"

#include "test_tr-model-data.hpp"
#include "test_tr-support.hpp"

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;
using namespace Optimization::Objective;
using namespace TestResources;
using namespace Runner;
using namespace Simulation;
using namespace std;

namespace {

class TrustRegionTest : public ::testing::Test,
                        public ::TestResourceOptimizer,
                        public ::TestResourceGrids {

 protected:
  TrustRegionTest() {}

  virtual ~TrustRegionTest() {}
  virtual void SetUp() {}

  TrustRegionOptimization *tr_dfo_;
  Optimization::Case *test_case_tr_dfo_probs_;
  VarPropContainer *varcont_tr_dfo_probs_;
  ::TrustRegionModelData tr_mdata;

  Settings::VerbParams vp_ = {};

  Optimization::Optimizer::TerminationCondition TC_NOT_FIN_ =
      Optimization::Optimizer::TerminationCondition::NOT_FINISHED;

  Optimization::Optimizer::TerminationCondition TC_MAX_ITERS_ =
      Optimization::Optimizer::TerminationCondition::MAX_ITERATIONS_REACHED;

  Optimization::Optimizer::TerminationCondition TC_MIN_STEP_ =
      Optimization::Optimizer::TerminationCondition::MINIMUM_STEP_LENGTH_REACHED;

  Optimization::Optimizer::TerminationCondition TC_OPT_CRIT =
      Optimization::Optimizer::TerminationCondition::OPTIMALITY_CRITERIA_REACHED;

  void SetUpOptimizer(::TrustRegionModelData::prob &prob,
                      double (*tr_dfo_prob)(VectorXd xs)) {
    VectorXd x0 = prob.xm.col(0);

    // Dummy var container based on initial point
    varcont_tr_dfo_probs_ = new VarPropContainer(vp_);
    QString base_varname = "BHP#PRODUCER#"; // dummy var name

    for (int i = 0; i < x0.rows(); ++i) {
      // Use initial point values to construct container
      auto *prop = new ContinuousProperty(x0(i));
      prop->setName(base_varname + QString::number(i));
      varcont_tr_dfo_probs_->AddVariable(prop);
    }

    // Set up base case using dummy var containter
    test_case_tr_dfo_probs_ = new Optimization::Case(
        QList<QPair<QUuid, bool>>(),
        QList<QPair<QUuid, int>>(),
        varcont_tr_dfo_probs_->GetContVarValues());
    TestResources::FindVarSequence(prob,
                                   *test_case_tr_dfo_probs_);

    VectorXd ordered_vec(test_case_tr_dfo_probs_->GetRealVarVector().size());
    auto vec =test_case_tr_dfo_probs_->GetRealVarVector();
    for (int i = 0; i < prob.idx.size(); i++) {
      ordered_vec(prob.idx[i]) = vec(i);
    }
    test_case_tr_dfo_probs_->SetRealVarValues(ordered_vec);

    // Use initial point (objf) from Matlab data
    test_case_tr_dfo_probs_->set_objf_value(tr_dfo_prob(x0));

    tr_dfo_ = new TrustRegionOptimization(
        settings_tr_opt_max_,
        test_case_tr_dfo_probs_,
        varcont_tr_dfo_probs_,
        grid_5spot_,
        logger_);
  }

  bool RunnerSubs(TestResources::TrustRegionModelData::prob prob,
                  double (*tr_dfo_prob)(VectorXd xs)){
    stringstream ss; ss << "[          ] " << FMAGENTA;
    double tol = 1e-06;
    int p_count = 0;

    // Start opt loop --------------------------------------
    while (tr_dfo_->IsFinished() == TC_NOT_FIN_) {

      auto next_case = tr_dfo_->GetCaseForEvaluation();
      while (next_case == nullptr) {
        if (tr_dfo_->IsFinished()) {
          break;
        } else {
          next_case = tr_dfo_->GetCaseForEvaluation();
        }
      }

      if (tr_dfo_->IsFinished()) {
        break;
      }

      // Compute obj.function value for case
      next_case->set_objf_value(
        tr_dfo_prob(next_case->GetRealVarVector()));

      // Override 2nd point
      if (p_count == 1 && prob.xm.cols() > 1) {
        TestResources::OverrideSecondPoint(prob, *next_case);
      }

      // Finish Runner
      tr_dfo_->SubmitEvaluatedCase(next_case);

      if (tr_dfo_->getTrustRegionModel()->isInitialized() && (p_count < 2)) {

      }
      p_count++;
    }

    stringstream sx;
    string cc;

    if (tr_dfo_->IsFinished() == TC_OPT_CRIT) {
      cc = "OPTIMALITY_CRITERIA_REACHED";

    } else if (tr_dfo_->IsFinished() == TC_MIN_STEP_) {
      cc = "MINIMUM_STEP_LENGTH_REACHED";

    } else if (tr_dfo_->IsFinished() == TC_MAX_ITERS_) {
      cc = "MAX_ITERATIONS_REACHED";
    }

    sx << setw(12) << scientific << right << setprecision(6)
       << "---------------------------------------------" << endl
       << "x* = " << tr_dfo_->getTrustRegionModel()->getCurrentPoint().transpose() << endl
       << "f* = " << tr_dfo_->getTrustRegionModel()->getCurrentFval() << endl
       << "tc: " << cc.c_str() << endl;
    cout << sx.str();

    return true;
  }
};

TEST_F(TrustRegionTest, trHS1) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
         << "[ HS1 ] "
         << " 100*pow(x(1) - pow(x(0), 2), 2) + pow(1 - x(0), 2)" << END << endl;

    SetUpOptimizer(tr_mdata.prob_hs1, hs1);
    EXPECT_TRUE(RunnerSubs(tr_mdata.prob_hs1, hs1));
}

TEST_F(TrustRegionTest, trDfoProb1) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob1 ] "
       << " f = @(x) (1 - x(0))^2; x0=[-1.2 2.0]" << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob1, tr_dfo_prob1);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob1, tr_dfo_prob1));
}

TEST_F(TrustRegionTest, trDfoProb2) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob2 ] "
       << "f = @(x) log1p(x(0)^2) + x(1)^2; x0=[2.0 2.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob2, tr_dfo_prob2);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob2, tr_dfo_prob2));
}

TEST_F(TrustRegionTest, trDfoProb3) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob3 ] "
       << "f = @(x) sin(pi*x(0)/12) * cos(pi*x(1)/16); x0=[0.0 0.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob3, tr_dfo_prob3);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob3, tr_dfo_prob3));
}

TEST_F(TrustRegionTest, trDfoProb4) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob4 ] "
       << "f = @(x) 0.01*(x(1) - 1)^2 + (x(2) - x(1)^2)^2; x0=[2.0 2.0 2.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob4, tr_dfo_prob4);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob4, tr_dfo_prob4));
}

TEST_F(TrustRegionTest, trDfoProb5) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob5 ] "
       << "f = @(x) (x(0)-x(1))^2 + (x(1) - x(2))^4; x0=[-2.6 2.0 2.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob5, tr_dfo_prob5);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob5, tr_dfo_prob5));
}

TEST_F(TrustRegionTest, trDfoProb6) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob6 ] "
       << "f = @(x) (x(0) + x(1))^2 + (x(1) + x(2))^2; x0=[-4.0 1.0 1.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob6, tr_dfo_prob6);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob6, tr_dfo_prob6));
}

TEST_F(TrustRegionTest, trDfoProb7) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob7 ] "
       << "f = @(x) log1p(x(0)^2) + log1p((x(0) - x(1))^2) " << endl
       << " + log1p((x(1) - x(2))^2) + log1p((x(2) - x(3))^2); "
          "x0=[2.0 2.0 2.0 2.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob7, tr_dfo_prob7);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob7, tr_dfo_prob7));
}

TEST_F(TrustRegionTest, trDfoProb8) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob8 ] "
       << "f = @(x) (x(0)*x(1)*x(2)*x(3))^2; x0=[0.8 0.8 0.8 0.8]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob8, tr_dfo_prob8);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob8, tr_dfo_prob8));
}

TEST_F(TrustRegionTest, trDfoProb9) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob9 ] "
       << "f = @(x) (x(0)-1)^2 + (x(1)-2)^2 + (x(2)-3)^2 + (x(3)-4)^2; x0=[1.0 1.0 1.0 1.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob9, tr_dfo_prob9);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob9, tr_dfo_prob9));
}

TEST_F(TrustRegionTest, trDfoProb10) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob10 ] "
       << "f = @(x) (x(1) - x(2))^2 + (x(2) - x(3))^2 + " << endl
       << "(x(2) - x(3))^4 + (x(3) - x(4))^4; x0=[2.0 sqrt(2) -1.0 2-sqrt(2) 0.5]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob10, tr_dfo_prob10);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob10, tr_dfo_prob10));
}

TEST_F(TrustRegionTest, trDfoProb11) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob11 ] "
       << "f = @(x) sum(2*x./(x.*x + 1));; x0=[1.0 1.0 1.0 1.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizer(tr_mdata.prob11, tr_dfo_prob11);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob11, tr_dfo_prob11));
}

}
