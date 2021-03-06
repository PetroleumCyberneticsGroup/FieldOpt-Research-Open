/***********************************************************
Created by bellout on 08.03.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
test_dftr built from test_tr-dfo.cpp
Copyright (C) 2018
Mathias Bellout <chakibbb.pcg@gmail.com>
Thiago Lima Silva <thiagolims@gmail.com>
Caio Giuliani <caiogiuliani@gmail.com>
--

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
#include <objective/objective.h>

#include "Runner/tests/test_resource_runner.hpp"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Settings/tests/test_resource_settings.hpp"

#include "Optimization/optimizers/trust_region/TrustRegionOptimization.h"
#include "Optimization/optimizers/dftr/DFTR.h"
#include "Runner/runners/ensemble_helper.h"

#include "Utilities/printer.hpp"
#include "Utilities/colors.hpp"

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

class DFTRTest : public ::testing::Test,
                 public ::TestResourceOptimizer,
                 public ::TestResourceGrids {

 protected:
  DFTRTest() {}

  virtual ~DFTRTest() {}
  virtual void SetUp() {}

  TrustRegionOptimization *tr_dfo_;
  DFTR *dftr_;
  Optimization::Case *test_case_tr_dfo_probs_;
  VarPropContainer *varcont_tr_dfo_probs_;
  ::TrustRegionModelData tr_mdata;

  Settings::VerbParams vp_ = {};
  double eps = std::numeric_limits<double>::epsilon();

  Optimization::Optimizer::TerminationCondition TC_NOT_FIN_ =
    Optimization::Optimizer::TerminationCondition::NOT_FINISHED;

  Optimization::Optimizer::TerminationCondition TC_MAX_ITERS_ =
    Optimization::Optimizer::TerminationCondition::MAX_ITERS_REACHED;

  Optimization::Optimizer::TerminationCondition TC_MIN_STEP_ =
    Optimization::Optimizer::TerminationCondition::MIN_STEP_LENGTH_REACHED;

  Optimization::Optimizer::TerminationCondition TC_OPT_CRIT =
    Optimization::Optimizer::TerminationCondition::OPT_CRITERIA_REACHED;

  // _________________________________________________________
  // TRDFO/DFTR OPTIMIZERS
  void SetUpOptimizers(::TrustRegionModelData::prob &prob,
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

    settings_tr_opt_max_->setTRProbName(prob.name);
    settings_dftr_max_->setTRProbName(prob.name);

    tr_dfo_ = new TrustRegionOptimization(
      settings_tr_opt_max_, test_case_tr_dfo_probs_,
      varcont_tr_dfo_probs_, grid_5spot_, logger_);

    dftr_ = new DFTR(
      settings_dftr_max_, test_case_tr_dfo_probs_,
      varcont_tr_dfo_probs_, grid_5spot_, logger_);
  }

  bool RunnerSubs(TestResources::TrustRegionModelData::prob prob,
                  double (*tr_dfo_prob)(VectorXd xs)){
    stringstream ss; ss << "[          ] " << FMAGENTA;
    double tol = 1e-06;
    int p_count_dftr = 0, p_count_trdfo = 0;

    std::vector<VectorXd> vX_dftr, vX_trdfo, vXA, vXB;
    std::vector<double> vF_dftr, vF_trdfo, vFA, vFB;
    bool DX = false, DR = false, DF = false;
    std::vector<string> lbl = { "A", "B" };
    VectorXd data;

    // Start DFTR loop ------------------------------------
    while (dftr_->IsFinished() == TC_NOT_FIN_) {
      auto nxt_cs_dftr = dftr_->GetCaseForEvaluation();

      while (nxt_cs_dftr == nullptr) {
        if (dftr_->IsFinished()) { break; }
        else { nxt_cs_dftr = dftr_->GetCaseForEvaluation(); }
      }
      if (dftr_->IsFinished()) { break; }

      nxt_cs_dftr->set_objf_value(
        tr_dfo_prob(nxt_cs_dftr->GetRealVarVector()));

      // Override 2nd point
      if (p_count_dftr == 1 && prob.xm.cols() > 1) {
        TestResources::OverrideSecondPoint(prob, *nxt_cs_dftr);
      }

      // Finish Runner
      dftr_->SubmitEvaluatedCase(nxt_cs_dftr);
      p_count_dftr++;

      if (dftr_->getTRMod()->isModInit()) {
        data.conservativeResize(dftr_->getTRMod()->getCurrPt().rows()+1, 1);
        data << dftr_->getTRMod()->getCurrPt(), dftr_->getTRMod()->getRad();
        vX_dftr.push_back(data);
        vF_dftr.push_back(dftr_->getTRMod()->getCurrFval());
      }
    }

    // Start TRDFO loop ------------------------------------
    while (tr_dfo_->IsFinished() == TC_NOT_FIN_) {
      auto nxt_cs_dfo = tr_dfo_->GetCaseForEvaluation();

      while (nxt_cs_dfo == nullptr) {
        if (tr_dfo_->IsFinished()) { break; }
        else { nxt_cs_dfo = tr_dfo_->GetCaseForEvaluation(); }
      }
      if (tr_dfo_->IsFinished()) { break; }

      // Compute obj.function value for case
      nxt_cs_dfo->set_objf_value(
        tr_dfo_prob(nxt_cs_dfo->GetRealVarVector()));

      // Override 2nd point
      if (p_count_trdfo == 1 && prob.xm.cols() > 1) {
        TestResources::OverrideSecondPoint(prob, *nxt_cs_dfo);
      }

      // Finish Runner
      tr_dfo_->SubmitEvaluatedCase(nxt_cs_dfo);
      p_count_trdfo++;

      if (tr_dfo_->getTrustRegionModel()->isInitialized()) {
        data.conservativeResize(tr_dfo_->getTrustRegionModel()->getCurrentPoint().rows()+1, 1);
        data << tr_dfo_->getTrustRegionModel()->getCurrentPoint(), tr_dfo_->getTrustRegionModel()->getRadius();
        vX_trdfo.push_back(data);
        vF_trdfo.push_back(tr_dfo_->getTrustRegionModel()->getCurrentFval());
      }
    }

    // Start comparison ------------------------------------
    stringstream sx;
    string ccA, ccB, cc_trdfo, cc_dftr;

    if (tr_dfo_->IsFinished() == TC_OPT_CRIT) {
      cc_trdfo = "OPT_CRITERIA_REACHED";

    } else if (tr_dfo_->IsFinished() == TC_MIN_STEP_) {
      cc_trdfo = "MIN_STEP_LENGTH_REACHED";

    } else if (tr_dfo_->IsFinished() == TC_MAX_ITERS_) {
      cc_trdfo = "MAX_ITERS_REACHED";
    }

    if (dftr_->IsFinished() == TC_OPT_CRIT) {
      cc_dftr += "OPT_CRITERIA_REACHED";

    } else if (dftr_->IsFinished() == TC_MIN_STEP_) {
      cc_dftr += "MIN_STEP_LENGTH_REACHED";

    } else if (dftr_->IsFinished() == TC_MAX_ITERS_) {
      cc_dftr += "MAX_ITERS_REACHED";
    }

    if (vX_dftr.size() >= vX_trdfo.size()) {
      vXA = vX_dftr;
      vXB = vX_trdfo;
      vFA = vF_dftr;
      vFB = vF_trdfo;
      lbl = { "DFTR__", "TR-DFO" };
      ccA = cc_dftr;
      ccB = cc_trdfo;
    } else {
      vXB = vX_dftr;
      vXA = vX_trdfo;
      vFB = vF_dftr;
      vFA = vF_trdfo;
      lbl = { "TR-DFO", "DFTR__" };
      ccB = cc_dftr;
      ccA = cc_trdfo;
    }

    bool prntDbg = false;
    stringstream so;
    if (prntDbg) {
      so << FYELLOW << setw(12) << scientific << right << setprecision(6);
      so << "" << endl;
    }

    for (int ii=0; ii < vXA.size(); ii++) {
      if (prntDbg) {
        so << lbl[0] << " x = " << vXA[ii].transpose() << " | f = " << vFA[ii] << endl;
      }
      if(ii < vXB.size()) {
        DX = ((vXA[ii].head(vXA[ii].rows()-1) - vXB[ii].head(vXB[ii].rows()-1)).array().abs() < 2*eps).all();
        DR = ((vXA[ii].tail(1) - vXB[ii].tail(1)).array().abs() < 2*eps).all();
        DF = (abs(vFA[ii] - vFB[ii]) < 2*eps);

        if (DX && DR && DF && prntDbg) {
          so << lbl[1] << " x = " << vXB[ii].transpose() << " | f = " << vFB[ii];
          so << " -> DX=" << DX << " -- DF=" << DF << " -- DR=" << DR << endl;
        } else if (prntDbg) {
          so << FRED << lbl[1] << " x = " << vXB[ii].transpose() << " | f = " << vFB[ii] << AEND;
          if (!DX) { so << FRED << " -> DX=" << DX << AEND; } else { so << " -> DX=" << DX; }
          if (!DF) { so << FRED << " -- DF=" << DF << AEND; } else { so << " -> DF=" << DF; }
          if (!DR) { so << FRED << " -- DR=" << DR << AEND; } else { so << " -> DR=" << DR; }
          so << endl;
        }
      }
      so << AEND;
    }
    cout << so.str();

    EXPECT_TRUE(DX);
    EXPECT_TRUE(DF);

    sx << setw(12) << scientific << right << setprecision(6)
       << FYELLOW << "---------------------------------------------" << endl
       << lbl[0] << ": x* = " << vXA.back().transpose() << "  |  f* = " << vFA.back() << endl
       << "tc: " << ccA << endl
       << lbl[1] << ": x* = " << vXB.back().transpose() << "  |  f* = " << vFB.back() << endl
       << "tc: " << ccB << endl;
    sx << AEND;
    cout << sx.str();

    return true;
  }
};

TEST_F(DFTRTest, trHS1) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ HS1 ] "
       << " 100*pow(x(1) - pow(x(0), 2), 2) + pow(1 - x(0), 2)" << END << endl;

  SetUpOptimizers(tr_mdata.prob_hs1, hs1);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob_hs1, hs1));
}

TEST_F(DFTRTest, Prob1) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob1 ] "
       << " f = @(x) (1 - x(0))^2; x0=[-1.2 2.0]" << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob1, tr_dfo_prob1);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob1, tr_dfo_prob1));
}

TEST_F(DFTRTest, trDfoProb2) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob2 ] "
       << "f = @(x) log1p(x(0)^2) + x(1)^2; x0=[2.0 2.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob2, tr_dfo_prob2);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob2, tr_dfo_prob2));
}

TEST_F(DFTRTest, trDfoProb3) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob3 ] "
       << "f = @(x) sin(pi*x(0)/12) * cos(pi*x(1)/16); x0=[0.0 0.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob3, tr_dfo_prob3);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob3, tr_dfo_prob3));
}

TEST_F(DFTRTest, trDfoProb4) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob4 ] "
       << "f = @(x) 0.01*(x(1) - 1)^2 + (x(2) - x(1)^2)^2; x0=[2.0 2.0 2.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob4, tr_dfo_prob4);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob4, tr_dfo_prob4));
}

TEST_F(DFTRTest, trDfoProb5) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob5 ] "
       << "f = @(x) (x(0)-x(1))^2 + (x(1) - x(2))^4; x0=[-2.6 2.0 2.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob5, tr_dfo_prob5);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob5, tr_dfo_prob5));
}

TEST_F(DFTRTest, trDfoProb6) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob6 ] "
       << "f = @(x) (x(0) + x(1))^2 + (x(1) + x(2))^2; x0=[-4.0 1.0 1.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob6, tr_dfo_prob6);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob6, tr_dfo_prob6));
}

TEST_F(DFTRTest, trDfoProb7) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob7 ] "
       << "f = @(x) log1p(x(0)^2) + log1p((x(0) - x(1))^2) " << endl
       << " + log1p((x(1) - x(2))^2) + log1p((x(2) - x(3))^2); "
          "x0=[2.0 2.0 2.0 2.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob7, tr_dfo_prob7);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob7, tr_dfo_prob7));
}

TEST_F(DFTRTest, trDfoProb8) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob8 ] "
       << "f = @(x) (x(0)*x(1)*x(2)*x(3))^2; x0=[0.8 0.8 0.8 0.8]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob8, tr_dfo_prob8);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob8, tr_dfo_prob8));
}

TEST_F(DFTRTest, trDfoProb9) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob9 ] "
       << "f = @(x) (x(0)-1)^2 + (x(1)-2)^2 + (x(2)-3)^2 + (x(3)-4)^2; x0=[1.0 1.0 1.0 1.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob9, tr_dfo_prob9);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob9, tr_dfo_prob9));
}

TEST_F(DFTRTest, trDfoProb10) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob10 ] "
       << "f = @(x) (x(1) - x(2))^2 + (x(2) - x(3))^2 + " << endl
       << "(x(2) - x(3))^4 + (x(3) - x(4))^4; x0=[2.0 sqrt(2) -1.0 2-sqrt(2) 0.5]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob10, tr_dfo_prob10);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob10, tr_dfo_prob10));
}

TEST_F(DFTRTest, trDfoProb11) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob11 ] "
       << "f = @(x) sum(2*x./(x.*x + 1));; x0=[1.0 1.0 1.0 1.0]"
       << END << endl;

  // -------------------------------------------------------
  SetUpOptimizers(tr_mdata.prob11, tr_dfo_prob11);
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob11, tr_dfo_prob11));
}

}