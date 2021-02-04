/***********************************************************
Created by bellout on 5/7/20.
Copyright (C) 2020- Mathias Bellout
<chakibbb-pcg@gmail.com>

Modified 2020 Amanda Machado
<amanda.automacaoufsc@gmail.com>

Modified 2020 Thiago Lima Silva
<thiagolims@gmail.com>

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

#include "trust_region/TrustRegionOptimization.h"
#include "Runner/runners/ensemble_helper.h"
#include "Simulation/simulator_interfaces/eclsimulator.h"
#include "Optimization/objective/weightedsum.h"

#include "Utilities/math.hpp"
#include "Utilities/printer.hpp"
#include "Utilities/debug.hpp"
#include "Utilities/colors.hpp"
#include "Utilities/stringhelpers.hpp"

#include "test_tr-model-data.hpp"
#include "test_tr-support.hpp"

// namespace TestResources {
//   void PrintCaseData(Optimization::Case &c);
// }

// ---------------------------------------------------------
using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;
using namespace Optimization::Objective;
using namespace TestResources;
using namespace Runner;
using namespace Simulation;
using namespace std;

namespace {

// =========================================================
class EnTRTest : public ::testing::Test,
                 public ::TestResourceOptimizer,
                 public ::TestResourceGrids {

 protected:
  EnTRTest() {}

  virtual ~EnTRTest() {}
  virtual void SetUp() {}

  // -------------------------------------------------------
  TrustRegionOptimization *tr_dfo_en_5spot_;
  Runner::EnsembleHelper ensemble_helper_;
  Simulation::Simulator *simulator_en_5spot_;
  Objective *objective_function_;
  Optimization::Case *base_case_en_5spot_;

  std::vector<int> simulation_times_;
  const double sentinel_value_ = 0.0001;
  const int timeout_value_ = 10000;

  Optimization::Optimizer::TerminationCondition TC_NOT_FIN_ =
      Optimization::Optimizer::TerminationCondition::NOT_FINISHED;

  Optimization::Optimizer::TerminationCondition TC_MAX_ITERS_ =
      Optimization::Optimizer::TerminationCondition::MAX_ITERATIONS_REACHED;

  Optimization::Optimizer::TerminationCondition TC_MIN_STEP_ =
      Optimization::Optimizer::TerminationCondition::MINIMUM_STEP_LENGTH_REACHED;

  Optimization::Optimizer::TerminationCondition TC_OPT_CRIT =
      Optimization::Optimizer::TerminationCondition::OPTIMALITY_CRITERIA_REACHED;

  Optimization::Case::CaseState::EvalStatus CS_EVAL_CURRENT_ =
      Optimization::Case::CaseState::EvalStatus::E_CURRENT;

  Optimization::Case::CaseState::EvalStatus CS_EVAL_DONE_ =
      Optimization::Case::CaseState::EvalStatus::E_DONE;

  Optimization::Case::CaseState::EvalStatus CS_EVAL_FAILED_ =
      Optimization::Case::CaseState::EvalStatus::E_FAILED;

  Optimization::Case::CaseState::EvalStatus CS_EVAL_TIMEOUT_ =
      Optimization::Case::CaseState::EvalStatus::E_TIMEOUT;

  Optimization::Case::CaseState::ErrorMessage CS_ERR_MSG_SIM_ =
      Optimization::Case::CaseState::ErrorMessage::ERR_SIM;

  Optimization::Case::CaseState::ErrorMessage CS_ERR_MSG_WIC_ =
      Optimization::Case::CaseState::ErrorMessage::ERR_WIC;

  string DBG_fn_ensrnr_ ;
  TrustRegionModel *tr_dfo_mod_;
  string SD60_ = string(60, '-');
  string DD60_ = string(60, '=');

  // =======================================================
  void InitRunnerSubs() {

    settings_en_5spot_full_->paths().ShowPaths();

    // -----------------------------------------------------
    ensemble_helper_ = EnsembleHelper(
        settings_en_5spot_full_->simulator()->get_ensemble(),
        settings_en_5spot_full_->optimizer()->parameters().rng_seed);

    // -----------------------------------------------------
    simulator_en_5spot_ = new ECLSimulator(settings_en_5spot_full_,
                                           model_en_5spot_);

    // -----------------------------------------------------
    simulator_en_5spot_->Evaluate(
        ensemble_helper_.GetBaseRealization(),
        timeout_value_, 1);

    // -----------------------------------------------------
    objective_function_ = new WeightedSum(
        settings_en_5spot_full_->optimizer(),
        simulator_en_5spot_->results(),
        model_en_5spot_);

    // -----------------------------------------------------
    base_case_en_5spot_ = new Optimization::Case(
        model_en_5spot_->variables()->GetBinaryVariableValues(),
        model_en_5spot_->variables()->GetDiscreteVariableValues(),
        model_en_5spot_->variables()->GetContinousVariableValues());

    model_->wellCost(settings_en_5spot_full_->optimizer());
    base_case_en_5spot_->set_objective_function_value(
        objective_function_->value());

    // -----------------------------------------------------
    auto opt_par = settings_en_5spot_full_->optimizer()->parameters();
    DBG_fn_ensrnr_ = "FOEx_EnsRnr_" + opt_par.tr_prob_name + ".txt";

  }

  // =======================================================
  void SetUpOptimizer() {
    tr_dfo_en_5spot_ = new TrustRegionOptimization(
        settings_en_5spot_full_->optimizer(),
        base_case_en_5spot_,
        model_en_5spot_->variables(),
        model_en_5spot_->grid(),
        logger_);

    // dbg
    tr_dfo_mod_ = tr_dfo_en_5spot_->getTrustRegionModel();
  }

  // =======================================================
  void DBG_prntSoln(stringstream &sx) {

    // stringstream sx;
    string cc;

    // ---------------------------------------------------
    if (tr_dfo_en_5spot_->IsFinished() == TC_OPT_CRIT) {
      cc = "OPTIMALITY_CRITERIA_REACHED";

    } else if (tr_dfo_en_5spot_->IsFinished() == TC_MIN_STEP_) {
      cc = "MINIMUM_STEP_LENGTH_REACHED";

    } else if (tr_dfo_en_5spot_->IsFinished() == TC_MAX_ITERS_) {
      cc = "MAX_ITERATIONS_REACHED";
    }

    sx << setw(12) << scientific << right << setprecision(6)
       << "---------------------------------------------" << endl
       << "x* = " << tr_dfo_mod_->getCurrentPoint().transpose() << endl
       << "f* = " << tr_dfo_mod_->getCurrentFval() << endl
       << "tc: " << cc.c_str();
    // cout << sx.str();

  }

  // =======================================================
  bool RunnerSubs() {

    // -----------------------------------------------------
    stringstream ss; ss << "[          ] " << FMAGENTA;
    stringstream sd; // dbg
    double tol = 1e-06;
    int p_count = 0;
    bool is_ensemble_run_ = true;

    // -----------------------------------------------------
    while (tr_dfo_en_5spot_->IsFinished() == TC_NOT_FIN_) {
      sd.str("");
      sd << SD60_ << " p_count=" << p_count  << endl; // dbg

      Optimization::Case *next_case;

        // -------------------------------------------------
        // RESET LIST
        if (ensemble_helper_.IsCaseDone()) {

          sd << "rlz list for ens.case *EMPTY* -- "
          << "get new case from optimizer"<< endl;
          next_case = tr_dfo_en_5spot_->GetCaseForEvaluation();
          while (next_case == nullptr) {
            if (tr_dfo_en_5spot_->IsFinished()) {
              break;
            } else {
              next_case = tr_dfo_en_5spot_->GetCaseForEvaluation();
            }
          }

          if (tr_dfo_en_5spot_->IsFinished()) {
            break;
          }

          ensemble_helper_.SetActiveCase(next_case);

        } else {
          sd << "rlz list for ens.case *NOT YET EMPTY*" << endl;
        }

        if (VERB_RUN >=2) {
          TestResources::PrintCaseData(*next_case); // dbg
        }

        // -------------------------------------------------
        // CONTINUE LIST -> NEXT RLZ
        // Making temp case for particular rlz [pass ref.:]
        // next_case = new Optimization::Case(base_case_en_5spot_);
        ensemble_helper_.GetCaseForEval(*next_case);
        // [using old.method:]
        // next_case = ensemble_helper_.GetCaseForEval();

        // -------------------------------------------------
        // APPLY RLZ
        string en_rlz = next_case->GetEnsembleRealization().toStdString();
        sd << "new rlz tag =: " << en_rlz << endl; // dbg

        // Update grid path for model from current en_rlz
        model_en_5spot_->set_grid_path(
            ensemble_helper_.GetRealization(en_rlz).grid());

      // Skipped bookkeeper process
      // ---------------------------------------------------
      // REGULAR COMPUTATION OF CASE
      try {
        bool simulation_success = true;
        next_case->state.eval = CS_EVAL_CURRENT_;

        // TestResources::PrintCaseData(*next_case); // dbg
        model_en_5spot_->ApplyCase(next_case);

        auto start = QDateTime::currentDateTime();

        auto nxt_rlz = ensemble_helper_.GetRealization(
            next_case->GetEnsembleRealization().toStdString());

        simulation_success = simulator_en_5spot_->Evaluate(
            nxt_rlz, timeout_value_, rts_en_5spot_->threads_per_sim());
        sd << "[1] " << ensemble_helper_.GetStateString() << endl; // dbg

        // -----------------------------------------------
        auto end = QDateTime::currentDateTime();
        int sim_time = time_span_seconds(start, end);

        // -----------------------------------------------
        if ( simulation_success ) {
          model_en_5spot_->wellCost(
              settings_en_5spot_full_->optimizer());

          next_case->set_objective_function_value(
              objective_function_->value());
          sd << "objf.r" << nxt_rlz.alias() << ": "
          << next_case->objective_function_value() << endl; // dbg

          next_case->state.eval = CS_EVAL_DONE_;

          next_case->SetSimTime(sim_time);
          simulation_times_.push_back((sim_time));

          } else {
            next_case->set_objective_function_value(sentinel_value_);

            next_case->state.eval = CS_EVAL_FAILED_;

            next_case->state.err_msg = CS_ERR_MSG_SIM_;

            if ( sim_time >= timeout_value_ ) {
              next_case->state.eval = CS_EVAL_TIMEOUT_;
            }
        }

        // -------------------------------------------------
      } catch (std::runtime_error e) {
        Printer::ext_warn(
            "Exception thrown while applying/simulating case: "
                + std::string(e.what()) + ". Setting fval to sentinel value.",
            "dc_1Runner", "SerialRunner");

        // Reapply later:
        // next_case->set_objective_function_value(sentinel_value_);
        next_case->state.eval = CS_EVAL_FAILED_;
        next_case->state.err_msg = CS_ERR_MSG_WIC_;
      }

      // ---------------------------------------------------
      if ( is_ensemble_run_ ) {
        ensemble_helper_.SubmitEvaluatedRealization(next_case);

        if ( ensemble_helper_.IsCaseDone()) {
          tr_dfo_en_5spot_->SubmitEvaluatedCase(
              ensemble_helper_.GetEvaluatedCase());

          sd << DD60_ << "\n[2] "
          << ensemble_helper_.GetStateString() << "\nobjf.mean: "
          << ensemble_helper_.GetEvaluatedCase()->GetEnsembleAverageOfv()
          << endl << endl; // dbg
        }
      } else {
        tr_dfo_en_5spot_->SubmitEvaluatedCase(next_case);
      }

      p_count++;

      // ---------------------------------------------------
      // dbg
      tr_dfo_mod_->DBG_printToFile(DBG_fn_ensrnr_, sd.str());
    }

    // FinalizeRun
    model_en_5spot_->ApplyCase(tr_dfo_en_5spot_->GetTentativeBestCase());
    simulator_en_5spot_->WriteDriverFilesOnly();

    // PrintCompletionMessage();
    model_en_5spot_->Finalize();
    logger_->FinalizePostrunSummary();

    // dbg
    DBG_prntSoln(sd);
    tr_dfo_mod_->DBG_printToFile(DBG_fn_ensrnr_, sd.str());
  }

};

// =========================================================
TEST_F(EnTRTest, en_5spot_prob1) {
  InitRunnerSubs();
  SetUpOptimizer();
  RunnerSubs();
}


}
