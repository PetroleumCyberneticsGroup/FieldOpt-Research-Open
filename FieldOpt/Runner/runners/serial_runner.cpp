/***********************************************************
Copyright (C) 2015+2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2021 Mathias Bellout
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

#include <Utilities/time.hpp>
#include "serial_runner.h"
#include "Utilities/printer.hpp"
#include "Model/model.h"

namespace Runner {

SerialRunner::SerialRunner(Runner::RuntimeSettings *runtime_settings)
  : AbstractRunner(runtime_settings) {
  InitLogger();
  InitSettings();
  InitModel();
  InitSimulator();
  EvalBaseModel();
  InitObjF();
  InitBaseCase();
  InitializeOptimizer();
  InitBookkeeper();
  FinalizeInit(true);
}

void SerialRunner::Execute() {
  while (optimizer_->IsFinished() == TC::NOT_FINISHED) {

    Optimization::Case *new_case;

    // MAKE NEW CASE
    if (is_ensemble_run_) { // ENSEMBLE RUN
      if (ensemble_helper_.IsCaseDone()) {
        ensemble_helper_.SetActiveCase(optimizer_->GetCaseForEvaluation());
      }
      if (vp_.vRUN >= 3) { ext_info("Getting ensemble case.", md_, cl_); }

      new_case = ensemble_helper_.GetCaseForEval();
      auto en_rlz = new_case->GetEnsembleRlz().toStdString();
      model_->set_grid_path(ensemble_helper_.GetRlz(en_rlz).grid());

    } else { // SINGLE RUN
      if (vp_.vRUN >= 3) { ext_info("Getting case from Optimizer.", md_, cl_); }
      new_case = optimizer_->GetCaseForEvaluation();
      if (vp_.vRUN >= 3) { ext_info("Got case from Optimizer.", md_, cl_); }

      while (new_case == nullptr) {
        if (optimizer_->IsFinished()) { break;
        } else {
          new_case = optimizer_->GetCaseForEvaluation();
        }
      }
      if (optimizer_->IsFinished()) { break; }
    }

    // RUN NEW CASE -> CHECK IF BOOKKEEPED
    if (!is_ensemble_run_ && bookkeeper_->IsEvaluated(new_case, true)) {
      if (vp_.vRUN >= 3) { ext_info("Bookkeeped case.", md_, cl_); }
      new_case->state.eval = ES::E_BOOKKEEPED;

    } else { // RUN NEW CASE
      try {
        bool sim_success = true;
        new_case->state.eval = ES::E_CURRENT;

        if (vp_.vRUN >= 3) { ext_info("Applying case to model.", md_, cl_); }
        model_->ApplyCase(new_case);
        auto start = QDateTime::currentDateTime();

        // SIMULATE CASE (SINGLE RUN, NOT TIMEOUT)
        if (!is_ensemble_run_ && (sim_times_.empty() || rts_->sim_timeout() == 0)) {
          if (vp_.vRUN >= 3) { ext_info("Simulating case.", md_, cl_); }
          simulator_->Evaluate();

        } else { // SIMULATE CASE (ENSEMBLE RUN, WITH TIMEOUT)
          if (is_ensemble_run_) {
            if (vp_.vRUN >= 3) { ext_info("Simulating ensemble case.", md_, cl_); }
            auto en_rlz = ensemble_helper_.GetRlz(new_case->GetEnsembleRlz().toStdString());
            sim_success = simulator_->Evaluate(en_rlz, timeoutVal(), rts_->threads_per_sim()
            );

          } else { // SIMULATE CASE (SINGLE RUN, WITH TIMEOUT)
            if (vp_.vRUN >= 3) { ext_info("Simulating case.", md_, cl_); }
            sim_success = simulator_->Evaluate(timeoutVal(), rts_->threads_per_sim());
          }
        }
        if (vp_.vRUN >= 3) { ext_info("Done simulating case.", md_, cl_); }
        auto end = QDateTime::currentDateTime();
        int sim_time = time_span_seconds(start, end);
        if (sim_success) {

          model_->wellCost(settings_->optimizer());
          if (settings_->optimizer()->objective().type == OT::Augmented) {
            new_case->set_objf_value(objf_->value(false));
          } else {
            new_case->set_objf_value(objf_->value());
          }

          string tm = "Objective function value set to ";
          tm += num2str(new_case->objf_value(), 8, 1);
          if (vp_.vRUN >= 1) { info(tm,  vp_.lnw); }

          new_case->state.eval = ES::E_DONE;
          new_case->SetSimTime(sim_time);
          sim_times_.push_back((sim_time));

        } else {
          new_case->set_objf_value(sentinelValue());
          new_case->state.eval = ES::E_FAILED;
          new_case->state.err_msg = EM::ERR_SIM;
          if (sim_time >= timeoutVal()) {
            new_case->state.eval = ES::E_TIMEOUT;
          }
        }
      } catch (runtime_error &e) {
        wm_ = "Exception thrown while applying/simulating case: " + string(e.what());
        wm_ += ". Setting obj. fun. value to sentinel value.";
        ext_warn(wm_, md_, cl_);
        new_case->set_objf_value(sentinelValue());
        new_case->state.eval = ES::E_FAILED;
        new_case->state.err_msg = EM::ERR_WIC;
      }
    }
    if (is_ensemble_run_) {
      ensemble_helper_.SubmitEvaluatedRealization(new_case);
      if  (ensemble_helper_.IsCaseDone()) {
        optimizer_->SubmitEvaluatedCase(ensemble_helper_.GetEvaluatedCase());
      }
    }
    else {
      if (vp_.vRUN >= 3) { ext_info("Submitting evaluated case to Optimizer.", md_, cl_); }
      optimizer_->SubmitEvaluatedCase(new_case);
    }
  }
  FinalizeRun(true);
}

}
