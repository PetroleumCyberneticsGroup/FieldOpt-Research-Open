/***********************************************************
Copyright (C) 2021
Mathias Bellout <chakibbb-pcg@gmail.com>

bellout - Sat Aug 28 2021 14:16:50 week 34 CET+0200

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
#include <dftr/DFTR.h>
#include "bilevel_runner.h"
#include "Utilities/printer.hpp"
#include "Model/model.h"

namespace Runner {

BilevelRunner::BilevelRunner(Runner::RuntimeSettings *runtime_settings)
  : AbstractRunner(runtime_settings) {
  InitializeLogger();
  InitializeSettings();
  InitializeModel();
  InitializeSimulator();
  EvaluateBaseModel();
  InitializeObjectiveFunction();
  InitializeBaseCase();
  InitializeOptimizer();
  InitializeBookkeeper();
  FinalizeInitialization(true);
}

void BilevelRunner::prntDbg(int lvl, Optimization::Optimizer *optmzr, Optimization::Case *cs) {

  string str;
  if (vp_.vOPT >= 3 && lvl == 1) {
    str = "[@upr_optmzr] |";
  } else if (vp_.vOPT >= 3 && lvl == 2) {
    str = "[@lwr_optmzr] |";
  }

  im_ = str;
  im_ += optmzr->GetStatusStringHeader().toStdString() + "|";
  im_ += optmzr->GetStatusString().toStdString();
  ext_info(im_, md_, cl_);

  im_ = str;
  im_ += cs->StringRepresentation(model_->variables());
  ext_info(im_, md_, cl_);

  im_ = str;
  im_ = "CaseHandler status: " + optmzr->case_handler()->status();
  im_ += " (QueuedCases()..empty()) => ";
  im_ += optmzr->case_handler()->QueuedCases().empty() ? "true" : "false";
  ext_info(im_, md_, cl_);

}

void BilevelRunner::Execute() {
  double fval_best = base_case_->objf_value();
  double bl_ps_fdiff = settings_->optimizer()->parameters().bl_ps_fdiff;
  im_ = "f diff in pattern threshold [fraction]: " + num2str(bl_ps_fdiff);
  ext_info(im_, md_, cl_);

  // bl_ps_fdiff = -std::numeric_limits<double>::infinity(); // override

  // note: ensemble functionality removed
  while (optimizer_->IsFinished() == TC::NOT_FINISHED) {

    // get case from main optimizer
    Optimization::Case *c_upper;
    c_upper = optimizer_->GetCaseForEvaluation();
    while (c_upper == nullptr) {
      if (optimizer_->IsFinished()) { break; } else {
        c_upper = optimizer_->GetCaseForEvaluation();
      }
    }
    if (optimizer_->IsFinished()) { break; }



    // run new case -> check if bookeeped
    if (bookkeeper_->IsEvaluated(c_upper, true)) {
      if (vp_.vRUN >= 3) { ext_info("Bookkeeped case.", md_, cl_); }
      c_upper->state.eval = ES::E_BOOKKEEPED;

    } else {
      bool upper_sim_success;
      c_upper->state.eval = ES::E_CURRENT;

      // compute fval for testcase0
      QDateTime upper_sim_start, upper_sim_end;
      model_->ApplyCase(c_upper);
      upper_sim_start = QDateTime::currentDateTime();
      upper_sim_success = simulator_->Evaluate(timeoutVal(), rts_->threads_per_sim());
      upper_sim_end = QDateTime::currentDateTime();
      int upper_sim_time = time_span_seconds(upper_sim_start, upper_sim_end);

      prntDbg(1,optimizer_, c_upper);

      if (upper_sim_success) {
        model_->wellCost(settings_->optimizer());
        if (settings_->optimizer()->objective().type == OT::Augmented) {
          c_upper->set_objf_value(objf_->value(false));
        } else {
          c_upper->set_objf_value(objf_->value());
        }

        string tm = "Objective function value set to ";
        tm += num2str(c_upper->objf_value(), 8, 1);
        if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }

        c_upper->state.eval = ES::E_DONE;
        c_upper->SetSimTime(upper_sim_time);
        sim_times_.push_back((upper_sim_time));

        // start of optmzr-lower-level scope ---------------
        if ((c_upper->objf_value() - fval_best)/fval_best > bl_ps_fdiff) {

          std::unique_ptr<Optimization::Optimizers::DFTR>
            optmzr_lwr_(new Optimization::Optimizers::DFTR(
            settings_->optimizer(), c_upper,
            model_->variables(), model_->grid(),
            log_lwr_, nullptr,
            model_->constraintHandler()));

          bool lower_sim_success;
          QDateTime lower_sim_start, lower_sim_end;
          Optimization::Case *c_lower = nullptr;

          while (optmzr_lwr_->IsFinished() == TC::NOT_FINISHED) {

            try {
              c_lower = optmzr_lwr_->GetCaseForEvaluation();
              c_lower->state.eval = ES::E_CURRENT;
              model_->ApplyCase(c_lower);

              prntDbg(2, optimizer_, c_upper);

              lower_sim_start = QDateTime::currentDateTime();
              lower_sim_success = simulator_->Evaluate(timeoutVal(), rts_->threads_per_sim());
              lower_sim_end = QDateTime::currentDateTime();
              int lower_sim_time = time_span_seconds(lower_sim_start, lower_sim_end);

              if (lower_sim_success) {
                c_lower->set_objf_value(objf_->value(false));
                c_lower->state.eval = ES::E_DONE;
                c_lower->SetSimTime(lower_sim_time);
                sim_times_.push_back((lower_sim_time));
              }

            } catch (runtime_error &e) {
              wm_ = "Exception while simulating case: " + string(e.what());
              ext_warn(wm_, md_, cl_);
            }

            optmzr_lwr_->SubmitEvaluatedCase(c_lower);
          }
          c_upper->CopyCaseVals(optmzr_lwr_->GetTentativeBestCase());
          // optimizer_->set_c_evald(optmzr_lwr_->get_c_evald());
          optimizer_->case_handler()->AddToNumberSimulated(optmzr_lwr_->get_c_evald());
        } // -------------------------------------------------

      } else { // set status for failed sim
        c_upper->set_objf_value(sentinelValue());
        c_upper->state.eval = ES::E_FAILED;
        c_upper->state.err_msg = EM::ERR_SIM;
        if (upper_sim_time >= timeoutVal()) {
          c_upper->state.eval = ES::E_TIMEOUT;
        }
      }
    }

    // main loop
    optimizer_->SubmitEvaluatedCase(c_upper);
    fval_best = optimizer_->GetTentativeBestCase()->objf_value();
  }
  FinalizeRun(true);
}

}
