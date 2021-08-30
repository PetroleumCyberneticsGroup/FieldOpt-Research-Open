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
  // InitSecndOptmzr(); // <-

  InitializeBookkeeper();
  FinalizeInitialization(true);
}

void BilevelRunner::Execute() {
  while (optimizer_->IsFinished() == TC::NOT_FINISHED) {

    // get case from main optimizer
    Optimization::Case *c_upper;
    c_upper = optimizer_->GetCaseForEvaluation();

    {
      std::unique_ptr<Optimization::Optimizers::DFTR>
        secndoptmzr_(new Optimization::Optimizers::DFTR(
        settings_->optimizer(), c_upper,
        model_->variables(),model_->grid(),
        logger_,nullptr,
        model_->constraintHandler()));

      // secndoptmzr_->updateTentativeBestCase(c_upper);

      bool lower_sim_success;
      QDateTime lower_sim_start, lower_sim_end;
      Optimization::Case *c_lower = nullptr;

      while (secndoptmzr_->IsFinished() == TC::NOT_FINISHED) {

        try {
          c_lower = secndoptmzr_->GetCaseForEvaluation();
          c_lower->state.eval = ES::E_CURRENT;
          model_->ApplyCase(c_lower);

          lower_sim_start = QDateTime::currentDateTime();
          lower_sim_success = simulator_->Evaluate(timeoutVal(), rts_->threads_per_sim());
          lower_sim_end = QDateTime::currentDateTime();
          int lower_sim_time = time_span_seconds(lower_sim_start, lower_sim_end);

          if(lower_sim_success) {
            c_lower->set_objf_value(objf_->value(false));
            c_lower->state.eval = ES::E_DONE;
            c_lower->SetSimTime(lower_sim_time);
            sim_times_.push_back((lower_sim_time));
          }

        } catch (runtime_error &e) {
          wm_ = "Exception while simulating case: " + string(e.what());
          ext_warn(wm_, md_, cl_);
        }

        secndoptmzr_->SubmitEvaluatedCase(c_lower);
      }
      c_upper->CopyCaseVals(secndoptmzr_->GetTentativeBestCase());
    }

    optimizer_->SubmitEvaluatedCase(c_upper);
    // secndoptmzr_->ResetOptimizer();



    // if (bookkeeper_->IsEvaluated(c_upper, true)) {
    //   if (vp_.vRUN >= 3) { ext_info("Bookkeeped case.", md_, cl_); }
    //   c_upper->state.eval = ES::E_BOOKKEEPED;
    //
    // } else {
    //   try {
    //     bool sim_success = true;
    //     c_upper->state.eval = ES::E_CURRENT;
    //
    //     model_->ApplyCase(c_upper);
    //     auto start = QDateTime::currentDateTime();
    //     simulator_->Evaluate();
    //
    //     auto end = QDateTime::currentDateTime();
    //     int sim_time = time_span_seconds(start, end);
    //
    //     // compute objective
    //     if (sim_success) {
    //       c_upper->set_objf_value(objf_->value(false));
    //
    //       c_upper->state.eval = ES::E_DONE;
    //       c_upper->SetSimTime(sim_time);
    //       sim_times_.push_back((sim_time));
    //     }
    //
    //   } catch (runtime_error &e) {}
    // }
    // optimizer_->SubmitEvaluatedCase(c_upper);
  }
  FinalizeRun(true);
}

}
