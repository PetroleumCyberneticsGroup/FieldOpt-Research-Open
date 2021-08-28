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
  InitSecndOptmzr();
  InitializeBookkeeper();
  FinalizeInitialization(true);
}

void BilevelRunner::Execute() {
  while (optimizer_->IsFinished() == TC::NOT_FINISHED) {

    // get case from main optimizer
    Optimization::Case *new_case;
    new_case = optimizer_->GetCaseForEvaluation();

    // evaluate case
    if (bookkeeper_->IsEvaluated(new_case, true)) {
      if (vp_.vRUN >= 3) { ext_info("Bookkeeped case.", md_, cl_); }
      new_case->state.eval = ES::E_BOOKKEEPED;

    } else {
      try {
        bool sim_success = true;
        new_case->state.eval = ES::E_CURRENT;

        model_->ApplyCase(new_case);
        auto start = QDateTime::currentDateTime();
        simulator_->Evaluate();
        auto end = QDateTime::currentDateTime();
        int sim_time = time_span_seconds(start, end);

        // compute objective
        if (sim_success) {
          new_case->set_objf_value(objf_->value(false));

          new_case->state.eval = ES::E_DONE;
          new_case->SetSimTime(sim_time);
          sim_times_.push_back((sim_time));
        }

      } catch (runtime_error &e) {}
    }
    optimizer_->SubmitEvaluatedCase(new_case);
  }
  FinalizeRun(true);
}

}
