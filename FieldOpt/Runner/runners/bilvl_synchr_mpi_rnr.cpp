/***********************************************************
Copyright (C) 2021
Mathias Bellout <chakibbb-pcg@gmail.com>

bellout - Tue Aug 31 2021 15:22:17 week 35 CET+0200

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

#include <dftr/DFTR.h>
#include "bilvl_synchr_mpi_rnr.h"

namespace Runner {
namespace MPI {

BilevelSynchrMPIRunner::BilevelSynchrMPIRunner(RuntimeSettings *rts) : MPIRunner(rts) {
  assert(world_.size() >= 2 && "BilevelSynchrMPIRunner requires at least two MPI processes.");

  if (world_.rank() == 0) {
    InitializeSettings("rank" + QString::number(rank()));
    InitializeLogger();
    InitializeModel();
    InitializeSimulator();
    EvaluateBaseModel();
    InitializeObjectiveFunction();
    InitializeBaseCase();
    InitializeOptimizer();
    InitializeBookkeeper();
    overseer_ = new MPI::Overseer(this);
    overseer_->addToBestF(base_case_->objf_value());
    FinalizeInitialization(true);

  } else {
    InitializeLogger("rank" + QString::number(rank()));
    InitializeSettings("rank" + QString::number(rank()));
    InitializeModel();
    InitializeSimulator();
    InitializeObjectiveFunction();
    worker_ = new MPI::Worker(this);
    FinalizeInitialization(false);
  }

  bl_ps_fdiff_ = settings_->optimizer()->parameters().bl_ps_fdiff;
  im_ = "f diff in pattern threshold [fraction]: " + num2str(bl_ps_fdiff_);
  ext_info(im_, md_, cl_);
}

// --
void BilevelSynchrMPIRunner::initialDistribution() {
  while (optimizer_->nr_queued_cases() > 0
    && overseer_->NumberOfFreeWorkers() > 1) { // Leave one free worker
    overseer_->AssignCase(optimizer_->GetCaseForEvaluation());
  }
}

void BilevelSynchrMPIRunner::Execute() {

  // --
  auto handle_new_case = [&]() mutable {
    Optimization::Case *c_upper;
    printMessage("Getting new case from optimizer.", 2);
    c_upper = optimizer_->GetCaseForEvaluation();

    if (bookkeeper_->IsEvaluated(c_upper, true)) {
      printMessage("Case found in bookkeeper");
      c_upper->state.eval = ES::E_BOOKKEEPED;
      optimizer_->SubmitEvaluatedCase(c_upper);
    } else {
      overseer_->AssignCase(c_upper);
      printMessage("New case assigned to worker.", 2);
    }
  };

  // --
  auto wait_for_evaluated_case = [&]() mutable {
    printMessage("Waiting to receive evaluated case...", 2);

    // TODO: This is a duplicate case that wont get deleted, i.e. a MEMORY LEAK.
    Optimization::Case *evaluated_case = overseer_->RecvEvaluatedCase();

    printMessage("Evaluated case received.", 2);
    if (overseer_->last_case_tag == MPIRunner::MsgTag::CASE_EVAL_SUCCESS) {
      printMessage("Setting state for evaluated case.", 2);
      evaluated_case->state.eval = ES::E_DONE;
      printMessage("Setting timings for evaluated case.", 2);
      if (!is_ensemble_run_ && optimizer_->GetSimulationDuration(evaluated_case) > 0) {
        printMessage("Setting timings for evaluated case.", 2);
        sim_times_.push_back(optimizer_->GetSimulationDuration(evaluated_case));
      }
    }

    optimizer_->SubmitEvaluatedCase(evaluated_case);
    printMessage("Submitted evaluated case to optimizer.", 2);
  };

  if (rank() == 0) { // Overseer
    printMessage("Performing initial distribution...", 2);
    initialDistribution();
    printMessage("Initial distribution done.", 2);

    while (optimizer_->IsFinished() == false) {

      // if (is_ensemble_run_) {
      //   printMessage(ensemble_helper_.GetStateString(), 2);
      // }

      // if (is_ensemble_run_ && ensemble_helper_.IsCaseAvailableForEval()) {
      //   printMessage("Queued realization cases available.", 2);
      //   if (overseer_->NumberOfFreeWorkers() > 0) { // Free workers available
      //     printMessage("Free workers available. Handling next case.", 2);
      //     handle_new_case();
      //   } else { // No workers available
      //     printMessage("No free workers available. Waiting for an evaluated case.", 2);
      //     wait_for_evaluated_case();
      //   }
      //
      // } else if (is_ensemble_run_ && !ensemble_helper_.IsCaseDone()) {
      //   printMessage("Not all ensemble realizations have been evaluated.", 2);
      //   wait_for_evaluated_case();
      //
      // } else
      //
      if (optimizer_->nr_queued_cases() > 0) { // Queued cases in optimizer
        printMessage("Queued cases available.", 2);
        if (overseer_->NumberOfFreeWorkers() > 0) { // Free workers available
          printMessage("Free workers available. Handling next case.", 2);
          handle_new_case();
        } else { // No workers available
          printMessage("No free workers available. Waiting for an evaluated case.", 2);
          wait_for_evaluated_case();
        }
      } else { // No queued cases in optimizer
        printMessage("No queued cases available.", 2);
        if (overseer_->NumberOfBusyWorkers() == 0) { // All workers are free
          printMessage("No workers are busy. Starting next iteration.", 2);
          handle_new_case();
        } else { // Some workers are performing simulations
          printMessage("Some workers are still evaluating cases from this iteration. Waiting for evaluated cases.", 2);
          wait_for_evaluated_case();
        }
      }

    } // optimizer overseer finished

    FinalizeRun(true);
    overseer_->TerminateWorkers();
    printMessage("Terminating workers.", 2);
    overseer_->EnsureWorkerTermination();
    env_.~environment();
    return;


  } else { // Worker
    printMessage("Waiting to receive initial unevaluated case...", 2);
    worker_->RecvUnevaluatedCase();
    printMessage("Received initial unevaluated case.", 2);

    while (worker_->GetCurrentCase() != nullptr) {
      MPIRunner::MsgTag tag = MPIRunner::MsgTag::CASE_EVAL_SUCCESS; // Tag to be sent along with the case.
      try {
        model_update_done_ = false;
        simulation_done_ = false;
        logger_->AddEntry(this);

        QDateTime upper_sim_start, upper_sim_end;
        bool upper_sim_success = true;

        // if (is_ensemble_run_) {
        //   printMessage("Updating grid path.", 2);
        //   model_->set_grid_path(ensemble_helper_.GetRlz(worker_->GetCurrentCase()->GetEnsembleRlz().toStdString()).grid());
        // }

        // compute fval for testcase0
        printMessage("Applying case to model.", 2);
        model_->ApplyCase(worker_->GetCurrentCase());
        model_update_done_ = true;
        logger_->AddEntry(this);


        upper_sim_start = QDateTime::currentDateTime();

        if (rts_->sim_timeout() == 0 && settings_->simulator()->max_minutes() < 0) {
          printMessage("Starting model evaluation.", 2);
          simulator_->Evaluate();
        } else if (sim_times_.empty() && settings_->simulator()->max_minutes() > 0) {
          // if (!is_ensemble_run_) {
          printMessage("Starting model evaluation with timeout.", 2);
          upper_sim_success = simulator_->Evaluate(settings_->simulator()->max_minutes() * 60,
                                                   rts_->threads_per_sim());
          // } else {
          //   printMessage("Starting ensemble model evaluation with timeout.", 2);
          //   simulation_success = simulator_->Evaluate(ensemble_helper_.GetRlz(worker_->GetCurrentCase()->GetEnsembleRlz().toStdString()),
          //                                             settings_->simulator()->max_minutes() * 60,
          //                                             rts_->threads_per_sim());
          // }
        } else {
          // if (!is_ensemble_run_) {
          printMessage("Starting model evaluation with timeout.", 2);
          upper_sim_success = simulator_->Evaluate(timeoutVal(), rts_->threads_per_sim());
          // } else {
          //   printMessage("Starting ensemble model evaluation with timeout.", 2);
          //   simulation_success = simulator_->Evaluate(ensemble_helper_.GetRlz(worker_->GetCurrentCase()->GetEnsembleRlz().toStdString()),
          //                                             settings_->simulator()->max_minutes() * 60,
          //                                             rts_->threads_per_sim());
          // }
        }

        simulation_done_ = true;
        logger_->AddEntry(this);
        upper_sim_end = QDateTime::currentDateTime();
        int upper_sim_time = time_span_seconds(upper_sim_start, upper_sim_end);

        // prntDbg(1,optimizer_, worker_->GetCurrentCase());


        if (upper_sim_success) {
          tag = MPIRunner::MsgTag::CASE_EVAL_SUCCESS;
          printMessage("Setting objective function value.", 2);
          model_->wellCost(settings_->optimizer());

          if (settings_->optimizer()->objective().type == OT::Augmented) {
            worker_->GetCurrentCase()->set_objf_value(objf_->value(false));
          } else {
            worker_->GetCurrentCase()->set_objf_value(objf_->value());
          }

          string tm = "Objective function value set to ";
          tm += num2str(worker_->GetCurrentCase()->objf_value(), 8, 1);
          if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }

          worker_->GetCurrentCase()->state.eval = ES::E_DONE;
          worker_->GetCurrentCase()->SetSimTime(upper_sim_time);
          sim_times_.push_back(upper_sim_time);

          // start of optmzr-lower-level scope ---------------
          // if (worker_->GetCurrentCase()->fdiff > bl_ps_fdiff_) {
          if (true) {
            im_ = "worker_->GetCurrentCase()->state.fdiff =>";
            im_ += num2str(worker_->GetCurrentCase()->getFDiff(), 8, 1);
            im_ += " >? " + num2str(bl_ps_fdiff_, 3);
            if (vp_.vRUN >= 1) { info(im_, vp_.lnw); }

            std::unique_ptr<Optimization::Optimizers::DFTR>
              optmzr_lwr_(new Optimization::Optimizers::DFTR(
              settings_->optimizer(), worker_->GetCurrentCase(),
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

                // prntDbg(2, optmzr_lwr_.get(), c_lower);

                lower_sim_start = QDateTime::currentDateTime();
                if (rts_->sim_timeout() == 0 && settings_->simulator()->max_minutes() < 0) {
                  lower_sim_success = true;
                  simulator_->Evaluate();
                } else {
                  lower_sim_success = simulator_->Evaluate(timeoutVal(), rts_->threads_per_sim());
                }
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

            optmzr_lwr_->GetTentativeBestCase()->setNrSims(optmzr_lwr_->get_c_evald());
            worker_->GetCurrentCase()->CopyCaseVals(optmzr_lwr_->GetTentativeBestCase());
          } // -------------------------------------------------


        } else {
          cout << "WARNING: CASE_EVAL_TIMEOUT!" << endl;
          tag = MPIRunner::MsgTag::CASE_EVAL_TIMEOUT;
          printMessage("Timed out. Setting objective function value to SENTINEL VALUE.", 2);
          worker_->GetCurrentCase()->state.eval = ES::E_TIMEOUT;
          worker_->GetCurrentCase()->state.err_msg = EM::ERR_SIM;
          worker_->GetCurrentCase()->set_objf_value(sentinelValue());
        }

      } catch (std::runtime_error e) {
        std::cout << e.what() << std::endl;
        tag = MPIRunner::MsgTag::CASE_EVAL_INVALID;
        worker_->GetCurrentCase()->state.eval = ES::E_FAILED;
        worker_->GetCurrentCase()->state.err_msg = EM::ERR_WIC;
        printMessage("Invalid case. Setting objective function value to SENTINEL VALUE.", 2);
        worker_->GetCurrentCase()->set_objf_value(sentinelValue());
      }
      printMessage("Sending back evaluated case.", 2);
      worker_->SendEvaluatedCase(tag);

      printMessage("Waiting to receive an unevaluated case...", 2);
      worker_->RecvUnevaluatedCase();

      if (worker_->GetCurrentTag() == TERMINATE) {
        printMessage("Received termination message. Breaking.", 2);
        break;
      } else {
        printMessage("Received an unevaluated case.", 2);
      }
    }
    FinalizeRun(false);
    printMessage("Finalized on worker.", 2);
    worker_->ConfirmFinalization();
    env_.~environment();
    return;
  }




}







Loggable::LogTarget BilevelSynchrMPIRunner::GetLogTarget() {
  return STATE_RUNNER;
}

map<string, string> BilevelSynchrMPIRunner::GetState() {
  map<string, string> statemap;
  statemap["case-desc"] = worker_->GetCurrentCase()->StringRepresentation(model_->variables());
  statemap["mod-update-done"] = model_update_done_ ? "yes" : "no";
  statemap["sim-done"] = simulation_done_ ? "yes" : "no";
  statemap["last-update"] = timestamp_string();
  return statemap;
}

QUuid BilevelSynchrMPIRunner::GetId() {
  return QUuid();
}






}
}