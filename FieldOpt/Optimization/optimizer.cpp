/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb.pcg@gmail.com>

Modified 2019-2020 Thiago Lima Silva
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

#include <Utilities/time.hpp>
#include "optimizer.h"
#include <time.h>
#include <cmath>

namespace Optimization {

using Printer::ext_info;
using Printer::num2str;

using Constraints::ConstraintHandler;

Optimizer::Optimizer(Settings::Optimizer *opt_settings,
                     Case *base_case,
                     Model::Properties::VarPropContainer *variables,
                     Reservoir::Grid::Grid *grid,
                     Logger *logger,
                     CaseHandler *case_handler,
                     Constraints::ConstraintHandler *constraint_handler) {

  // Verify that the base case has been evaluated.
  try {
    base_case->objf_value();
  } catch (ObjectiveFunctionException) {
    throw OptimizerInitializationException(
      "The objective function value of the base case "
      "must be set before initializing an Optimizer.");
  }

  max_evaluations_ = opt_settings->parameters().max_evaluations;
  tentative_best_case_iteration_ = 0;
  seconds_spent_in_iterate_ = 0;

  // CONSTRAINT HANDLER ------------------------------------
  // Moved constraint handliner initialization to Model
  // if (constraint_handler == nullptr) {
  //   constraint_handler_ = new ConstraintHandler(opt_settings,
  //                                               variables, grid);
  // } else {
  //   ext_info("Using shared ConstraintHandler.", md_, cl_);
  // }
  if (constraint_handler != nullptr) {
    constraint_handler_ = constraint_handler;
  } else {
    im_ = "constraint_handler == nullptr in constructor call ";
    im_ += "(ok for tests).";
    ext_warn(im_, md_, cl_, vp_.lnw);
  }

  // BASE CASE -> NEW CASE ---------------------------------
  tentative_best_case_ = base_case;

  if (case_handler == nullptr) {
    case_handler_ = new CaseHandler(tentative_best_case_);
  } else {
    ext_info("Using shared CaseHandler.", md_, cl_);
    case_handler_ = case_handler;
    is_hybrid_component_ = true;
  }

  // ITERATION ---------------------------------------------
  iteration_ = 0;
  evaluated_cases_ = 0;
  mode_ = opt_settings->mode();
  type_ = opt_settings->type();

  is_async_ = false;
  start_time_ = QDateTime::currentDateTime();
  logger_ = logger;
  enable_logging_ = true;

  vp_ = opt_settings->verbParams();

  // PENALIZATION ------------------------------------------
  penalize_ = opt_settings->objective().use_penalty_function;

  if (penalize_) {
    if (!normalizer_ofv_.is_ready()) {
      if (vp_.vOPT >= 1) {
        ext_info("Initializing normalizers", md_, cl_);
      }
      initializeNormalizers();

      // penalize the base case
      double org_ofv = tentative_best_case_->objf_value();
      double pen_ofv = PenalizedOFV(tentative_best_case_);
      tentative_best_case_->set_objf_value(pen_ofv);
      if (vp_.vOPT >= 1) {
        im_ = "Penalized base case. Orig.val=" + num2str(org_ofv);
        im_ += "; Pen.val=" + num2str(pen_ofv);
        ext_info(im_, md_, cl_);
      }
    }
  }
}

Settings::Optimizer::OptimizerType TR_DFO =
  Settings::Optimizer::OptimizerType::TrustRegionOptimization;

Settings::Optimizer::OptimizerType DFTR =
  Settings::Optimizer::OptimizerType::DFTR;

Case *Optimizer::GetCaseForEvaluation() {
  if (case_handler_->QueuedCases().empty()) {
    time_t start, end;
    time(&start);
    iterate();
    time(&end);
    seconds_spent_in_iterate_ = (int)difftime(end, start);
  }

  if (IsFinished() || (case_handler_->QueuedCases().empty())) {
    im_ = "(IsFinished() == NOT_FINISHED) => ";
    im_ += (IsFinished() == NOT_FINISHED) ? "true" : "false";
    im_ += " (case_handler_->QueuedCases()..empty() == 0) => ";
    im_ += case_handler_->QueuedCases().empty() ? "true" : "false";
    im_ += " Returning nullptr";

    if (type_ == TR_DFO || type_ == DFTR) {
      if (vp_.vOPT >= 3) { ext_info(im_, md_, cl_); }
      return nullptr;
    }
  }
  return case_handler_->GetNextCaseForEvaluation();
}

void Optimizer::SubmitEvaluatedCase(Case *c) {
  if (penalize_ && iteration_ > 0) {
    double penalized_ofv = PenalizedOFV(c);
    c->set_objf_value(penalized_ofv);
  }
  case_handler_->UpdateCaseObjectiveFunctionValue(c->id(), c->objf_value());
  case_handler_->SetCaseState(c->id(), c->state, c->GetWICTime(), c->GetSimTime());
  case_handler_->SetCaseEvaluated(c->id());
  handleEvaluatedCase(case_handler_->GetCase(c->id()));

  if(c->state.eval == Case::CaseState::E_DONE) {
    evaluated_cases_++;
  }

  if (enable_logging_) {
    logger_->AddEntry(case_handler_->GetCase(c->id()));
  }
}

Case *Optimizer::GetTentativeBestCase() const {
  return tentative_best_case_;
}

bool Optimizer::isImprovement(const Case *c) {
  return isBetter(c, tentative_best_case_);
}

bool Optimizer::isBetter(const Case *c1, const Case *c2) const {
  if (mode_ == Settings::Optimizer::OptimizerMode::Maximize) {
    if (c1->objf_value() > c2->objf_value())
      return true;
  }
  else if (mode_ == Settings::Optimizer::OptimizerMode::Minimize) {
    if (c1->objf_value() < c2->objf_value())
      return true;
  }
  return false;
}

QString Optimizer::GetStatusStringHeader() const {
  return QString("%1,%2,%3,%4,%5,%6\n")
    .arg("Iteration")
    .arg("EvaluatedCases")
    .arg("QueuedCases")
    .arg("RecentlyEvaluatedCases")
    .arg("TentativeBestCaseID")
    .arg("TentativeBestCaseOFValue");
}

QString Optimizer::GetStatusString() const {
  return QString("%1,%2,%3,%4,%5,%6\n")
    .arg(iteration_)
    .arg(nr_evaluated_cases())
    .arg(nr_queued_cases())
    .arg(nr_recently_evaluated_cases())
    .arg(tentative_best_case_->id().toString())
    .arg(tentative_best_case_->objf_value());
}

void Optimizer::EnableConstraintLogging(const QString& output_dir_path) {
  for (Constraints::Constraint *con : constraint_handler_->constraints()) {
    con->EnableLogging(output_dir_path);
  }
}

int Optimizer::GetSimulationDuration(Case *c) {
  auto cs = case_handler_->GetCase(c->id());
  if (cs->state.eval != Case::CaseState::EvalStatus::E_DONE) {
    return -1;
  }
  return c->GetSimTime();
}

Loggable::LogTarget Optimizer::GetLogTarget() {
  return Loggable::LogTarget::LOG_OPTIMIZER;
}

map<string, string> Optimizer::GetState() {
  return map<string, string>();
}

QUuid Optimizer::GetId() {
  return tentative_best_case_->GetId();
}

map<string, vector<double>> Optimizer::GetValues() {
  map<string, vector<double>> valmap;
  valmap["TimeEl"] = vector<double>{time_since_seconds(start_time_)};
  valmap["IterNr"] = vector<double>{iteration_};
  valmap["TimeIt"] = vector<double>{seconds_spent_in_iterate_};
  valmap["TotlNr"] = vector<double>{case_handler_->NumberTotal()};
  valmap["EvalNr"] = vector<double>{case_handler_->NumberSimulated()};
  valmap["BkpdNr"] = vector<double>{case_handler_->NumberBookkeeped()};
  valmap["TimONr"] = vector<double>{case_handler_->NumberTimeout()};
  valmap["FailNr"] = vector<double>{case_handler_->NumberFailed()};
  valmap["InvlNr"] = vector<double>{case_handler_->NumberInvalid()};
  valmap["CBOFnV"] = vector<double>{tentative_best_case_->objf_value()};
  return valmap;
}

Loggable::LogTarget Optimizer::Summary::GetLogTarget() {
  return LOG_SUMMARY;
}

map<string, string> Optimizer::Summary::GetState() {

  map<string, string> statemap  = ext_state_;
  statemap["Start"] = timestamp_string(opt_->start_time_);
  statemap["Duration"] = timespan_string(
    time_span_seconds(opt_->start_time_, QDateTime::currentDateTime())
                                        );
  statemap["End"] = timestamp_string(QDateTime::currentDateTime());

  switch (cond_) {
    case MAX_EVALS_REACHED: {
      statemap["Term. condition"] = "Reached max. sims";
      break;
    }
    case MIN_STEP_LENGTH_REACHED: {
      statemap["Term. condition"] = "Reached min. step length";
      break;
    }
    case MAX_ITERS_REACHED: {
      statemap["Term. condition"] = "Reached max. iterations";
      break;
    }
    default: {
      statemap["Term. condition"] = "Unknown";
    }
  }

  statemap["bc Best case found in iter"] = boost::lexical_cast<string>(opt_->tentative_best_case_iteration_);
  statemap["bc UUID"] = opt_->tentative_best_case_->GetId().toString().toStdString();
  statemap["bc Objective function value"] = boost::lexical_cast<string>(opt_->tentative_best_case_->objf_value());
  statemap["bc Constraint status"] = statemap["bc Constraint status"] = opt_->tentative_best_case_->GetState()["ConsSt"];
  statemap["bc Simulation time"] = timespan_string(opt_->tentative_best_case_->GetSimTime());
  return statemap;
}

QUuid Optimizer::Summary::GetId() {
  return opt_->tentative_best_case_->GetId();
}

map<string, vector<double>> Optimizer::Summary::GetValues() {
  map<string, vector<double>> valmap;
  valmap["generated"] = vector<double>{opt_->case_handler_->NumberTotal()};
  valmap["simulated"] = vector<double>{opt_->case_handler_->NumberSimulated()};
  valmap["invalid"] = vector<double>{opt_->case_handler_->NumberInvalid()};
  valmap["failed"] = vector<double>{opt_->case_handler_->NumberFailed()};
  valmap["timed out"] = vector<double>{opt_->case_handler_->NumberTimeout()};
  valmap["bookkeeped"] = vector<double>{opt_->case_handler_->NumberBookkeeped()};
  return valmap;
}

void Optimizer::updateTentativeBestCase(Case *c) {
  tentative_best_case_ = c;
  tentative_best_case_iteration_ = iteration_;
}

void Optimizer::initializeNormalizers() {
  initializeOfvNormalizer();
  if (constraint_handler_ != nullptr) { // All actual cases, i.e., not unit tests
    constraint_handler_->InitializeNormalizers(case_handler_->AllCases());
  }
}

void Optimizer::initializeOfvNormalizer() {
  if (case_handler_->EvaluatedCases().empty()
    || normalizer_ofv_.is_ready()) {
    em_ = "Unable to initialize normalizer with no evaluated cases available.";
    throw runtime_error(em_);
  }

  vector<double> abs_ofvs;
  for (auto c : case_handler_->EvaluatedCases()) {
    abs_ofvs.push_back(abs(c->objf_value()));
  }
  long double max_ofv = *max_element(abs_ofvs.begin(), abs_ofvs.end());

  normalizer_ofv_.set_max(1.0L);
  normalizer_ofv_.set_midpoint(max_ofv);
  normalizer_ofv_.set_steepness(1.0L / max_ofv);
}

double Optimizer::PenalizedOFV(Case *c) {
  long double norm_ofv = normalizer_ofv_.normalize(c->objf_value());
  long double penalty = constraint_handler_->GetWeightedNormalizedPenalties(c);
  long double norm_pen_ovf = norm_ofv - penalty;
  double denormalized_ofv = normalizer_ofv_.denormalize(norm_pen_ovf);

  if (vp_.vOPT >= 3) {
    im_ = "Penalized case " + c->id().toString().toStdString()  + ". ";
    im_ += "Initial OFV: " + num2str(c->objf_value()) + "; ";
    im_ += "Normalized OFV :" + num2str(norm_ofv) + "; ";
    ext_info(im_, md_, cl_);
  }

  if (norm_pen_ovf <= 0.0L) {
    cout << "RETURNING ZERO OFV" << endl;
    return 0.0;
  } else {
    return denormalized_ofv;
  }
}

void Optimizer::DisableLogging() {
  enable_logging_ = false;
}
}

