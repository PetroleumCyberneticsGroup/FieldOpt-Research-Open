/***********************************************************
Created by einar on 11/3/16.
Copyright (C) 2015-2019
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 02/19/19 Einar J. M. Baumann
<einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
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

#include <iostream>
#include "GSS.h"
#include "Utilities/math.hpp"
#include "gss_patterns.hpp"
#include "Utilities/stringhelpers.hpp"
#include "Utilities/printer.hpp"

namespace Optimization {
namespace Optimizers {

using Printer::ext_info;
using Printer::ext_warn;
using Printer::num2str;
using Printer::DBG_prntVecXdSz;

GSS::
GSS(Settings::Optimizer *settings,
    Case *base_case,
    Model::Properties::VarPropContainer *variables,
    Reservoir::Grid::Grid *grid,
    Logger *logger,
    CaseHandler *case_handler,
    Constraints::ConstraintHandler *constraint_handler)
  : Optimizer(settings, base_case, variables, grid,
              logger, case_handler, constraint_handler) {
  vp_ = settings->verbParams();

  contr_fac_ = settings->parameters().contraction_factor;
  assert(contr_fac_ < 1.0);

  expan_fac_ = settings->parameters().expansion_factor;
  assert(expan_fac_ >= 1.0);

  int numRvars = base_case->GetRealVarVector().size();
  int numIvars = base_case->GetIntegerVarVector().size();
  num_vars_ = numRvars + numIvars;

  if (numRvars > 0 && numIvars > 0) {
    im_ = "CS cannot handle both continuous and discrete variables at the same time";
  } else {
    im_ = "#.Rvars=" + num2str(numRvars) + "; #.Ivars=" + num2str(numIvars);
  }
  if (vp_.vOPT >= 3) { ext_info(im_, md_, cl_); }

  directions_ = GSSPatterns::Compass(num_vars_);

  if (!settings->parameters().auto_step_lengths) {
    step_lengths_ = Eigen::VectorXd(directions_.size());
    step_lengths_.fill(settings->parameters().initial_step_length);
    step_tol_ = Eigen::VectorXd(directions_.size());
    step_tol_.fill(settings->parameters().minimum_step_length);

  } else {
    assert(constraint_handler_->HasBoundaryConstraints());
    auto lower_bound = constraint_handler_->GetLowerBounds(
      base_case->GetRealVarIdVector());
    auto upper_bound = constraint_handler_->GetUpperBounds(
      base_case->GetRealVarIdVector());

    auto difference = upper_bound - lower_bound;
    Eigen::VectorXd base = Eigen::VectorXd(directions_.size());
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < difference.size(); ++j) {
        base(i*num_vars_+j) = difference(j);
      }
    }
    step_lengths_ = base * settings->parameters().auto_step_init_scale;
    step_tol_ = base * settings->parameters().auto_step_conv_scale;
  }

  if (vp_.vOPT >= 4) {
    ext_info("Step lengths: " + DBG_prntVecXdSz(step_lengths_), md_, cl_);
    ext_info("Step tols: " + DBG_prntVecXdSz(step_tol_), md_, cl_);
  }
  assert(step_lengths_.size() == directions_.size());
  assert(step_lengths_.size() == step_tol_.size());
}

Optimizer::TerminationCondition GSS::IsFinished() {
  TerminationCondition tc = NOT_FINISHED;
  if (!case_handler_->CasesBeingEvaluated().empty()) {
    return tc;
  }

  if (evaluated_cases_ >= max_evaluations_) {
    tc = MAX_EVALS_REACHED;
  } else if (is_converged()) {
    tc = MIN_STEP_LENGTH_REACHED;
  } else if (is_sufficient_improvement()) {
        tc = SUFFICIENT_IMPROVEMENT;
  }

  if (tc != NOT_FINISHED) {
    if (enable_logging_) {
      logger_->AddEntry(this);
      logger_->AddEntry(new Summary(this, tc));
    }
  }
  return tc;
}

void GSS::expand(vector<int> dirs) {
  if (dirs[0] == -1) {
    step_lengths_ = step_lengths_ * expan_fac_;
  } else {
    for (int dir : dirs)
      step_lengths_(dir) = step_lengths_(dir) * expan_fac_;
  }
}

void GSS::contract(vector<int> dirs) {
  if (dirs[0] == -1) {
    step_lengths_ = step_lengths_ * contr_fac_;
  } else {
    for (int dir : dirs)
      step_lengths_(dir) = step_lengths_(dir) * contr_fac_;
  }
}

QList<Case *> GSS::generate_trial_points(vector<int> dirs) {
  auto trial_points = QList<Case *>();
  if (dirs[0] == -1) {
    dirs = range(0, (int)directions_.size(), 1);
  }

  VectorXi int_base = GetTentBestCase()->GetIntegerVarVector();
  VectorXd rea_base = GetTentBestCase()->GetRealVarVector();

  for (int dir : dirs) {
    auto trial_point = new Case(GetTentBestCase());

    if (int_base.size() > 0) {
      trial_point->SetIntegerVarValues(perturb(int_base, dir));
    } else if (rea_base.size() > 0) {
      trial_point->SetRealVarValues(perturb(rea_base, dir));
    }

    trial_point->set_origin_data(GetTentBestCase(),
                                 dir, step_lengths_(dir));
    trial_points.append(trial_point);
  }

  if (constraint_handler_ != nullptr) {
    for (Case *c : trial_points) {
      constraint_handler_->SnapCaseToConstraints(c);
    }
  } else if (vp_.vOPT > 3) {
    wm_ = "No constraints (constraint_handler_ == nullptr).";
    ext_warn(wm_, md_, cl_, vp_.lnw);
  }

  return trial_points;
}

template<typename T>
Matrix<T, Dynamic, 1> GSS::perturb(Matrix<T, Dynamic, 1> base, int dir) {
  Matrix<T, Dynamic, 1> dirc = directions_[dir].cast<T>();
  T sl = step_lengths_(dir);
  Matrix<T, Dynamic, 1> perturbation = base + dirc * sl;
  return perturbation;
}

bool GSS::is_converged() {
  for (int i = 0; i < step_lengths_.size(); ++i) {
    if (step_lengths_(i) >= step_tol_(i)) {
      return false;
    }
  }
  return true;
}

bool GSS::is_sufficient_improvement() {
  if ((iteration_ != 0) &&
      (tentative_best_case_->objf_value() != base_case_obj_) &&
      (sufficient_improvement_tol_ >=0) &&
      (abs((tentative_best_case_->objf_value()-base_case_obj_)/base_case_obj_) >= sufficient_improvement_tol_))
    return true;
  else
    return false;
}

void GSS::set_step_lengths(int dir_idx, double len) {

  if (constraint_handler_ != nullptr) { // All actual cases

    if (!constraint_handler_->HasBoundaryConstraints()) {
      if (vp_.vOPT >= 2) {
        im_ = "No boundary constraints. Setting all step ";
        im_ += "lengths to " + num2str(len);
        ext_warn(im_, md_, cl_, vp_.lnw);
      }
      for (int i = 0; i < step_lengths_.size(); ++i) {
        step_lengths_[i] = max(step_tol_[i], len);
      }
      step_lengths_.fill(len);
      return;
    }

  } else { // constraint_handler_ == nullptr in unit tests

    if (vp_.vOPT > 3) {
      wm_ = "constraint_handler_ == nullptr.";
      ext_warn(wm_, md_, cl_, vp_.lnw);
    }
    for (int i = 0; i < step_lengths_.size(); ++i) {
      step_lengths_[i] = max(step_tol_[i], len);
    }
    step_lengths_.fill(len);
    return;
  }

  // HasBoundaryConstraints == true
  Eigen::VectorXd lbs, ubs;
  lbs = constraint_handler_->GetLowerBounds(
    tentative_best_case_->GetRealVarIdVector());

  ubs = constraint_handler_->GetUpperBounds(
    tentative_best_case_->GetRealVarIdVector());

  int dir = dir_idx > num_vars_ - 1 ? dir_idx - num_vars_ : dir_idx;
  double S = (len - lbs[dir]) / (ubs[dir] - lbs[dir]);
  for (int i = 0; i < num_vars_; ++i) {
    step_lengths_[i] = max(step_tol_[i], lbs[i] + S * (ubs[i] - lbs[i]));
    step_lengths_[i+num_vars_] = max(step_tol_[i], lbs[i] + S * (ubs[i] - lbs[i]));
  }
}

Case * GSS::dequeue_case_with_worst_origin() {
  auto queued_cases = case_handler_->QueuedCases();
  std::sort(queued_cases.begin(), queued_cases.end(),
            [this](const Case *c1, const Case *c2) -> bool {
              return isBetter(c1, c2);
            });
  case_handler_->DequeueCase(queued_cases.last()->id());
  return queued_cases.last();
}

}
}
