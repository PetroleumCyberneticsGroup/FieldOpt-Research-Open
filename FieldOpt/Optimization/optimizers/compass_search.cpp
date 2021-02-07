/***********************************************************
Copyright (C) 2015-2017
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

#include <iostream>
#include "compass_search.h"
#include "gss_patterns.hpp"
#include "Utilities/verbosity.h"
#include "Utilities/printer.hpp"
#include "Utilities/stringhelpers.hpp"

namespace Optimization {
namespace Optimizers {

CompassSearch::CompassSearch(Settings::Optimizer *settings,
                             Case *base_case,
                             Model::Properties::VarPropContainer *variables,
                             Reservoir::Grid::Grid *grid,
                             Logger *logger,
                             CaseHandler *case_handler,
                             Constraints::ConstraintHandler *constraint_handler)
  : GSS(settings, base_case, variables, grid, logger, case_handler, constraint_handler) {
  vp_ = settings->verbParams();
  if (enable_logging_) {
    logger_->AddEntry(this);
  }
}

void CompassSearch::iterate() {
  if (vp_.vOPT >= 2) {
    im_ = "Compass search iteration " + num2str(iteration_) + ".| ";
    im_ += "Step lengths: " + eigenvec_to_str(step_lengths_);
    ext_info(im_);
  }

  if (enable_logging_) {
    logger_->AddEntry(this);
  }

  if (!is_successful_iteration() && iteration_ != 0) {
    contract();
  }

  case_handler_->AddNewCases(generate_trial_points());
  case_handler_->ClearRecentlyEvaluatedCases();
  iteration_++;
}

QString CompassSearch::GetStatusStringHeader() const {
  return QString("%1,%2,%3,%4,%5,%6,%7")
    .arg("Iteration")
    .arg("EvaluatedCases")
    .arg("QueuedCases")
    .arg("RecentlyEvaluatedCases")
    .arg("TentativeBestCaseID")
    .arg("TentativeBestCaseOFValue")
    .arg("StepLength");
}

QString CompassSearch::GetStatusString() const {
  return QString("%1,%2,%3,%4,%5,%6,%7")
    .arg(iteration_)
    .arg(nr_evaluated_cases())
    .arg(nr_queued_cases())
    .arg(nr_recently_evaluated_cases())
    .arg(GetTentativeBestCase()->id().toString())
    .arg(GetTentativeBestCase()->objf_value())
    .arg(step_lengths_[0]);
}

void CompassSearch::handleEvaluatedCase(Case *c) {
  if (isImprovement(c))
    updateTentativeBestCase(c);
}

bool CompassSearch::is_successful_iteration() {
  return case_handler_->RecentlyEvaluatedCases().contains(GetTentativeBestCase());
}

}}
