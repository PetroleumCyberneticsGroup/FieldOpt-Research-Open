/***********************************************************
Copyright (C) 2015-2016
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#include "rate_constraint.h"

namespace Optimization {
namespace Constraints {

RateConstraint::RateConstraint(Settings::Optimizer::Constraint const settings,
                               Model::Properties::VarPropContainer *variables,
                               Settings::VerbParams vp)
                               : Constraint(vp) {

  assert(!settings.wells.empty());
  assert(settings.min < settings.max);

  if (vp_.vOPT >= 1) {
    info("Adding Rate constraint for " + settings.well.toStdString());
  }

  affected_well_names_ = settings.wells;
  min_ = settings.min;
  max_ = settings.max;
  penalty_weight_ = settings.penalty_weight;

  for (auto &wname : affected_well_names_) {
    affected_real_variables_.append(variables->GetWellRateVars(wname));
  }
}

bool RateConstraint::CaseSatisfiesConstraint(Case *c) {
  for (auto var : affected_real_variables_) {
    double case_value = c->get_real_variable_value(var->id());
    if (case_value > max_ || case_value < min_)
      return false;
  }
  return true;
}

void RateConstraint::SnapCaseToConstraints(Case *c) {
  for (auto var : affected_real_variables_) {
    if (c->get_real_variable_value(var->id()) > max_)
      c->set_real_variable_value(var->id(), max_);
    else if (c->get_real_variable_value(var->id()) < min_)
      c->set_real_variable_value(var->id(), min_);
  }
}

Eigen::VectorXd RateConstraint::GetLowerBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd lbounds(id_vector.size());
  lbounds.fill(0);
  for (auto var : affected_real_variables_) {
    lbounds[id_vector.indexOf(var->id())] = min_;
  }
  return lbounds;
}

Eigen::VectorXd RateConstraint::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);
  for (auto var : affected_real_variables_) {
    ubounds[id_vector.indexOf(var->id())] = max_;
  }
  return ubounds;
}

bool RateConstraint::IsBoundConstraint() const {
  return true;
}

}
}
