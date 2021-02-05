/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>
Modified 2017/08/22 by Einar Einar J.M. Baumann

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
#include "bhp_constraint.h"

namespace Optimization {
namespace Constraints {

BhpConstraint::BhpConstraint(Settings::Optimizer::Constraint const settings,
                             Model::Properties::VarPropContainer *variables,
                             Settings::VerbParams vp)
                             : Constraint(vp) {

  assert(settings.wells.size() > 0);
  assert(settings.min < settings.max);

  if (vp_.vOPT >= 1) {
    info("Adding BHP constraint for " + settings.well.toStdString());
  }

  min_ = settings.min;
  max_ = settings.max;

  bhp_cnstrnd_well_nms_ = settings.wells;
  penalty_weight_ = settings.penalty_weight;

  if (vp_.vOPT >= 1) {
    im_ = "Adding BHP constraint with [min, max] = [";
    im_ += num2str(min_, 5) + ", " + num2str(max_, 5);
    im_ += "] for well " + settings.well.toStdString() + " with variables: ";
  }

  for (auto &wname : bhp_cnstrnd_well_nms_) {
    auto bhp_vars = variables->GetWellBHPVars(wname);
    for (auto var : bhp_vars) {
      var->setBounds(min_, max_);
      bhp_cnstrnd_uuid_vars_.push_back(var->id());
      im_ += var->name().toStdString() + ", ";
    }
    bhp_cnstrnd_real_vars_.append(bhp_vars);
  }

  if(settings.scaling_) {
    min_ = -1.0;
    max_ = 1.0;
  }

  ext_info(im_, md_, cl_, vp_.lnw);
}

bool BhpConstraint::CaseSatisfiesConstraint(Case *c) {
  for (auto var : bhp_cnstrnd_real_vars_) {
    // double case_value = c->real_variables()[var->id()];
    // double case_value = c->get_real_variable_value(var->id());
    if (var->value() > max_ || var->value() < min_)
      return false;
  }
  return true;
}

bool BhpConstraint::IsBoundConstraint() const {
  return true;
}

Eigen::VectorXd BhpConstraint::GetLowerBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd lbounds(id_vector.size());
  lbounds.fill(0);
  for (auto var : bhp_cnstrnd_real_vars_) {
    lbounds[id_vector.indexOf(var->id())] = min_;
  }
  return lbounds;
}

Eigen::VectorXd BhpConstraint::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);
  for (auto var : bhp_cnstrnd_real_vars_) {
    ubounds[id_vector.indexOf(var->id())] = max_;
  }
  return ubounds;
}

void BhpConstraint::SnapCaseToConstraints(Case *c) {
  for (auto var : bhp_cnstrnd_real_vars_) {
    if (c->get_real_variable_value(var->id()) > max_) {
      c->set_real_variable_value(var->id(), max_);
    } else if (c->get_real_variable_value(var->id()) < min_) {
      c->set_real_variable_value(var->id(), min_);
    }
  }
}

}
}
