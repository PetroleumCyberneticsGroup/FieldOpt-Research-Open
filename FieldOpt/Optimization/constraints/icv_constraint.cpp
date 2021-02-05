/***********************************************************
Copyright (C) 2015-2018
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

#include <Utilities/printer.hpp>
#include <Utilities/verbosity.h>
#include "icv_constraint.h"
#include "Utilities/stringhelpers.hpp"

namespace Optimization {
namespace Constraints {

ICVConstraint::ICVConstraint(const Settings::Optimizer::Constraint& settings,
                             Model::Properties::VarPropContainer *variables,
                             Settings::VerbParams vp)
    : Constraint(vp) {

  variables_ = variables;

  assert(settings.wells.size() > 0);
  assert(settings.min < settings.max);

  if (vp_.vOPT >= 1) {
    info("Adding ICV constraint for " + settings.well.toStdString());
  }

  min_ = settings.min;
  max_ = settings.max;

  icd_cnstrnd_well_nms_ = settings.wells;
  penalty_weight_ = settings.penalty_weight;

  if (vp_.vOPT >= 1) {
    im_ = "Adding ICV constraint with [min, max] = [";
    im_ += num2str(min_, 5) + ", " + num2str(max_, 5);
    im_ += "] for well " + settings.well.toStdString() + " with variables: ";
  }

  for (auto &wname : icd_cnstrnd_well_nms_) {
    auto icd_vars = variables->GetWellICDVars(wname);
    for (auto var : icd_vars) {
      var->setBounds(min_, max_);
      icd_cnstrnd_uuid_vars_.push_back(var->id());
      im_ += var->name().toStdString() + ", ";
    }
    icd_cnstrnd_real_vars_.append(icd_vars);
  }

  if(settings.scaling_) {
    min_ = -1.0;
    max_ = 1.0;
  }

  ext_info(im_, md_, cl_, vp_.lnw);
}

bool ICVConstraint::CaseSatisfiesConstraint(Optimization::Case *c) {
  for (auto var : icd_cnstrnd_real_vars_) {
    if (var->value() > max_ || var->value() < min_)
      return false;
  }
  return true;
}

void ICVConstraint::SnapCaseToConstraints(Optimization::Case *c) {
  string tm;
  if (vp_.vOPT >= 4) {
    tm = "Check bounds: [" + DBG_prntDbl(min_) + DBG_prntDbl(max_) + "] ";
    tm += "for case: " + c->id_stdstr();
    pad_text(tm, vp_.lnw );
    tm += "x: " + DBG_prntVecXd(c->GetRealVarVector());
    ext_info(tm, md_, cl_, vp_.lnw);
  }

  tm = "";
  for (auto id : icd_cnstrnd_uuid_vars_) {
    if (c->get_real_variable_value(id) > max_) {
      c->set_real_variable_value(id, max_);
      tm = "Upper bound active";
    } else if (c->get_real_variable_value(id) < min_) {
      c->set_real_variable_value(id, min_);
      tm = "Lower bound active";
    }
  }

  if (vp_.vOPT >= 4 && tm.size() > 0) {
    pad_text(tm, vp_.lnw );
    tm += "x: " + DBG_prntVecXd(c->GetRealVarVector());
    ext_info(tm, md_, cl_, vp_.lnw);
  }
}

bool ICVConstraint::IsBoundConstraint() const { return true; }

Eigen::VectorXd ICVConstraint::GetLowerBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd lbounds(id_vector.size());
  lbounds.fill(0);
  for (auto id : icd_cnstrnd_uuid_vars_) {
    lbounds[id_vector.indexOf(id)] = min_;
  }
  return lbounds;
}

Eigen::VectorXd ICVConstraint::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);
  for (auto id : icd_cnstrnd_uuid_vars_) {
    ubounds[id_vector.indexOf(id)] = max_;
  }
  return ubounds;
}

}
}
