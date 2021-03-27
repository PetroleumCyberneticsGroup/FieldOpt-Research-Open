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

ICVConstraint::ICVConstraint(SO& seto, VPC *vars, SV vp)
  : Constraint(seto, vars, vp) {

  assert(!seto.wells.empty());
  assert(seto.min < seto.max);

  string ws;
  if (seto.well.isEmpty()) {
    for (const auto& wn : seto.wells) { ws += wn.toStdString() + "; "; }
  } else { ws = seto.well.toStdString(); }

  if (vp_.vOPT >= 3) {
    im_ = "Adding ICV constraint for " + ws;
    ext_info(im_, md_, cl_, vp_.lnw);
  }

  min_ = seto.min;
  max_ = seto.max;

  icd_cnstrnd_well_nms_ = seto.wells;
  penalty_weight_ = seto.penalty_weight;

  if (vp_.vOPT >= 3) {
    im_ = "Adding ICV constraint with [min, max] = [";
    im_ += num2str(min_, 5) + ", " + num2str(max_, 5);
    im_ += "] for well " + ws + " with variables: ";
  }

  for (auto &wname : icd_cnstrnd_well_nms_) {
    auto icd_vars = vars->GetWellICDVars(wname);
    icd_cnstrnd_real_vars_.append(icd_vars);

    for (auto var : icd_vars) {
      var->setBounds(min_, max_);
      icd_cnstrnd_uuid_vars_.push_back(var->id());
      im_ += var->name().toStdString() + ", ";
    }
  }

  if(seto.scaling_) {
    min_ = -1.0;
    max_ = 1.0;
  }

  if (vp_.vOPT >= 3) { ext_info(im_, md_, cl_, vp_.lnw); }
}

// Keep for ref
// bool ICVConstraint::CaseSatisfiesConstraint(Case *c) {
//   for (auto id : affected_variables_) {
//     if (c->real_variables()[id] > max_ || c->real_variables()[id] < min_) {
//       return false;
//     }
//   }
//   return true;
// }

bool ICVConstraint::CaseSatisfiesConstraint(Case *c) {
  for (int ii=0; ii < icd_cnstrnd_real_vars_.size(); ii++ ) {
    auto var_id = icd_cnstrnd_uuid_vars_[ii];
    double c_val = c->get_real_variable_value(var_id);
    if (c_val > max_ || c_val < min_) { return false; }
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
