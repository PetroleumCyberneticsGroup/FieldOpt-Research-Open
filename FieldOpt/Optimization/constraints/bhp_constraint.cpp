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

BhpConstraint::BhpConstraint(SO& seto, VPC *vars, SV vp)
  : Constraint(seto, vars, vp) {

  min_ = seto_.min;
  max_ = seto_.max;
  assert(min_ < max_);
  assert(!seto_.wells.empty());

  bhp_cnstrnd_well_nms_ = seto_.wells;
  penalty_weight_ = seto_.penalty_weight;

  PrntWellInfo("BHP", 1);

  for (auto &wname : bhp_cnstrnd_well_nms_) {
    auto bhp_vars = vars->GetWellBHPVars(wname);
    bhp_cnstrnd_real_vars_.append(bhp_vars);

    for (auto var : bhp_vars) {
      var->setBounds(min_, max_);
      bhp_cnstrnd_uuid_vars_.push_back(var->id());
      im_ += var->name().toStdString() + ", ";
    }
  }

  if(!bhp_cnstrnd_real_vars_.empty()) { isEnabled_ = true; }

  if(seto.scaling_) {
    min_ = -1.0;
    max_ = 1.0;
  }

  if (vp_.vOPT >= 3) { ext_info(im_, md_, cl_, vp_.lnw); }
}

// Keep for ref
// bool BhpConstraint::CaseSatisfiesConstraint(Case *c) {
//   for (auto var : bhp_cnstrnd_real_vars_) {
//     double case_value = c->get_real_variable_value(var->id());
//     if (case_value > max_ || case_value < min_) {
//       return false;
//     }
//   }
//   return true;
// }

bool BhpConstraint::CaseSatisfiesConstraint(Case *c) {
  for (int ii=0; ii < bhp_cnstrnd_real_vars_.size(); ii++ ) {
    auto var_id = bhp_cnstrnd_uuid_vars_[ii];
    double c_val = c->get_real_variable_value(var_id);
    if (c_val > max_ || c_val < min_) { return false; }
  }
  return true;
}

void BhpConstraint::SnapCaseToConstraints(Case *c) {
  string tm;
  if (vp_.vOPT >= 4) {
    tm = "Check bounds: [" + DBG_prntDbl(min_) + DBG_prntDbl(max_) + "] ";
    tm += "for case: " + c->id_stdstr();
    pad_text(tm, vp_.lnw );
    tm += "x: " + DBG_prntVecXd(c->GetRealVarVector());
    ext_info(tm, md_, cl_, vp_.lnw);
  }

  tm = "";
  for (auto id : bhp_cnstrnd_uuid_vars_) {
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



}
}
