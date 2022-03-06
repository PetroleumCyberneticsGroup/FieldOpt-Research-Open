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

#include "packer_constraint.h"
#include <algorithm>
#include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>

namespace Optimization {
namespace Constraints {

using namespace Model::Properties;
using Printer::info;
using Printer::ext_info;

PackerConstraint::PackerConstraint(SO& seto, VPC *vars, SV vp)
  : Constraint(seto, vars, vp) {

  if (vp_.vOPT >= 1) {
    info("Adding Packer constraint for " + seto.well.toStdString());
  }

  for (auto var : *vars->GetContinuousVariables()) {
    auto lvar = var.second;
    if (lvar->propertyInfo().prop_type == Property::PropertyType::Packer
      && QString::compare(lvar->propertyInfo().parent_well_name, seto.well) == 0) {
      affected_variables_.push_back(lvar->id());
      affected_var_props_[lvar->id()] = lvar->propertyInfo();
    }
  }
  // Sort affected variables list by packer index
  std::sort(affected_variables_.begin(),
            affected_variables_.end(),
            [&](const QUuid &lhs, const QUuid &rhs) {
              return affected_var_props_[lhs].index < affected_var_props_[rhs].index;
            });
}

string PackerConstraint::name() {
  return "PackerConstraint";
}

bool PackerConstraint::CaseSatisfiesConstraint(Optimization::Case *c) {
  for (auto id : affected_variables_) {
    if (c->get_real_variable_value(id) > 1.0 || c->get_real_variable_value(id) < 0.0) {
      return false;
    }
  }
  return true;
}

void PackerConstraint::SnapCaseToConstraints(Optimization::Case *c) {
  // Snap to upper/lower bounds
  for (auto id : affected_variables_) {
    if (c->get_real_variable_value(id) > 1.0) {
      c->set_real_variable_value(id, 1.0);
      if (vp_.vOPT >= 2) {
        ext_info("Snapped value to upper bound.", md_, cl_, vp_.lnw);
      }
    } else if (c->get_real_variable_value(id) < 0.0) {
      c->set_real_variable_value(id, 0.0);
      if (vp_.vOPT >= 2) {
        ext_info("Snapped value to lower bound.", md_, cl_, vp_.lnw);
      }
    }
  }

  // Enforce packer-ordering
  for (int i = 1; i < affected_variables_.size(); ++i) {
    if (c->get_real_variable_value(affected_variables_[i]) < c->get_real_variable_value(affected_variables_[i-1])) {
      c->set_real_variable_value(affected_variables_[i], c->get_real_variable_value(affected_variables_[i-1]));
      if (vp_.vOPT >= 2) {
        ext_info("Enforced packer-ordering.", md_, cl_, vp_.lnw);
      }
    }
  }
}

bool PackerConstraint::IsBoundConstraint() const {
  return true;
}

Eigen::VectorXd PackerConstraint::GetLowerBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd lbounds(id_vector.size());
  lbounds.fill(0);
  for (auto id : affected_variables_) {
    lbounds[id_vector.indexOf(id)] = 0.0;
  }
  return lbounds;
}

Eigen::VectorXd PackerConstraint::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);
  for (auto id : affected_variables_) {
    ubounds[id_vector.indexOf(id)] = 1.0;
  }
  return ubounds;
}

}
}
