/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2021 Mathias Bellout
<chakibbb-pcg@gmail.com>

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

using namespace Model::Properties;
using Printer::info;
using Printer::ext_info;
using Printer::num2str;

ICVConstraint::ICVConstraint(const Settings::Optimizer::Constraint& settings,
                             Model::Properties::VarPropContainer *variables,
                             Settings::VerbParams vp)
  : Constraint(vp) {
  min_ = settings.min;
  max_ = settings.max;

  string tm = "Adding ICV constraint with [min, max] = [";
  tm += num2str(min_, 5) + ", " + num2str(max_, 5);
  tm += "] for well " + settings.well.toStdString() + " with variables: ";

  for (auto var : *variables->GetContinuousVariables()) {
    auto lvar = var.second;
    if (lvar->propertyInfo().prop_type == Property::PropertyType::ICD
      && QString::compare(lvar->propertyInfo().parent_well_name, settings.well) == 0) {
      affected_variables_.push_back(lvar->id());
      tm += lvar->name().toStdString() + ", ";
    }
  }

  ext_info(tm, md_, cl_, vp_.lnw);
}

bool ICVConstraint::CaseSatisfiesConstraint(Optimization::Case *c) {
  for (auto id : affected_variables_) {
    if (c->get_real_variable_value(id) > max_ || c->get_real_variable_value(id) < min_) {
      return false;
    }
  }
  return true;
}

void ICVConstraint::SnapCaseToConstraints(Optimization::Case *c) {
  if (vp_.vOPT >= 1) {
    string tm = "Checking case with vars { " + eigenvec_to_str(c->GetRealVarVector()) + " } ";
    tm += "against constraint [" + num2str(min_, 5) + ", " + num2str(max_, 5) + "]";
    ext_info(tm, md_, cl_, vp_.lnw);
  }

  for (auto id : affected_variables_) {
    if (c->get_real_variable_value(id) > max_) {
      c->set_real_variable_value(id, max_);
      if (vp_.vOPT >= 1) { ext_info("Snapped value to upper bound.", md_, cl_, vp_.lnw); }

    } else if (c->get_real_variable_value(id) < min_) {
      c->set_real_variable_value(id, min_);
      if (vp_.vOPT >= 1) { ext_info("Snapped value to lower bound.", md_, cl_, vp_.lnw); }
    }
  }
}

bool ICVConstraint::IsBoundConstraint() const { true; }

Eigen::VectorXd ICVConstraint::GetLowerBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd lbounds(id_vector.size());
  lbounds.fill(0);
  for (auto id : affected_variables_) {
    lbounds[id_vector.indexOf(id)] = min_;
  }
  return lbounds;
}

Eigen::VectorXd ICVConstraint::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);
  for (auto id : affected_variables_) {
    ubounds[id_vector.indexOf(id)] = max_;
  }
  return ubounds;
}

string ICVConstraint::name() {
  return "ICVConstraint";
}

}
}
