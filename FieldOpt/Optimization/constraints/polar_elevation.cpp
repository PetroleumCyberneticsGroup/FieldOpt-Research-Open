/***********************************************************
Copyright (C) 2019
Einar J.M. Baumann <einar.baumann@gmail.com>
Brage Strand Kristoffersen <brage_sk@hotmail.com>

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

#include "polar_elevation.h"
namespace Optimization {
namespace Constraints {

PolarElevation::PolarElevation(Settings::Optimizer::Constraint const settings,
                               Model::Properties::VarPropContainer *variables,
                               Settings::VerbParams vp)
                               : Constraint(vp) {

  if (vp_.vOPT >= 1) {
    info("Adding PolarElevation constraint for " + settings.well.toStdString());
  }

  min_elevation_ = settings.min;
  max_elevation_ = settings.max;
  for (auto var : *variables->GetContinuousVariables()){
    auto lvar = var.second;
    if (lvar->propertyInfo().polar_prop == Model::Properties::Property::PolarProp::Elevation
      && QString::compare(lvar->propertyInfo().parent_well_name, settings.well) == 0){
      affected_variable_ = lvar->id();
      break;
    }
  }
  if (affected_variable_.isNull()) {
    string em = "Affected variable null in PolarElevation constraint.";
    throw std::runtime_error(em);
  }
}

bool PolarElevation::CaseSatisfiesConstraint(Optimization::Case *c) {
  if (c->get_real_variable_value(affected_variable_) <= max_elevation_
    && c->get_real_variable_value(affected_variable_) >= min_elevation_){
    return true;
  } else {
    return false;
  }
}

void PolarElevation::SnapCaseToConstraints(Optimization::Case *c) {
  if (c->get_real_variable_value(affected_variable_) >= max_elevation_){
    c->set_real_variable_value(affected_variable_, max_elevation_);
  } else if (c->get_real_variable_value(affected_variable_) <= min_elevation_) {
    c->set_real_variable_value(affected_variable_, min_elevation_);
  }
}

bool PolarElevation::IsBoundConstraint() const {
  return true;
}

Eigen::VectorXd PolarElevation::GetLowerBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd lbounds(id_vector.size());
  lbounds.fill(0);
  lbounds[id_vector.indexOf(affected_variable_)] = min_elevation_;
  return lbounds;
}

Eigen::VectorXd PolarElevation::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);
  ubounds[id_vector.indexOf(affected_variable_)] = max_elevation_;
  return ubounds;
}

string PolarElevation::name() {
  return "PolarElevation";
}

}
}