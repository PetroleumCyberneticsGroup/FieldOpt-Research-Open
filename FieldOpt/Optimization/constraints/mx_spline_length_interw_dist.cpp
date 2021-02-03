/***********************************************************
Copyright (C) 2015-2017
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

#include "mx_spline_length_interw_dist.h"
#include <iostream>

namespace Optimization {
namespace Constraints {

MxSplineLengthInterwDist::
MxSplineLengthInterwDist(Settings::Optimizer::Constraint settings,
                         Model::Properties::VarPropContainer *variables,
                         Settings::VerbParams vp) : Constraint(vp) {

  max_iterations_ = settings.max_iterations;
  Settings::Optimizer::Constraint dist_constr_settings;
  dist_constr_settings.wells = settings.wells;
  dist_constr_settings.min = settings.min_distance;
  distance_constraint_ = new InterwDist(dist_constr_settings, variables, vp);

  if (vp_.vOPT >= 2){
    std::cout << "... ... initialized distance constraint for wells: ";
    for (QString wname : settings.wells) {
      std::cout << wname.toStdString() << ", ";
    }
    std::cout << std::endl;
  }

  length_constraints_ = QList<WSplineLength *>();
  for (QString wname : settings.wells) {
    Settings::Optimizer::Constraint len_constr_settings;
    len_constr_settings.well = wname;
    len_constr_settings.min = settings.min_length;
    len_constr_settings.max = settings.max_length;
    length_constraints_.append(new WSplineLength(len_constr_settings, variables, vp));
  }
}

bool MxSplineLengthInterwDist::CaseSatisfiesConstraint(Case *c) {
  if (!distance_constraint_->CaseSatisfiesConstraint(c))
    return false;
  for (WSplineLength *wsl : length_constraints_) {
    if (!wsl->CaseSatisfiesConstraint(c))
      return false;
  }
  return true;
}

void MxSplineLengthInterwDist::SnapCaseToConstraints(Case *c) {
  for (int i = 0; i < max_iterations_; ++i) {
    if (CaseSatisfiesConstraint(c)) {
      return;
    } else {
      distance_constraint_->SnapCaseToConstraints(c);
      for (WSplineLength *wsl : length_constraints_) {
        wsl->SnapCaseToConstraints(c);
      }
    }
  }
}

}
}
