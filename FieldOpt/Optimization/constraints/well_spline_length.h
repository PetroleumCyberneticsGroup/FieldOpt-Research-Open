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

#ifndef WELLSPLINELENGTH_H
#define WELLSPLINELENGTH_H

#include "constraint.h"
#include "well_spline_constraint.h"

namespace Optimization {
namespace Constraints {

/*!
 * \brief The WSplineLength class defines a constraint on
 * the maximum and minimum length of a well defined by a
 * WellSpline. It uses the WellIndexCalculation library.
 */
class WSplineLength : public Constraint, WellSplineConstraint
{
 public:
  WSplineLength(SO& seto, VPC *vars, SV vp);

  string name() override { return "WSplineLength"; };

  // Constraint interface
 public:
  bool CaseSatisfiesConstraint(Case *c) override;
  void SnapCaseToConstraints(Case *c) override;

  void InitializeNormalizer(QList<Case *> cases) override;
  double Penalty(Case *c) override;
  long double PenaltyNormalized(Case *c) override;

 private:
  double min_length_;
  double max_length_;
  Well affected_well_;

  string md_ = "Optimization::Constraints";
  string cl_ = "WSplineLength";

};

}
}

#endif // WELLSPLINELENGTH_H
