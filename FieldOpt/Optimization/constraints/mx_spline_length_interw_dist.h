/***********************************************************
Copyright (C) 2015-2017
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

#ifndef COMBINEDSPLINELENGTHINTERWELLDISTANCE_H
#define COMBINEDSPLINELENGTHINTERWELLDISTANCE_H

#include "constraint.h"
#include "well_spline_length.h"
#include "interw_dist.h"

namespace Optimization {
namespace Constraints {

/*!
 * \brief The MxSplineLengthInterwDist class
 * combines the WSplineLength and InterwDist
 * constraints. It instantiates one well length constraint
 * per well and _one_ distance constraint.
 * All the constraints are applied cyclically in a loop
 * until either all constraints are satisfied or the
 * maximum number of iterations is reached.
 */
class MxSplineLengthInterwDist : public Constraint
{
 public:
  MxSplineLengthInterwDist(SO& seto, VPC *vars, SV vp);

  string name() override { return "MxSplineLengthInterwDist"; }

  // Constraint interface
 public:
  bool CaseSatisfiesConstraint(Case *c) override;
  void SnapCaseToConstraints(Case *c) override;

 private:
  int max_iterations_;
  QList<WSplineLength *> length_constraints_;
  InterwDist *distance_constraint_;
};

}}

#endif // COMBINEDSPLINELENGTHINTERWELLDISTANCE_H
