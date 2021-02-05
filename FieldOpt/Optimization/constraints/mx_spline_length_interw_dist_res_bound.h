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

#ifndef COMBINED_SPLINE_LENGTH_INTERWELL_DISTANCE_RESERVOIR_BOUNDARY_H
#define COMBINED_SPLINE_LENGTH_INTERWELL_DISTANCE_RESERVOIR_BOUNDARY_H

#include "constraint.h"
#include "well_spline_length.h"
#include "interw_dist.h"
#include "reservoir_boundary.h"

namespace Optimization {
namespace Constraints {

/*!
 * \brief The MxSplineLengthInterwDistResBound class combines
 * the ResBoundary constraint, the WSplineLength as
 * well as the InterwDist constraint into a routine that checks
 * for these constraints in a sequential manner within a loop.
 *
 * The constraints are applied until all are satisfied
 *  or until a maximum number of iterations is reached.
 * The class instantiates one well length constraint per
 * well, _one_ distance constraint and _one_ reservoir
 * boundary constraint.
 */

class MxSplineLengthInterwDistResBound : public Constraint
{
 public:
  MxSplineLengthInterwDistResBound(
    Settings::Optimizer::Constraint settings,
    Model::Properties::VarPropContainer *variables,
    Reservoir::Grid::Grid *grid,
    Settings::VerbParams vp);

  bool IsBoundConstraint() const override;
  Eigen::VectorXd GetLowerBounds(QList<QUuid> id_vector) const override;
  Eigen::VectorXd GetUpperBounds(QList<QUuid> id_vector) const override;

  string name() override { return "MxSplineLengthInterwDistResBound"; }

  // Constraint interface
 public:
  bool CaseSatisfiesConstraint(Case *c);
  void SnapCaseToConstraints(Case *c);

 private:
  int max_iterations_;
  QList<WSplineLength *> length_constraints_;
  QList<ResBoundary *> boundary_constraints_;
  InterwDist *distance_constraint_;
};

}}
#endif // COMBINED_SPLINE_LENGTH_INTERWELL_DISTANCE_RESERVOIR_BOUNDARY_H
