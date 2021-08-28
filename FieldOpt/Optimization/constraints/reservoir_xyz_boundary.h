/***********************************************************
Created by Brage on 31/03/19.
Copyright (C) 2019
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

#include "ConstraintMath/well_constraint_projections/well_constraint_projections.h"
#include "constraint.h"
#include "Reservoir/grid/grid.h"
#include "well_spline_constraint.h"

#ifndef FIELDOPT_RESERVOIR_XYZ_BOUNDARY_H
#define FIELDOPT_RESERVOIR_XYZ_BOUNDARY_H

namespace Optimization {
namespace Constraints {

using Model::Properties::ContinuousProperty;

/*!
 * @brief This class implements the box-constraint with the
 * WellSpline well definition. The purpose is to contain all
 * points of the well within the specified box.
 *
 * Note: Works well with the 3d box intersection generator
 * found in ResInsight @ ResInsight.org
 */

class ReservoirXYZBoundary : public Constraint, WellSplineConstraint
{
 public:
  ReservoirXYZBoundary(SO& seto, VPC *vars,
                       Reservoir::Grid::Grid *grid, SV vp);

  string name() override { return cl_; }

  bool CaseSatisfiesConstraint(Case *c) override;
  void SnapCaseToConstraints(Case *c) override;

  bool IsBoundConstraint() const override { return true; }
  Eigen::VectorXd GetLowerBounds(QList<QUuid> id_vector) const override;
  Eigen::VectorXd GetUpperBounds(QList<QUuid> id_vector) const override;

 protected:

  //!< Constraint limits for the box
  double xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;

  //!< Index list for the cells in the reservoir that are within the box
  QList<int> index_list_;

  //!< Grid, as defined in Reservoir/grid.h
  Reservoir::Grid::Grid *grid_;

  Well box_xyz_cnstrnd_well_;

  QList<ContinuousProperty *> wspline_cnstrnd_real_vars_;
  QList<QUuid> wspline_cnstrnd_uuid_vars_;

  string md_ = "Optimizer::constraints";
  string cl_ = "ReservoirXYZBoundary";

};
}
}
#endif //FIELDOPT_RESERVOIR_XYZ_BOUNDARY_H
