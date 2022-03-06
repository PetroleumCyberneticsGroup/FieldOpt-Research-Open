/***********************************************************
Copyright (C) 2015-2016
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

#ifndef CONSTRAINTHANDLER_H
#define CONSTRAINTHANDLER_H
#include "constraint.h"
#include "bhp_constraint.h"
#include "well_spline_length.h"
#include "interw_dist.h"
#include "mx_spline_length_interw_dist.h"
#include "mx_spline_length_interw_dist_res_bound.h"
#include "polar_azimuth.h"
#include "polar_elevation.h"
#include "polar_well_length.h"
#include "polar_spline_boundary.h"
#include "polar_xyz_boundary.h"
#include "reservoir_boundary.h"
#include "pseudo_cont_boundary_2d.h"
#include "reservoir_boundary_toe.h"
#include "icv_constraint.h"
#include "packer_constraint.h"
#include "rate_constraint.h"
#include "Optimization/case.h"
#include "reservoir_xyz_boundary.h"
#include "Model/properties/var_prop_container.h"

#include "Settings/optimizer.h"
#include "Reservoir/grid/grid.h"

#include <QList>

namespace Optimization {
namespace Constraints {

/*!
 * \brief The ConstraintHandler class facilitates the
 * initialization and usage of multiple constraints.
 */
class ConstraintHandler
{
 public:
  ConstraintHandler(Settings::Optimizer *opt_settings,
                    Model::Properties::VarPropContainer *variables,
                    Reservoir::Grid::Grid *grid);

  //!< Check if a Case satisfies _all_ constraints.
  bool CaseSatisfiesConstraints(Case *c);

  //!< Snap all variables to _all_ constraints.
  void SnapCaseToConstraints(Case *c);

  QList<Constraint *> constraints() const {
    return constraints_;
  }

  /*!
   * @brief Check whether any of the constraints
   * within are boundary constraints.
   */
  bool HasBoundaryConstraints() const;

  /*!
   * @brief Initialize the normalizers for all constraints.
   * @param cases Cases to be used for determining parameters.
   */
  void InitializeNormalizers(QList<Case *> cases);

  /*!
   * @brief Get the sum of all normalized penalties
   * multiplied by their respective weights.
   * @param c The case to get the penalties for.
   * @return Weighted sum of all normalized penalties
   */
  long double GetWeightedNormalizedPenalties(Case *c);

  Eigen::VectorXd GetLowerBounds(QList<QUuid> id_vector) const;
  Eigen::VectorXd GetUpperBounds(QList<QUuid> id_vector) const;

 private:
  QList<Constraint *> constraints_;
  QList<Settings::Optimizer::Constraint> constraint_set_;
  Model::Properties::VarPropContainer *vars_;

  Settings::VerbParams vp_;
  string md_ = "Optimization/constraints";
  string cl_ = "ConstraintHandler";
  string im_ = "", wm_ = "", em_ = "";


};

}
}

#endif // CONSTRAINTHANDLER_H
