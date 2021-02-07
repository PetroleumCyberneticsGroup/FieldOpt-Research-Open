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

#ifndef BOXCONSTRAINT_H
#define BOXCONSTRAINT_H

#include "constraint.h"
#include <QPair>

namespace Optimization {
namespace Constraints {

using Model::Properties::ContinuousProperty;

/*!
 * \brief The BhpConstraint class Represents a is a simple
 * max/min constraint for BHP values for a well.
 */
class BhpConstraint : public Constraint
{
 public:
  //!< This class' constructor should not be used directly.
  //!< The constructors of subclasses should be used.
  BhpConstraint(::Settings::Optimizer::Constraint settings,
                ::Model::Properties::VarPropContainer *variables,
                Settings::VerbParams vp);

  string name() override { return cl_; }

  bool CaseSatisfiesConstraint(Case *c) override;
  void SnapCaseToConstraints(Case *c) override;

  bool IsBoundConstraint() const override;
  Eigen::VectorXd GetLowerBounds(QList<QUuid> id_vector) const override;
  Eigen::VectorXd GetUpperBounds(QList<QUuid> id_vector) const override;

 private:
  double min_, max_;
  QStringList bhp_cnstrnd_well_nms_;
  QList<ContinuousProperty *> bhp_cnstrnd_real_vars_;
  QList<QUuid> bhp_cnstrnd_uuid_vars_;

  Model::Properties::VarPropContainer *variables_;

  string md_ = "Optimizer::constraints";
  string cl_ = "BhpConstraint";
};

}
}

#endif // BOXCONSTRAINT_H
