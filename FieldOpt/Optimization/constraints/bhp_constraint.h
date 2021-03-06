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
  BhpConstraint(SO& seto, VPC *vars, SV vp);

  string name() override { return cl_; }

  bool CaseSatisfiesConstraint(Case *c) override;
  void SnapCaseToConstraints(Case *c) override;

  bool IsBoundConstraint() const override  { return true; };
  Eigen::VectorXd GetLowerBounds(QList<QUuid> id_vector) const override;
  Eigen::VectorXd GetUpperBounds(QList<QUuid> id_vector) const override;

 private:
  QStringList bhp_cnstrnd_well_nms_;
  QList<ContinuousProperty *> bhp_cnstrnd_real_vars_;
  QList<QUuid> bhp_cnstrnd_uuid_vars_;

  string md_ = "Optimizer::constraints";
  string cl_ = "BhpConstraint";
};

}
}

#endif // BOXCONSTRAINT_H
