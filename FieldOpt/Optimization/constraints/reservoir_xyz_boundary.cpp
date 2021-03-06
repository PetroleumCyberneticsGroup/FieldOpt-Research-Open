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

#include "reservoir_xyz_boundary.h"
#include "ConstraintMath/well_constraint_projections/well_constraint_projections.h"
#include <iomanip>
#include "Utilities/math.hpp"

namespace Optimization {
namespace Constraints {

using Printer::ext_warn;

ReservoirXYZBoundary::ReservoirXYZBoundary(SO& seto, VPC *vars,
                                           Reservoir::Grid::Grid *grid, SV vp)
  : Constraint(seto, vars, vp) {

  if (vp_.vOPT >= 3) {
    im_ = "Adding ReservoirXYZBoundary constraint for ";
    im_ += seto_.well.toStdString();
    ext_info(im_, md_, cl_, vp_.lnw);
  }

  xmin_ = seto_.box_xyz_xmin;
  xmax_ = seto_.box_xyz_xmax;
  ymin_ = seto_.box_xyz_ymin;
  ymax_ = seto_.box_xyz_ymax;
  zmin_ = seto_.box_xyz_zmin;
  zmax_ = seto_.box_xyz_zmax;

  assert(xmin_ < xmax_);
  assert(ymin_ < ymax_);
  assert(zmin_ < zmax_);
  assert(!seto_.wells.empty());

  grid_ = grid;
  penalty_weight_ = seto_.penalty_weight;

  if (!vars->GetWSplineVars(seto_.well).empty()) {
    auto wspline_vars = vars->GetWSplineVars(seto.well);

    for (auto var : wspline_vars) {
      if (var->propertyInfo().coord == Model::Properties::Property::x) {
        var->setBounds(xmin_, xmax_);
        wspline_cnstrnd_uuid_vars_.push_back(var->id());
        im_ += var->name().toStdString() + ", ";

      } else if (var->propertyInfo().coord == Model::Properties::Property::y) {
        var->setBounds(ymin_, ymax_);
        wspline_cnstrnd_uuid_vars_.push_back(var->id());
        im_ += var->name().toStdString() + ", ";

      } else if (var->propertyInfo().coord == Model::Properties::Property::z) {
        var->setBounds(zmin_, zmax_);
        wspline_cnstrnd_uuid_vars_.push_back(var->id());
        im_ += var->name().toStdString() + ", ";

      }
    }

    if(seto_.scaling_) {
      xmin_ = -1.0;
      xmax_ = 1.0;
      ymin_ = -1.0;
      ymax_ = 1.0;
      zmin_ = -1.0;
      zmax_ = 1.0;
    }

    box_xyz_cnstrnd_well_ = initWSplineConstraint(wspline_vars, vp);

  } else if (vp_.vOPT >= 4) {
    wm_ = "GetWSplineVars for well ";
    wm_ += seto_.well.toStdString() + " is empty.";
    ext_warn(wm_, md_, cl_, vp_.lnw);
  }

  if(!wspline_cnstrnd_uuid_vars_.empty()) { isEnabled_ = true; }
}

bool ReservoirXYZBoundary::CaseSatisfiesConstraint(Case *c) {

  double heel_x_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.heel.x);
  double heel_y_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.heel.y);
  double heel_z_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.heel.z);

  double toe_x_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.toe.x);
  double toe_y_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.toe.y);
  double toe_z_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.toe.z);

  bool heel_feasible = false;
  bool toe_feasible = false;

  for (int ii = 0; ii < index_list_.length(); ii++) {
    if (grid_->GetCell(index_list_[ii]).EnvelopsPoint(
      Eigen::Vector3d(heel_x_val, heel_y_val, heel_z_val))) {
      heel_feasible = true;
    }
    if (grid_->GetCell(index_list_[ii]).EnvelopsPoint(
      Eigen::Vector3d(toe_x_val, toe_y_val, toe_z_val))) {
      toe_feasible = true;
    }
  }

  return heel_feasible && toe_feasible;
}

void ReservoirXYZBoundary::SnapCaseToConstraints(Case *c) {
  // dbg
  string s1 = DBG_prntVecXd(c->GetVarVector(wspline_cnstrnd_uuid_vars_));
  DBG_SnapCase(1, c->id_stdstr(), s1);

  double heel_x_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.heel.x);
  double heel_y_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.heel.y);
  double heel_z_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.heel.z);

  double toe_x_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.toe.x);
  double toe_y_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.toe.y);
  double toe_z_val = c->get_real_variable_value(box_xyz_cnstrnd_well_.toe.z);

  Eigen::Vector3d projected_heel =
    WellConstraintProjections::well_domain_constraint_indices(
      Eigen::Vector3d(heel_x_val, heel_y_val, heel_z_val),
      grid_,index_list_);

  Eigen::Vector3d projected_toe =
    WellConstraintProjections::well_domain_constraint_indices(
      Eigen::Vector3d(toe_x_val, toe_y_val, toe_z_val),
      grid_,index_list_);

  c->set_real_variable_value(box_xyz_cnstrnd_well_.heel.x, projected_heel(0));
  c->set_real_variable_value(box_xyz_cnstrnd_well_.heel.y, projected_heel(1));
  c->set_real_variable_value(box_xyz_cnstrnd_well_.heel.z, projected_heel(2));

  c->set_real_variable_value(box_xyz_cnstrnd_well_.toe.x, projected_toe(0));
  c->set_real_variable_value(box_xyz_cnstrnd_well_.toe.y, projected_toe(1));
  c->set_real_variable_value(box_xyz_cnstrnd_well_.toe.z, projected_toe(2));

}

Eigen::VectorXd ReservoirXYZBoundary::GetLowerBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd lbounds(id_vector.size());
  lbounds.fill(0);

  int ind_heel_x = id_vector.indexOf(box_xyz_cnstrnd_well_.heel.x);
  int ind_heel_y = id_vector.indexOf(box_xyz_cnstrnd_well_.heel.y);
  int ind_heel_z = id_vector.indexOf(box_xyz_cnstrnd_well_.heel.z);
  int ind_toe_x = id_vector.indexOf(box_xyz_cnstrnd_well_.toe.x);
  int ind_toe_y = id_vector.indexOf(box_xyz_cnstrnd_well_.toe.y);
  int ind_toe_z = id_vector.indexOf(box_xyz_cnstrnd_well_.toe.z);

  cout << "ind_heel_x: " << ind_heel_x << endl;
  cout << "ind_heel_y: " << ind_heel_y << endl;
  cout << "ind_heel_z: " << ind_heel_z << endl;

  lbounds(ind_heel_x) = xmin_;
  lbounds(ind_toe_x) = xmin_;
  lbounds(ind_heel_y) = ymin_;
  lbounds(ind_toe_y) = ymin_;
  lbounds(ind_heel_z) = zmin_;
  lbounds(ind_toe_z) = zmin_;

  return lbounds;
}

Eigen::VectorXd ReservoirXYZBoundary::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);

  int ind_heel_x = id_vector.indexOf(box_xyz_cnstrnd_well_.heel.x);
  int ind_heel_y = id_vector.indexOf(box_xyz_cnstrnd_well_.heel.y);
  int ind_heel_z = id_vector.indexOf(box_xyz_cnstrnd_well_.heel.z);
  int ind_toe_x = id_vector.indexOf(box_xyz_cnstrnd_well_.toe.x);
  int ind_toe_y = id_vector.indexOf(box_xyz_cnstrnd_well_.toe.y);
  int ind_toe_z = id_vector.indexOf(box_xyz_cnstrnd_well_.toe.z);

  ubounds(ind_heel_x) = xmax_;
  ubounds(ind_toe_x) = xmax_;
  ubounds(ind_heel_y) = ymax_;
  ubounds(ind_toe_y) = ymax_;
  ubounds(ind_heel_z) = zmax_;
  ubounds(ind_toe_z) = zmax_;
  return ubounds;
}

}
}
