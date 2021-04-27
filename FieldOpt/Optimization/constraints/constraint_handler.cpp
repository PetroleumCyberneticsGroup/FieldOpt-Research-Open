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

#include "constraint_handler.h"
#include <iostream>
#include <Utilities/printer.hpp>
#include <Utilities/verbosity.h>

namespace Optimization {
namespace Constraints {

// using Printer::info;
// using Printer::ext_info;
// using Printer::ext_warn;
// using Printer::num2str;

ConstraintHandler::ConstraintHandler(Settings::Optimizer *opt_settings,
                                     Model::Properties::VarPropContainer *variables,
                                     Reservoir::Grid::Grid *grid) {

  constraint_set_ = opt_settings->constraints();
  vp_ = opt_settings->verbParams();
  vars_ = variables;

  for (Settings::Optimizer::Constraint &constraint : constraint_set_) {

    switch (constraint.type) {

      //! -- ICD CONSTRAINTS
      case Settings::Optimizer::ConstraintType::ICVConstraint: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new ICVConstraint(cons, vars_, vp_));
        }
        break;
      }

        //! -- BHP
      case Settings::Optimizer::ConstraintType::BHP: {
        constraints_.append(new BhpConstraint(constraint, vars_, vp_));
        break;
      }

        //! -- RATE
      case Settings::Optimizer::ConstraintType::Rate: {
        constraints_.append(new RateConstraint(constraint, vars_, vp_));
        break;
      }


      // !----


        //! -- RESERVOIR XYZ BOUNDARY
      case Settings::Optimizer::ConstraintType::ReservoirXYZBoundary: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new ReservoirXYZBoundary(cons, vars_, grid, vp_));
        }
        break;
      }

        //! -- RESERVOIR BOUNDARY
      case Settings::Optimizer::ConstraintType::ReservoirBoundary: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new ResBoundary(cons, variables, grid, vp_));
        }
        break;
      }



        //! -- WELL SPLINE LENGTH
      case Settings::Optimizer::ConstraintType::WSplineLength: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new WSplineLength(cons, variables, vp_));
        }
        break;
      }

        //! -- WELLSPLINE INTERWELL DISTANCE
      case Settings::Optimizer::ConstraintType::WSplineInterwDist: {
        constraints_.append(new InterwDist(constraint, variables, vp_));
        break;
      }

        //! -- WSPLINE + WSPLINE LENGTH + INTERWELL DISTANCE
      case Settings::Optimizer::ConstraintType::MxWSplineLengthInterwDist: {
        constraints_.append(new MxSplineLengthInterwDist(constraint, vars_, vp_));
        break;
      }

        //! -- WSPLINE + WSPLINE LENGTH + INTERW DISTANCE + RESERVOIR BOUND
      case Settings::Optimizer::ConstraintType::
        MxWSplineLengthInterwDistResBound: {
        constraints_.append(new MxSplineLengthInterwDistResBound(constraint, vars_, grid, vp_));
        break;
      }














        // !----


        //! -- POLAR AZIMUTH
      case Settings::Optimizer::ConstraintType::PolarAzimuth: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new PolarAzimuth(cons, vars_, vp_));
        }
        break;
      }

        //! -- POLAR ELEVATION
      case Settings::Optimizer::ConstraintType::PolarElevation: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new PolarElevation(cons, vars_, vp_));
        }
        break;
      }

        //! -- POLAR WELL LENGTH
      case Settings::Optimizer::ConstraintType::PolarWellLength: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new PolarWellLength(cons, vars_, vp_));
        }
        break;
      }

        //! -- POLAR SPLINE BOUNDARY
      case Settings::Optimizer::ConstraintType::PolarSplineBoundary: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new PolarSplineBoundary(cons, vars_, grid, vp_));
        }
        break;
      }

        //! -- POLAR XYZ BOUNDARY
      case Settings::Optimizer::ConstraintType::PolarXYZBoundary: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new PolarXYZBoundary(cons, vars_, grid, vp_));
        }
        break;
      }

        //! -- PACKER CONSTRAINTS
      case Settings::Optimizer::ConstraintType::PackerConstraint: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new PackerConstraint(cons, vars_, vp_));
        }
        break;
      }

        //! -- RESERVOIR BOUNDARY TOE
      case Settings::Optimizer::ConstraintType::ReservoirBoundaryToe: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new ReservoirBoundaryToe(cons, vars_, grid, vp_));
        }
        break;
      }

        //! -- PSEUDO CONTINUOUS BOUNDARY 2D
      case Settings::Optimizer::ConstraintType::PseudoContBoundary2D: {
        for (auto &wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          constraints_.append(new PseudoContBoundary2D(cons, vars_, grid, vp_));
        }
        break;
      }

      default:
        ext_warn("Constraint type not recognized.", md_, cl_);
    }
  }

  if (opt_settings->ScaleVars() && !constraint_set_.empty()) {
    variables->ScaleVariables();
  }

}

bool ConstraintHandler::CaseSatisfiesConstraints(Case *c) {
  for (Constraint *constraint : constraints_) {
    if (!constraint->CaseSatisfiesConstraint(c)) {
      c->state.cons = Case::CaseState::ConsStatus::C_INFEASIBLE;
      return false;
    }
  }
  c->state.cons = Case::CaseState::ConsStatus::C_FEASIBLE;
  return true;
}

void ConstraintHandler::SnapCaseToConstraints(Case *c) {
  // \todo Case c may contain other vars besides real vars
  auto vec_before = c->GetRealVarVector();

  for (Constraint *con : constraints_) {
    if(con->isEnabled()) {
      con->SnapCaseToConstraints(c);
    }
  }

  if (vec_before != c->GetRealVarVector()) {
    c->state.cons = Case::CaseState::ConsStatus::C_PROJECTED;
  } else {
    c->state.cons = Case::CaseState::ConsStatus::C_FEASIBLE;
  }
}

bool ConstraintHandler::HasBoundaryConstraints() const {
  for (int i = 0; i < constraints_.size(); ++i) {
    if (constraints_[i]->IsBoundConstraint()) {
      return true;
    }
  }
  return false;
}

Eigen::VectorXd ConstraintHandler::GetLowerBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd lbounds(id_vector.size());
  lbounds.fill(0);
  for (auto con : constraints_) {
    if (con->IsBoundConstraint() && con->isEnabled()) {
      lbounds = lbounds + con->GetLowerBounds(id_vector);
    }
  }
  return lbounds;
}

Eigen::VectorXd ConstraintHandler::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);
  for (auto con : constraints_) {
    if (con->IsBoundConstraint() && con->isEnabled()) {
      ubounds = ubounds + con->GetUpperBounds(id_vector);
    }
  }
  return ubounds;
}

void ConstraintHandler::InitializeNormalizers(QList<Case *> cases) {
  if (vp_.vOPT >= 1) {
    info("ConstraintHandler Initializing constraint normalizers.");
  }
  for (auto con : constraints_) {
    con->InitializeNormalizer(cases);
  }
}

long double ConstraintHandler::GetWeightedNormalizedPenalties(Case *c) {
  long double wnp = 0.0L;

  for (auto con : constraints_) {
    long double pen = con->PenaltyNormalized(c);
    if (vp_.vOPT >= 3 && pen > 0) {
      im_ = "Penalty from constraint " + con->name() + ": " + num2str(pen);
      ext_info(im_, md_, cl_);
    }
    wnp += pen * con->GetPenaltyWeight();
  }

  if (vp_.vOPT >= 2) {
    im_ = "Weighted, normalized penalty for case ";
    im_ += c->id().toString().toStdString() + ": " + num2str(wnp);
    ext_info(im_, md_, cl_);
  }
  return wnp;
}

}
}
