/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2021 Mathias Bellout
<chakibbb-pcg@gmail.com>

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

using Printer::info;
using Printer::ext_info;

ConstraintHandler::ConstraintHandler(Settings::Optimizer *opt_settings,
                                     Model::Properties::VarPropContainer *variables,
                                     Reservoir::Grid::Grid *grid) {

  constraint_set_ = opt_settings->constraints();
  vp_ = opt_settings->verbParams();

  for (Settings::Optimizer::Constraint constraint : constraint_set_) {
    switch (constraint.type) {
      case Settings::Optimizer::ConstraintType::BHP: {
        if (VERB_OPT >= 1 || vp_.vOPT >= 1) info("Adding BHP constraint for " + constraint.well.toStdString());
        constraints_.append(new BhpConstraint(constraint, variables, vp_));
        break;
      }
      case Settings::Optimizer::ConstraintType::Rate: {
        if (VERB_OPT >= 1) Printer::info("Adding Rate constraint for " + constraint.well.toStdString());
        constraints_.append(new RateConstraint(constraint, variables, vp_));
        break;
      }
      case Settings::Optimizer::ConstraintType::WellSplineLength: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding WSplineLength constraint for " + cons.well.toStdString());
          constraints_.append(new WSplineLength(cons, variables, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::PolarAzimuth: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding PolarAzimuth constraint for " + cons.well.toStdString());
          constraints_.append(new PolarAzimuth(cons, variables, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::PolarElevation: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding PolarElevation constraint for " + cons.well.toStdString());
          constraints_.append(new PolarElevation(cons, variables, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::PolarWellLength: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding PolarWellLength constraint for " + cons.well.toStdString());
          constraints_.append(new PolarWellLength(cons, variables, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::PolarSplineBoundary: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding PolarSplineBoundary constraint for " + cons.well.toStdString());
          constraints_.append(new PolarSplineBoundary(cons, variables, grid, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::ICVConstraint: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding ICV constraint for " + cons.well.toStdString());
          constraints_.append(new ICVConstraint(cons, variables, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::PackerConstraint: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding Packer constraint for " + cons.well.toStdString());
          constraints_.append(new PackerConstraint(cons, variables, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::WellSplineInterwellDistance: {
        if (VERB_OPT >= 1) Printer::info("Adding WellSplineInterwellDistance constraint.");
        constraints_.append(new InterwDist(constraint, variables, vp_));
        break;
      }
      case Settings::Optimizer::ConstraintType::CombinedWellSplineLengthInterwellDistance: {
        if (VERB_OPT >= 1) Printer::info("Adding CombinedWellSplineLengthInterwellDistance constraint.");
        constraints_.append(new MxSplineLengthInterwDist(constraint, variables, vp_));
        break;
      }
      case Settings::Optimizer::ConstraintType::
        CombinedWellSplineLengthInterwellDistanceReservoirBoundary: {
        if (VERB_OPT >= 1) Printer::info("Adding CombinedWellSplineLengthInterwellDistanceReservoirBoundary constraint.");
        constraints_.append(new MxSplineLengthInterwDistResBound(constraint, variables, grid, vp_));
        break;
      }
      case Settings::Optimizer::ConstraintType::ReservoirBoundary: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding ResBoundary constraint for " + cons.well.toStdString());
          constraints_.append(new ResBoundary(cons, variables, grid, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::PolarXYZBoundary: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding PolarXYZBoundary constraint for " + cons.well.toStdString());
          constraints_.append(new PolarXYZBoundary(cons, variables, grid, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::ReservoirXYZBoundary: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding ReservoirXYZBoundary constraint for " + cons.well.toStdString());
          constraints_.append(new ReservoirXYZBoundary(cons, variables, grid, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::ReservoirBoundaryToe: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding ReservoirBoundaryToe constraint for " + cons.well.toStdString());
          constraints_.append(new ReservoirBoundaryToe(cons, variables, grid, vp_));
        }
        break;
      }
      case Settings::Optimizer::ConstraintType::PseudoContBoundary2D: {
        for (auto wname : constraint.wells) {
          auto cons = Settings::Optimizer::Constraint(constraint);
          cons.well = wname;
          if (VERB_OPT >= 1) Printer::info("Adding PseudoContBoundary2D constraint for " + cons.well.toStdString());
          constraints_.append(new PseudoContBoundary2D(cons, variables, grid, vp_));
        }
        break;
      }
      default:
        Printer::ext_warn("Constraint type not recognized.", "Optimization", "ConstraintHandler");
    }
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
  auto vec_before = c->GetRealVarVector();
  for (Constraint *constraint : constraints_) {
    constraint->SnapCaseToConstraints(c);
  }
  if (vec_before != c->GetRealVarVector()) {
    c->state.cons = Case::CaseState::ConsStatus::C_PROJECTED;
  }
  else {
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
    if (con->IsBoundConstraint()) {
      lbounds = lbounds + con->GetLowerBounds(id_vector);
    }
  }
  return lbounds;
}

Eigen::VectorXd ConstraintHandler::GetUpperBounds(QList<QUuid> id_vector) const {
  Eigen::VectorXd ubounds(id_vector.size());
  ubounds.fill(0);
  for (auto con : constraints_) {
    if (con->IsBoundConstraint()) {
      ubounds = ubounds + con->GetUpperBounds(id_vector);
    }
  }
  return ubounds;
}

void ConstraintHandler::InitializeNormalizers(QList<Case *> cases) {
  cout << "ConstraintHandler Initializing constraint normalizers" << endl;
  for (auto con : constraints_) {
    con->InitializeNormalizer(cases);
  }
}

long double ConstraintHandler::GetWeightedNormalizedPenalties(Case *c) {
  long double wnp = 0.0L;
  for (auto con : constraints_) {
    long double pen = con->PenaltyNormalized(c);
    if (VERB_OPT >= 3 && pen > 0) {
      Printer::ext_info("Penalty from constraint " + con->name()
                          + ": " + Printer::num2str(pen), "Optimization", "ConstraintHandler");
    }
    wnp += pen * con->GetPenaltyWeight();
  }
  if (VERB_OPT >= 2) {
    Printer::ext_info("Weighted, normalized penalty for case " + c->id().toString().toStdString()
                        + ": " + Printer::num2str(wnp), "Optimization", "ConstraintHandler");
  }
  return wnp;
}

}
}
