/***********************************************************
Sun May 02 2021 01:07:08 week 18 CET+0200
Copyright (C) 2021-
Mathias Bellout <chakibbb.pcg@gmail.com>

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

#include "SQP_SNOPT.h"

namespace Optimization {
namespace Optimizers {

SQP_SNOPT::SQP_SNOPT(Settings::Optimizer *seto,
  Case *bscs,
  Model::Model *model,
  Simulation::Simulator *simulator,
  Logger *logger, CaseHandler *case_handler,
  ConstraintHandler *constraint_handler)
  : Optimizer(seto, bscs, model->variables(), model->grid(), logger,
    case_handler, constraint_handler) {

  model_ = model;
  simulator_ = simulator;
  vars_ = model->variables();
  bscs_ = bscs;

  SNOPTSolver_ = new SNOPTSolver(); // subproblem solver

  // runner_ = new Runner::EmbeddedRunner(model, simulator);

  if (enable_logging_) {
    logger_->AddEntry(this);
  }
}

// _________________________________________________________
// SQP SNOPT SOLVER
tuple<VectorXd, double, int>
SQP_SNOPT::solver(VectorXd x, VectorXd bl, VectorXd bu) {

  auto xdim = x.rows();

  Eigen::VectorXd x0(getXDim(), 1);
  double fval;
  int exitflag = -1;

  double c;
  Eigen::VectorXd g(xdim);
  Eigen::MatrixXd H(xdim, xdim);

  // Set prob name (incl. spec) to be solved
  Optimizer::EmbeddedProblem prob;
  prob.tcase_ = new Case(bscs_);
  prob.model_ = model_;
  prob.simulator_ = simulator_;

  // prob.setProbName("TRMod_A");
  // prob.setProbName("Rosenbrock"); // dbg
  prob.setProbName("SQP");

  // Set prob dims
  prob.setNunVars(xdim);
  prob.setNumLinConst(0);
  prob.setNunNnlConst(0);

  // Set init point, bounds on vars
  prob.setXInit(x0);
  prob.setXUb(bu);
  prob.setXLb(bl);

  VectorXd FUb(1 + prob.getNunNnlConst());
  VectorXd FLb(1 + prob.getNunNnlConst());

  // Bounds on objective of tr prob
  FUb(0) = 1e200;
  FLb(0) = -1e200;

  // Bounds on nonlinear c imposed on tr prob
  // NOTE: not active since prob.getNunNnlConst() = 0
  // FUb(1) = radius;
  // FLb(1) = 0;

  prob.setFUb(FUb);
  prob.setFLb(FLb);

  // TODO Missing spec for linear constraints

  SNOPTSolver_->setUpSNOPTSolver(prob);

  x = prob.getXSol();
  fval = prob.getFSol()(0);

  // Should be zero -- double check
  if (prob.getSNOPTExitCode() == 1) { exitflag = 0; }
  return make_tuple(x, fval, exitflag);
}

}
}
