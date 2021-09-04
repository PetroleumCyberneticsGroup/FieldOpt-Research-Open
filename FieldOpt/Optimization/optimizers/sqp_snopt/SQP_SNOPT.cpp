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
  Optimization::Objective::Objective *objf,
  Logger *logger, CaseHandler *case_handler,
  ConstraintHandler *constraint_handler)
  : Optimizer(seto, bscs, model->variables(), model->grid(), logger,
    case_handler, constraint_handler) {

  model_ = model;
  simulator_ = simulator;
  seto_ = seto;
  objf_ = objf;
  vars_ = model->variables();
  bscs_ = bscs;

  if (vp_.vOPT >= 1) { info("new SNOPTSolver()", vp_.lnw); }
  SNOPTSolver_ = new SNOPTSolver();

  if (enable_logging_) {
    logger_->AddEntry(this);
  }

  ub_.conservativeResize(bscs_->GetRealVarIdVector().size());
  lb_.conservativeResize(bscs_->GetRealVarIdVector().size());
  ub_.fill(seto->parameters().sqp_upper_bnd);
  lb_.fill(seto->parameters().sqp_lower_bnd);

  solver();
}

TC SQP_SNOPT::IsFinished() {
  // TC tc = NOT_FINISHED;
  // return tc;
}

void SQP_SNOPT::iterate() {
  if (enable_logging_) {
    logger_->AddEntry(this);
  }
  solver();
  iteration_++;
}

// _________________________________________________________
// SQP SNOPT SOLVER
tuple<VectorXd, double, int> SQP_SNOPT::solver() {

  auto xdim = bscs_->GetRealVarIdVector().size();

  Eigen::VectorXd x0(xdim, 1);
  double fval;
  VectorXd xstar;
  int exitflag = -1;

  double c;
  Eigen::VectorXd g(xdim);
  Eigen::MatrixXd H(xdim, xdim);

  // Set prob name (incl. spec) to be solved
  Optimizer::EmbeddedProblem prob;
  prob.tcase_ = new Case(bscs_);
  prob.model_ = model_;
  prob.simulator_ = simulator_;
  prob.objf_ = objf_;

  prob.setSeto(seto_);

  // prob.setProbName("TRMod_A");
  // prob.setProbName("Rosenbrock"); // dbg
  prob.setProbName("SQP");

  // Set prob dims
  prob.setNunVars(xdim);
  prob.setNumLinConst(0);
  prob.setNunNnlConst(0);

  // Set init point, bounds on vars
  cout << "bscs_->GetRealVarVector():\n" << bscs_->GetRealVarVector() << endl;
  prob.setXInit(bscs_->GetRealVarVector());
  prob.setXUb(ub_);
  prob.setXLb(lb_);

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

  xstar = prob.getXSol();
  fval = prob.getFSol()(0);

  // Should be zero -- double check
  if (prob.getSNOPTExitCode() == 1) { exitflag = 0; }
  return make_tuple(xstar, fval, exitflag);
}

}
}
