/***********************************************************
Created by bellout on 25.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
TRSolve built from TrustRegionModel.cpp
Copyright (C) 2018
Thiago Lima Silva <thiagolims@gmail.com>
Caio Giuliani <caiogiuliani@gmail.com>
Mathias Bellout <chakibbb.pcg@gmail.com>
--

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

// _________________________________________________________
// POINTNEW
tuple<MatrixXd, RowVectorXd, bool>
TRFrame::ptNew(poly polynomial, VectorXd tr_center_point,
                  double radius_used, VectorXd bl, VectorXd bu,
                  double pivot_threshold) {

  VectorXd new_point_min(getXDim(), 1);
  VectorXd new_point_max(getXDim(), 1);
  MatrixXd new_points(getXDim(), 2);

  double pivot_min, pivot_max;
  RowVectorXd new_pivot_values(new_points.cols());

  int exitflag_min, exitflag_max;
  bool point_found;

  tie(new_point_min, pivot_min, exitflag_min) =
    minimizeTr(polynomial, tr_center_point,
               radius_used, bl, bu);

  auto polynomial_max = multiplyPoly(polynomial, -1);

  tie(new_point_max, pivot_max, exitflag_max) =
    minimizeTr(polynomial_max, tr_center_point,
               radius_used, bl, bu);

  if ((exitflag_min >= 0)
    && (abs(pivot_min) >= pivot_threshold)) {


    if ((exitflag_max >= 0)
      && (abs(pivot_max) >= abs(pivot_min))) {
      new_points.col(0) = new_point_max;
      new_points.col(1) = new_point_min;
      new_pivot_values << pivot_max, pivot_min;

    } else if ((exitflag_max >= 0)
      && (abs(pivot_max) >= pivot_threshold)) {
      new_points.col(0) = new_point_min;
      new_points.col(1) = new_point_max;
      new_pivot_values << pivot_min, pivot_max;

    } else {
      new_points.conservativeResize(dim_,1);
      new_points.col(0) = new_point_min;

      new_pivot_values.conservativeResize(new_points.cols());
      new_pivot_values << pivot_min;
    }
    point_found = true;

  } else if ((exitflag_max >= 0)
    && (abs(pivot_max) >= pivot_threshold)) {
    new_points.conservativeResize(dim_,1);
    new_points.col(0) = new_point_max;

    new_pivot_values.conservativeResize(new_points.cols());
    new_pivot_values << pivot_max;
    point_found = true;

  } else {
    point_found = false;
    new_points.conservativeResize(0, 0);
    new_pivot_values.conservativeResize(1);
    new_pivot_values << 0;
  }

  return make_tuple(new_points, new_pivot_values, point_found);
}

// _________________________________________________________
// MINIMIZETR
tuple<VectorXd, double, int>
  TRFrame::minimizeTr(poly p, VectorXd x_tr_center,
                      double radius, VectorXd bl, VectorXd bu) {

  VectorXd x(getXDim(), 1);
  VectorXd bl_tr(getXDim(), 1), bu_tr(getXDim(), 1);
  VectorXd bl_mod(getXDim(), 1), bu_mod(getXDim(), 1);

  Eigen::VectorXd x0(getXDim(), 1);
  double fval;

  double tol_norm, tol_arg, tol_tr, tol_eps;
  int exitflag = -1;

  // TR tol
  // CG: tol_arg = max(1, norm(x_tr_center, inf));
  tol_tr = 10*std::numeric_limits<double>::epsilon();
  if (x_tr_center.lpNorm<Infinity>() > 2.0) {
    // Fix to match matlab code
    tol_tr *= 2;
  }

  // TR bounds
  bl_tr = x_tr_center.array() - radius;
  bu_tr = x_tr_center.array() + radius;

  bl_mod = bl.cwiseMax(bl_tr);
  bl_mod = bl_mod.array() + tol_tr;

  bu_mod = bu.cwiseMin(bu_tr);
  bu_mod = bu_mod.array() - tol_tr;

  // Restoring feasibility at TR center
  for (int ii=0; ii < bl.rows(); ii++) {
    if (x_tr_center(ii) <= bl(ii)) {
      bl_mod(ii) = bl(ii);
    }
    if (x_tr_center(ii) >= bu(ii)) {
      bu_mod(ii) = bu(ii);
    }
  }

  x0 = x_tr_center;

  double c;
  Eigen::VectorXd g(p.dim);
  Eigen::MatrixXd H(p.dim, p.dim);

  dbg_->prntPolys("minimizeTr", p);
  tie(c, g, H) = coeffsToMatrices(p.dim, p.coeffs);

  // Set prob name (incl. spec) to be solved
  Optimizer::EmbeddedProblem prob;
  prob.setProbName("TRMod_A");
  // prob.setProbName("Rosenbrock"); // dbg

  prob.setSeto(settings_);

  // Set prob dims
  prob.setNunVars(getXDim());
  prob.setNumLinConst(0);
  prob.setNunNnlConst(0);

  // Set init point, bounds on vars
  prob.setXInit(x0);
  prob.setXUb(bu_mod);
  prob.setXLb(bl_mod);

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

  // Transfer tr model data
  prob.setTrC(c);
  prob.setTrG(g);
  prob.setTrH(H);

  SNOPTSolver_->setUpSNOPTSolver(prob);

  x = prob.getXSol();
  fval = prob.getFSol()(0);

  // Should be zero -- double check
  if (prob.getSNOPTExitCode() == 1) { exitflag = 0; }
  return make_tuple(x, fval, exitflag);
}

// _________________________________________________________
// SOLVETRSUBPROB
tuple<VectorXd, double> TRFrame::solveTrSubprob() {
  auto x_tr_center = pts_abs_.col(tr_center_);
  auto p = getModPolys()[0];
  dbg_->prntPolys("solveTrSubproblem-00", p);

  auto ps = shiftPoly(p, -x_tr_center); //!<Shift to origin>
  dbg_->prntPolys("solveTrSubproblem-01", ps);

  VectorXd trial_point(getXDim(), 1);
  double trial_fval;
  int exit_flag;

  std::tie(trial_point, trial_fval, exit_flag) = minimizeTr(ps, x_tr_center, radius_, lb_, ub_);
  double current_fval = fvals_(tr_center_);
  double trial_decrease = current_fval - trial_fval;

  if (current_fval <= trial_fval) {
    // Printer::ext_info("current_fval <= trial_fval", "Optimization", "TrustRegionModel"); //TODO: commenting this prints to clean up the mess in the output
  }

  VectorXd tol_interp(1);
  tol_interp(0) = 1e-8;
  tol_interp.cwiseMax(eps_);
  double val, error_interp;

  for (int ii=0; ii < getNPts(); ii++) {
    val = evaluatePoly(ps, pts_abs_.col(ii));
    error_interp = abs(val - fvals_[ii]);

    if (error_interp > tol_interp(0)) {
      // TODO: removed this error to avoid messing up with output
      // Printer::ext_info("error_interp > tol_interp", "Optimization", "TrustRegionModel");

    }
  }
  return make_tuple(trial_point, trial_decrease);
}
































































