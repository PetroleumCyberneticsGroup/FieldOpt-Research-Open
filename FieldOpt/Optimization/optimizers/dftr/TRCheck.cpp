/***********************************************************
Created by bellout on 21.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
TRStatus built from TrustRegionModel.cpp
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
// FINDBESTPOINT
int TRFrame::findBestPt() {

  int dim = (int)pts_abs_.rows();
  int n_points = (int)pts_abs_.cols();
  int best_i = 0;
  auto min_f = std::numeric_limits<double>::infinity();


  for (int k = 0; k < n_points; k++) {
    VectorXd point = pts_abs_.col(k);
    if (((point - lb_).minCoeff() >= 0) && ((ub_ - point).minCoeff() >= 0)) {
      auto val = fvals_(k);
      if (val < min_f) {
        min_f = val;
        best_i = k;
      }
    }
  }
  return best_i;
}

// _________________________________________________________
// MOVETOBESTPOINT
void TRFrame::moveToBestPt() {
  auto best_i = findBestPt();
  if (best_i != tr_center_) {
    tr_center_ = best_i;
  }
}

// _________________________________________________________
// CRITSTEP
Optimization::Optimizers::critExecStat TRFrame::critStep(double rad_bf_crit_step) {
  // tr_crit_mu_:    factor b/e radius & criticality measure
  // tr_crit_omega_: factor used to reduce radius
  // tr_crit_beta_:  ensure final radius reduction is not drastic
  // tr_tol_rad_:    tolerance of TR algorithm
  int exit_flag = 0; // use?

  //! Run ensureImpr
  if (!isLambdaPoised() || isModOld()) {
    exit_flag = ensureImpr();
    dbg_->prntCritStep(1, exit_flag, impr_cases_.size(), 0, areImprPtsCd());
    if (!areImprPtsCd() && exit_flag < 3) {
      setIsImprNeeded(true);
      return Optimization::Optimizers::ONGOING;
    }
  }

  //! Get gradient (g) at TR center
  // Note: assuming center polynomial is at position "0"
  dbg_->prntPolys("criticalityStep", modeling_polys_[0]);
  auto mMatrix = getModMatrix(0);

  ///! Project gradient
  auto x_center = pts_abs_.col(tr_center_);
  auto g_proj = ub_.cwiseMin(lb_.cwiseMax(x_center - mMatrix.g)) - x_center;
  auto crit_measure = tr_crit_mu_ * g_proj;

  if (radius_ > crit_measure.minCoeff()) {
    radius_ *= tr_crit_omega_;

    //! Re-run ensureImpr
    if (!isLambdaPoised() || isModOld()) {
      exit_flag = ensureImpr();
      dbg_->prntCritStep(2, exit_flag, impr_cases_.size(), 0, areImprPtsCd());
      if (!areImprPtsCd() && exit_flag < 3) {
        setIsImprNeeded(true);
        // revert to prev rad, reimpl once returning fval had been evaluated
        radius_ = rad_bf_crit_step;
        return Optimization::Optimizers::ONGOING;
      }
    }

    //! Check radius against tols
    auto C = (radius_ < tr_tol_rad_);
    auto E = (tr_crit_beta_ * crit_measure.norm() < tr_tol_f_);
    auto F = (radius_ < 100 * tr_tol_rad_);
    dbg_->prntCritStep(3, 0, 0, 0, C, E, F);
    if (C || (E && F)) {
      // Note: condition structured as: (C || (E && F))
      // whereas in CG it is: (C || E && F) (CLion complains)

      // CG comment: Better break. Not the end of algorithm,
      // but satisfies stopping condition for outer algorithm
      // anyway...
      wm_ = "[@critStep] Stopping condition for outer algorithm.";
      ext_warn(wm_, md_, cl_);
    }
  }

  // Final radius is increased to avoid a drastic reduction
  double fin_rad = max(radius_, tr_crit_beta_ * crit_measure.norm());
  radius_ = min(rad_bf_crit_step, fin_rad);
  dbg_->prntCritStep(4, 0, 0, 0, F, F, F, rad_bf_crit_step, fin_rad, radius_);
  if (abs(rad_bf_crit_step - radius_) < eps_) {
    return Optimization::Optimizers::FAILED;
  } else {
    return Optimization::Optimizers::SUCCESS;
  }
}

// _________________________________________________________
// ISLAMBDAPOISED
bool TRFrame::isLambdaPoised() {

  int dim = (int)pts_abs_.rows();
  int points_num = (int)pts_abs_.cols();

  double pivot_threshold = tr_piv_thld_;
  bool result = false;

  if (!settings_->parameters().tr_basis.compare("dummy")) {
    result = true;

  } else if (points_num >= dim+1) {
    //!<Fully linear, already>
    result = true;
    //!<but lets double check>
//        if (pivot_values_.lpNorm<Infinity>() > settings_->parameters().tr_piv_thld) {
//            Printer::ext_warn("Low pivot values.", "Optimization", "TrustRegionOptimization");
//        }
  } else {
    result = false;
  }
  return result;
}

// _________________________________________________________
// ISMODCOMPLT
bool TRFrame::isModComplete() {
  int dim = (int)pts_abs_.rows();
  int points_num = (int)pts_abs_.cols();
  int max_terms = ((dim+1)*(dim+2))/2;
  auto is_mod_cmpl = (points_num >= max_terms);

  if (points_num > max_terms) {
    ext_warn("Too many points in TR model.", md_, cl_);
  }
  return is_mod_cmpl;
}

// _________________________________________________________
// ISMODOLD
bool TRFrame::isModOld() {
  VectorXd distance(pts_abs_.rows());
  distance = pts_abs_.col(0) - pts_abs_.col(tr_center_);

  auto dist = distance.lpNorm<Infinity>();
  auto tr_rad_fac = settings_->parameters().tr_rad_fac;
  auto is_old = (dist > tr_rad_fac);

  return is_old;
}

// _________________________________________________________
// UNSHIFTPT
VectorXd TRFrame::unshiftPt(VectorXd &x) {

  auto shift_center = pts_abs_.col(0);
  auto bl_shifted = lb_ - shift_center;
  auto bu_shifted = ub_ - shift_center;

  dbg_->prntVecXd(shift_center, "shift_cntr: ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(bl_shifted,   "bl_shifted: ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(bu_shifted,   "bu_shifted: ", "% 10.3e ", dbg_->fn_mdat_);

  Eigen::VectorXd shifted_x = x + shift_center;
  dbg_->prntVecXd(shifted_x,   "shifted_x0: ", "% 10.3e ", dbg_->fn_mdat_);

  shifted_x = shifted_x.cwiseMin(bu_shifted + shift_center).eval();
  dbg_->prntVecXd(shifted_x,   "shifted_x1: ", "% 10.3e ", dbg_->fn_mdat_);

  shifted_x = shifted_x.cwiseMax(bl_shifted + shift_center).eval();
  dbg_->prntVecXd(shifted_x,   "shifted_x2: ", "% 10.3e ", dbg_->fn_mdat_);

  return shifted_x;
}

// _________________________________________________________
// CHECKINTERP
double TRFrame::checkInterp() {

  double tol1, tol2;
  if (radius_ < 1e-3) {
    tol1 = 100 * eps_;
  } else {
    tol1 = 10 * eps_;
  }
  tol2 = 10 * sqrt(eps_);

  // Remove shift center from all points
  int points_num = (int)pts_abs_.cols();

  MatrixXd hh = pts_abs_;
  for (int ii=points_num-1; ii >= 0; ii--) {
    hh.col(ii) = hh.col(ii) - pts_abs_.col(tr_center_);
  }

  double cval, cdiff, cA;

  auto p = getModPolys();

  double max_diff = -1.0;
  for (int kk = 0; kk < fvals_.rows(); kk++) {
    for (int ii = 0; ii < points_num; ii++) {
      VectorXd point = hh.col(ii);

      cval = evaluatePoly(p[kk], point);

      cdiff = std::abs(fvals_(ii) - cval);
      if (cdiff > max_diff) {
        max_diff = cdiff;
      }

      // This should be written using only Eigen methods...
      cA = std::max(tol1 * fvals_.array().abs().maxCoeff(), tol2);
      if (std::abs(cdiff) > cA) {
        // TODO: removed this print just to avoid messing up output
        // Printer::ext_warn("cmg:tr_interpolation_error",
        //                  "Interpolation error");
      }
    }
  }
  return max_diff;
}

// _________________________________________________________
// MEASURECRIT
VectorXd TRFrame::measureCrit() {
  // g is just the gradient, measured on the tr_center
  dbg_->prntPolys("measureCrit", modeling_polys_[0]);
  auto mMatrix = getModMatrix(0);

  // Projected gradient
  VectorXd tr_center_pt = pts_abs_.col(tr_center_);
  VectorXd xmax = lb_.cwiseMax(tr_center_pt - mMatrix.g);
  VectorXd xmin = ub_.cwiseMin(xmax);
  VectorXd xdiff = xmin - tr_center_pt;

  dbg_->prntVecXd(mMatrix.g,    "mMat.g: ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(tr_center_pt, "ctr_pt: ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(xmax,         "xmax:   ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(xmin,         "xmin:   ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(xdiff,        "xdiff:  ", "% 10.3e ", dbg_->fn_mdat_);

  return xdiff;
}