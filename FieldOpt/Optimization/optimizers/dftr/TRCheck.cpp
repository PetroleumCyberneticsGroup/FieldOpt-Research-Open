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
// GETCRITGRAD
VectorXd TRFrame::getCritGrad() {
  dbg_->prntPolys("getCriticality [grad]", modeling_polys_[0]);
  auto mMat = getModMatrix(0);

  // Projected gradient? -> get page # from DFBook
  VectorXd ctr_pt = pts_abs_.col(tr_center_);
  VectorXd xmax = lb_.cwiseMax(ctr_pt - mMat.g);
  VectorXd xmin = ub_.cwiseMin(xmax);
  VectorXd xdiff = xmin - ctr_pt;

  dbg_->prntVecXd(mMat.g, "mMat.g: ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(ctr_pt, "ctr_pt: ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(xmax,   "xmax:   ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(xmin,   "xmin:   ", "% 10.3e ", dbg_->fn_mdat_);
  dbg_->prntVecXd(xdiff,  "xdiff:  ", "% 10.3e ", dbg_->fn_mdat_);

  return xdiff;
}

// _________________________________________________________
// CRITSTEP
bool TRFrame::critStep(double rad_bf_crit_step) {

  // factor b/e radius & criticality measure
  double mu = settings_->parameters().tr_crit_mu;

  // factor used to reduce radius
  double omega = settings_->parameters().tr_crit_omega;

  // Ensure the final radius reduction is not drastic
  double beta = settings_->parameters().tr_crit_beta;

  // Tolerance of TR algorithm
  double tol_radius = settings_->parameters().tr_tol_radius;
  double tol_f = settings_->parameters().tr_tol_f;

  int exit_flag = 0;

  if (!isLambdaPoised() || isModOld()) {
    exit_flag = ensureImpr();
    if (!areImprPtsCd()) {
      setIsImprNeeded(true);
      return true;
    }
  }

  // Get gradient (g) at TR center
  // Note: assuming center polynomial is at position "0"
  dbg_->prntPolys("criticalityStep", modeling_polys_[0]);
  auto mMatrix = getModMatrix(0);

  // Project gradient
  auto x_center = pts_abs_.col(tr_center_);
  auto g_proj = ub_.cwiseMin(lb_.cwiseMax(x_center - mMatrix.g)) - x_center;
  auto crit_measure = mu * g_proj;

  if (radius_ > crit_measure.minCoeff()) {
    radius_ *= omega;

    if (!isLambdaPoised() || isModOld()) {

      exit_flag = ensureImpr();
      if (!areImprPtsCd()) {
        setIsImprNeeded(true);
        return true;
      }
    }

    if ((radius_ < tol_radius) ||
      (beta * crit_measure.norm() < tol_f) &&
        (radius_ < 100 * tol_radius)) {

      // Note: condition structured as: ((A || C) && B)
      // whereas in CG it is: (A || C && B) (CLion complains)

      // CG comment: Better break. Not the end of algorithm,
      // but satisfies stopping condition for outer algorithm
      // anyway...

    }
  }

  // The final radius is increased to avoid a drastic reduction
  radius_ = std::min(rad_bf_crit_step,
                     std::max(radius_, beta * crit_measure.norm()));
  return false;
}

// _________________________________________________________
// ISLAMBDAPOISED
bool TRFrame::isLambdaPoised() {

  int dim = (int)pts_abs_.rows();
  int points_num = (int)pts_abs_.cols();

  double pivot_threshold = tr_pivot_thld_;
  bool result = false;

  if (!settings_->parameters().tr_basis.compare("dummy")) {
    result = true;

  } else if (points_num >= dim+1) {
    //!<Fully linear, already>
    result = true;
    //!<but lets double check>
//        if (pivot_values_.lpNorm<Infinity>() > settings_->parameters().tr_pivot_thld) {
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
  auto tr_rad_fac = settings_->parameters().tr_radius_fac;
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