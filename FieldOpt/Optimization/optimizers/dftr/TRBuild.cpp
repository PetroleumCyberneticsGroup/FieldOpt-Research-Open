/***********************************************************
Created by bellout on 25.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
TRBuild built from TrustRegionModel.cpp
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
// CHANGETRCENTER
int TRFrame::changeTrCenter(VectorXd new_pt, double new_fval) {

  int pt_i = 0;
  bool point_added = false;
  bool point_exchanged = false;
  int exit_flag = 0;

  if (!isModComplete()) { // Add point
    exit_flag = addPt(new_pt, new_fval, tr_pivot_thld_);
    if(exit_flag > 0){ point_added = true; }
  }

  if (point_added) {
    // addPoint introduced point at last column.
    // Last column now set as new TR center.
    tr_center_ = (int)pts_abs_.cols() - 1;
    exit_flag = 1;

  } else {
    tie(point_exchanged, pt_i) =
      exchangePt(new_pt, new_fval, tr_pivot_thld_);

    if (point_exchanged) {
      tr_center_ = pt_i;
      exit_flag = 2;
    } else {
      // AddPoint and exchangePoint failed; however, we stil
      // need to add this new point as TR center. Model then
      // needs rebuilding
      int nc = (int)pts_abs_.cols();
      int nr = (int)pts_abs_.rows();
      pts_abs_.conservativeResize(nr, nc+1);
      pts_abs_.col(nc) = new_pt;

      int nc_f = (int)fvals_.cols();
      fvals_.conservativeResize(nc_f+1);
      fvals_(nc_f) = new_fval;

      tr_center_ = (int)pts_abs_.cols(); // Last point
      bool model_changed = rebuildModel();
      exit_flag = 4;
    }
  }
  return exit_flag;
}

// _________________________________________________________
// ADDPOINT
int TRFrame::addPt(VectorXd new_pt, double new_fval,
                   double rel_piv_thld) {
  int dim = (int)pts_abs_.rows();
  int last_pt = (int)pts_abs_.cols() - 1;
  double piv_thld = min(1.0, radius_)*rel_piv_thld;

  VectorXd new_pt_shftd(dim);
  VectorXd shift_center(dim);

  double piv_val = 0.0;
  bool success = false;
  int exit_flag;

  int polys_num = pivot_polys_.size();
  Matrix<double,Dynamic,Dynamic> pts_shftd_old = pts_shftd_;

  if (last_pt >= 0) {
    shift_center = pts_abs_.col(0);
    new_pt_shftd = new_pt - shift_center;

    pts_shftd_temp_ = pts_shftd_;
    int nc = pts_shftd_temp_.cols();
    pts_shftd_temp_.conservativeResize(
      pts_shftd_temp_.rows(), nc+1);
    pts_shftd_temp_.col(nc) = new_pt_shftd;

  } else {
    wm_ = "TR model had no point. Should this ever happen ?";
    ext_warn(wm_, md_, cl_);
  }

  int next_pos = last_pt + 1;
  if (next_pos == 0) {
    tr_center_ = 0; // Should be rebuilding model!
    // pivot_polys_ = bandPrioritizingBasis(dim);
    // TODO: implement a band prioritizing basis
    exit_flag = 1;
    im_ = "TR model had no point. This should never happen";
    ext_info(im_, md_, cl_);

  } else if (next_pos < polys_num) {

    int block_beg = 0;
    int block_end = 0;
    if (next_pos <= dim)  { // Add to linear block
      block_beg = 1;
      block_end = dim;
    } else { // Add to quadratic block
      block_beg = dim + 1;
      block_end = polys_num - 1;
    }

    tie(piv_val, success) =
      choosePivPoly(next_pos, block_end, piv_thld);

    if (success) {
      // Normalize poly value
      pivot_polys_[next_pos] =
        normalizePoly(next_pos, new_pt_shftd);

      // Re-orthogonalize
      pivot_polys_[next_pos] =
        orthogzToOthrPolys(next_pos, next_pos-1);

      // Orthogonalize polynomials on present
      // block (deferring subsequent ones)
      orthogzBlock(new_pt_shftd, next_pos,
                   block_beg, next_pos-1);
      exit_flag = 1;

    } else { exit_flag = 0; }
  } else {
    wm_ = "Full TR Model (next_pos >= polys_num). ";
    wm_ += "Remove another pt to add this one.";
    ext_warn(wm_, md_, cl_);
    exit_flag = 0;
  }

  if (exit_flag > 0) { // Update point matrices / vectors
    updatePts(new_pt, new_pt_shftd,
              new_fval, piv_val, next_pos);
  }
  return exit_flag;
}

// _________________________________________________________
// UPDATEPTS
void TRFrame::updatePts(VectorXd &new_pt, VectorXd &new_pt_shftd,
                        double &new_fval, double &piv_val, int &next_pos) {
  pts_shftd_temp_.col(next_pos) = new_pt_shftd;
  pts_shftd_.conservativeResize(pts_shftd_temp_.rows(),
                                pts_shftd_temp_.cols());
  pts_shftd_ = pts_shftd_temp_;

  fvals_.conservativeResize(fvals_.cols()+1);
  fvals_(next_pos) = new_fval;

  int nr = (int)pts_abs_.rows();
  int nc = (int)pts_abs_.cols();
  pts_abs_.conservativeResize(nr, nc+1);
  pts_abs_.col(next_pos) = new_pt;

  // dbg
  im_ = "mod_polys_.empty() - 00: % 10.0e ";
  dbg_->prntDbl((double)modeling_polys_.empty(),im_, dbg_->fn_pcfs_);
  im_ = "mod_polys_ vector b/f clear() 00:       ";
  dbg_->prntVecXd(modeling_polys_[0].coeffs,im_, "% 10.3e ", dbg_->fn_pcfs_);
  dbg_->prntPolys("addPoint b/f clear()             -- 00", modeling_polys_[0]); // dbg

  modeling_polys_.clear();

  im_ = "mod_polys_.empty() - 01: % 10.0e ";
  dbg_->prntDbl((double)modeling_polys_.empty(), im_, dbg_->fn_pcfs_);

  // bug: should not be resized:
  // pivot_values_.conservativeResize(pivot_values_.cols() + 1);
  pivot_values_(next_pos) = piv_val;
}

// _________________________________________________________
// EXCHANGEPT
tuple<bool,int> TRFrame::exchangePt(VectorXd new_pt,
                                    double new_fval,
                                    double rel_piv_thld) {

  bool succeeded = false;
  int pt_i = 0;
  int pivot_thld = min(double(1), radius_)*rel_piv_thld;
  int dim = (int)pts_abs_.rows();
  int last_p = (int)pts_abs_.cols()-1;
  int center_i = tr_center_;

  if (last_p >= 1) {
    auto shift_center = pts_abs_.col(0);
    auto new_pt_shftd = new_pt - shift_center;

    int block_beg = 0;
    int block_end = 0;
    if (last_p <= dim) {
      block_beg = 1;
      block_end = min(dim, last_p);
    } else {
      block_beg = dim+1;
      block_end = last_p;
    }

    dbg_->prntPivotPolys("exchangePoint-01");

    double max_val = 0;
    int max_poly_i = 0;
    for (int poly_i = block_end; poly_i >= block_beg; poly_i--) {
      if (poly_i != center_i) {
        double val = pivot_values_(poly_i);
        val *= evaluatePoly(pivot_polys_[poly_i], new_pt_shftd);
        if (abs(max_val) < abs(val)) {
          max_val = val;
          max_poly_i = poly_i;
        }
      }
    }

    double new_pivot_val = max_val;
    if (abs(new_pivot_val) > pivot_thld) {
      pts_shftd_.col(max_poly_i) = new_pt_shftd;

      // Normalize polynomial value
      pivot_polys_[max_poly_i] = normalizePoly(max_poly_i, new_pt_shftd);

      // Re-orthogonalize
      pivot_polys_[max_poly_i] = orthogzToOthrPolys(max_poly_i, last_p);

      // Orthogonalize polynomials on present block (until end)
      orthogzBlock(new_pt_shftd, max_poly_i, block_beg, last_p);

      // Update model
      int nc_cp = (int)cached_pts_.cols();
      cached_pts_.conservativeResize(dim, nc_cp+1);
      cached_pts_.col(nc_cp) = pts_abs_.col(max_poly_i);

      int nc_fv = (int)cached_fvals_.cols();
      cached_fvals_.conservativeResize(nc_fv+1);
      cached_fvals_(nc_fv) = fvals_(max_poly_i);

      pts_abs_.col(max_poly_i) = new_pt;
      fvals_(max_poly_i) = new_fval;
      pivot_values_(max_poly_i) = new_pivot_val;
      modeling_polys_.clear();
      dbg_->prntPolys("exchangePoint-02", modeling_polys_[0]);

      succeeded = true;
      pt_i = max_poly_i;
    } else {
      succeeded = false;
      pt_i = 0;
    }
    dbg_->prntFunctionData("exchangePoint", "", new_pt, shift_center);

  } else {
    succeeded = false;
    pt_i = 0;
    ext_info("cmg:runtime", md_, cl_);
  }
  return make_tuple(succeeded, pt_i);
}

// _________________________________________________________
// ORDERPTSPREBUILD
tuple<int, int> TRFrame::orderPtsPreBuild() {
  // Assemble all points available, incl. from cache
  all_pts_.conservativeResize(pts_abs_.rows(),
                              pts_abs_.cols() + cached_pts_.cols());
  if (cached_pts_.size() == 0) { all_pts_ = pts_abs_;
  } else { all_pts_ << pts_abs_, cached_pts_; }

  int dim = (int)all_pts_.rows();
  int n_points = (int)all_pts_.cols();

  // Assemble all function values, incl. from cache
  all_fvals_.conservativeResize(fvals_.cols() + cached_fvals_.cols());
  if (cached_fvals_.size() == 0) { all_fvals_ = fvals_;
  } else { all_fvals_ << fvals_, cached_fvals_; }

  if (tr_center_ != 0) { // Place center pt at 1st column
    all_pts_.col(0).swap(all_pts_.col(tr_center_));
    auto center_fvalue = all_fvals_(tr_center_);
    all_fvals_(tr_center_) = all_fvals_(0);
    all_fvals_(0) = center_fvalue;
  }

  // Calculate distances from TR center
  pts_shftd_.conservativeResize(dim, n_points);
  distances_.conservativeResize(n_points);
  pts_shftd_.setZero(dim, n_points);
  distances_.setZero(n_points);

  // Shift all points to TR center
  for (int i = 1; i < n_points; i++) { // Compute distances
    pts_shftd_.col(i) << all_pts_.col(i) - all_pts_.col(0);
    // Distances in infinity or 2-norm
    distances_(i) = pts_shftd_.col(i).lpNorm<Infinity>();
  }

  // Reorder points based on their distances to the TR center
  index_vector_.setLinSpaced(distances_.size(),
                             0, distances_.size() - 1);

  sort(index_vector_.data(),
       index_vector_.data() + index_vector_.size(),
       bind(comp,
            std::placeholders::_1,
            std::placeholders::_2,
            distances_.data()));

  // All points
  sortMatrixByIndex(pts_shftd_, index_vector_);
  sortMatrixByIndex(all_pts_, index_vector_);

  // All fvalues
  sortVectorByIndex(distances_, index_vector_);
  sortVectorByIndex(all_fvals_, index_vector_);

  return make_tuple(dim, n_points);
}

// _________________________________________________________
// REBUILDMODEL
bool TRFrame::rebuildModel() {
  int dim, n_points;
  tie(dim, n_points) = orderPtsPreBuild();
  dbg_->prntPivotPolys("rebuildModel");

  nfpBasis(dim); // build nfp polynomial basis>

  // Starting rowPivotGaussianElimination
  double pivot_thld = tr_pivot_thld_ * std::fmin(1, radius_);

  int polynomials_num = pivot_polys_.size();
  pivot_values_.conservativeResize(polynomials_num);
  pivot_values_.setZero(polynomials_num);

  // Constant term
  int last_pt_included = 0;
  int poly_i = 1;
  pivot_values_(0) = 1;

  // Gaussian elimination (using previous points)
  dbg_->prntModelData("rowPivotGaussianElim-00");

  for (int iter = 1; iter < polynomials_num; iter++) {

    pivot_polys_[poly_i] = orthogzToOthrPolys(poly_i, last_pt_included);

    double max_layer;
    auto farthest_point = (double)distances_(distances_.size()-1);
    double distance_farthest_point = (farthest_point/radius_);

    int block_beginning;
    int block_end;

    if (poly_i <= dim) {
      block_beginning = 1;
      block_end = dim;

      max_layer = std::min(2 * tr_rad_fac_, distance_farthest_point);
      if (iter > dim) {
        //!< We already tested all linear terms.
        //!< We do not have points to build a FL model.
        //!< How did this happen??? see Comment [1]>
        break;
      }

    } else { //!<Quadratic block -- being more carefull>
      max_layer = min(tr_rad_fac_, distance_farthest_point);
      block_beginning = dim + 1;
      block_end = polynomials_num - 1;
    }

    max_layer = std::fmax(1.0, max_layer);

    VectorXd all_layers;
    all_layers.setLinSpaced(ceil(max_layer), 1.0, max_layer);

    double max_absval = 0.0;
    double pt_max = 0.0;
    for (int i = 0; i < all_layers.size(); i++) {

      auto layer = all_layers(i);
      double dist_max = layer * radius_;
      for (int n = last_pt_included + 1; n < n_points; n++) {
        if (distances_(n) > dist_max) {
          break; //!<for n>
        }

        auto val = evaluatePoly(pivot_polys_[poly_i], pts_shftd_.col(n));
        val = val / dist_max; //!<minor adjustment>
        if (abs(max_absval) < abs(val)) {
          max_absval = val;
          pt_max = n;
        }
        if (abs(max_absval) > pivot_thld) {
          break; //!<for(layer)>
        }
      }
    }

    if (abs(max_absval) > pivot_thld) {
      //!<Points accepted>
      int pt_next = last_pt_included + 1;
      if (pt_next != pt_max) {
        pts_shftd_.col(pt_next).swap(pts_shftd_.col(pt_max));
        all_pts_.col(pt_next).swap(all_pts_.col(pt_max));
        std::swap(all_fvals_[pt_next], all_fvals_[pt_max]);
        std::swap(distances_[pt_next], distances_[pt_max]);
      }

      pivot_values_(pt_next) = max_absval;

      //!<Normalize polynomial value>
      auto point = pts_shftd_.col(pt_next);
      auto poly_normalized = normalizePoly(poly_i, point);

      pivot_polys_[poly_i] = poly_normalized;

      //!<Re-orthogonalize (just to make sure it still assumes 0 in previous points).
      //!< Unnecessary in infinite precision>
      pivot_polys_[poly_i] = orthogzToOthrPolys(poly_i, last_pt_included);

      //!<Orthogonalize polynomials on present block (deferring subsequent ones)>
      orthogzBlock(pts_shftd_.col(poly_i), poly_i, block_beginning, poly_i);

      last_pt_included = pt_next;
      poly_i++;
    } else {
      //!<These points don't render good pivot value for this>
      //!<specific polynomial. Exchange some polynomials>
      //!<and try to advance moving this polynomial>
      //!<to the end of the block>

      shiftPolyToEndBlock(poly_i, block_end);

      //!<Comment [1]: If we are on the linear block,>
      //!<this means we won't be able to build a Fully Linear model>
    }
  }

  dbg_->prntModelData("rowPivotGaussianElimination-01");

  tr_center_ = 0;
  pts_abs_ = all_pts_.leftCols(last_pt_included + 1).eval();
  pts_shftd_ = pts_shftd_.leftCols(last_pt_included + 1).eval();
  fvals_ = all_fvals_.head(last_pt_included + 1).eval();

  dbg_->prntModelData("rowPivotGaussianElimination-02");


  double cache_size = std::min(double(n_points - last_pt_included - 1), 3 * pow(dim, 2));

  // Consider removing, vector seems cleared already
  modeling_polys_.clear();
  //  DBG_printPolynomials("rowPivotGaussianElimination-02", modeling_polynomials_[0]);

  //!<Points not included>
  if (cache_size > 0) {
    cached_pts_ = all_pts_.middleCols(last_pt_included+1, cache_size).eval();
    cached_fvals_ = all_fvals_.segment(last_pt_included+1, cache_size).eval();

  } else {
    cached_pts_.conservativeResize(0, 0);
    cached_fvals_.conservativeResize(0);
  }

  //!<Clean auxiliary objects>
  all_pts_.conservativeResize(0, 0);
  all_fvals_.conservativeResize(0);

  dbg_->prntModelData("rowPivotGaussianElim-03");

  return last_pt_included < n_points; //!<model has changed>
}

// __________________________________________________________
// TRYTOADDPOINT
int TRFrame::tryToAddPt(VectorXd new_point, double new_fvalue) {
  int exit_flag = 0;
  bool point_added = false;

  if (!isModComplete()) { // Add this point
    exit_flag = addPt(new_point, new_fvalue, tr_add_thld_);
    if (exit_flag > 0){
      point_added = true;
    }
  }

  if (!point_added) { // Keep information about the new point
    int nc = cached_pts_.cols();
    int nr = max(cached_pts_.rows(), pts_abs_.rows());

    MatrixXd P(nr, nc+1);
    if (cached_pts_.cols() > 0) {
      P << cached_pts_, new_point;
    } else {
      P << new_point;
    }
    cached_pts_.conservativeResize(P.rows(), P.cols());
    cached_pts_.swap(P);

    int nc_f = cached_fvals_.cols();
    RowVectorXd F(nc_f+1);

    if (cached_fvals_.size() > 0) {
      F << cached_fvals_, new_fvalue;
    } else {
      F << new_fvalue;
    }

    cached_fvals_.conservativeResize(F.size());
    cached_fvals_.swap(F);

    exit_flag = 0;

  } else {
    exit_flag = 1;
  }
  return exit_flag;
}

// _________________________________________________________
// CHOOSENREPLACEPT
bool TRFrame::chooseNReplacePt() {
  double pivot_threshold = tr_exch_thld_;
  double radius = radius_;

  int dim = (int)pts_shftd_.rows();
  int points_num = (int)pts_shftd_.cols();
  auto points_shifted = pts_shftd_;

  int tr_center = tr_center_;
  auto shift_center = pts_abs_.col(0);
  auto tr_center_x = pts_shftd_.col(tr_center);

  auto pivot_values = pivot_values_;
  auto pivot_polynomials = pivot_polys_;

  int linear_terms = dim+1;

  auto bl_shifted = lb_ - shift_center;
  auto bu_shifted = ub_ - shift_center;
  bool success = false;

  dbg_->prntPivotPolys("chooseAndReplacePoint-00");
  dbg_->prntModelData("chooseAndReplacePoint-00");

  Eigen::Matrix<bool, 1, Dynamic> f_succeeded;
  f_succeeded.conservativeResize(1);
  f_succeeded.fill(false);

  //!<Reorder points based on their pivot_values>
  pivot_order_.setLinSpaced(pivot_values_.size(),
                            0, pivot_values_.size() - 1);
  dbg_->prntPivotVals("chRplcPnt");

  std::sort(pivot_order_.data(),
            pivot_order_.data() + pivot_order_.size(),
            std::bind(compAbs,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      pivot_values_.data()));
  dbg_->prntPivotVals("chRplcPnt");

  int polynomials_num = pivot_polys_.size();
  dbg_->prntPivotPolys("chooseAndReplacePoint-01");
  dbg_->prntModelData("chooseAndReplacePoint-01");

  //!<Could iterate through all pivots, but will
  //!< try just dealing with the worst one>
  auto pos = pivot_order_(0);

  if ((pos == 0) || (pos == tr_center) || (pos <= linear_terms && points_num > linear_terms)) {
    //!<Better to just rebuild model>
    success = false;
  } else {
    double new_pivot_value = 1;
    auto current_pivot_value = pivot_values_(pos);

    if (!areReplPtsCd()) {
      tie(repl_new_pts_shftd_, repl_new_pivots_, repl_pt_found_) =
        ptNew(pivot_polys_[pos], tr_center_x,
              radius_, bl_shifted, bu_shifted, pivot_threshold);
      dbg_->prntModelData("chooseAndReplacePoint-02");
    }

    if (repl_pt_found_) {
      for (int found_i = 0; found_i < repl_new_pts_shftd_.cols(); found_i++) {

        if (!areReplPtsCd()) {

          dbg_->prntModelData("chooseAndReplacePoint-03");
          dbg_->prntDbl((double)found_i, "found_i: % 10.0e ", dbg_->fn_mdat_);

          repl_new_pt_shftd_ = repl_new_pts_shftd_.col(found_i);
          auto new_pivot_value = repl_new_pivots_(found_i);

          repl_new_pt_abs_ = unshiftPt(repl_new_pt_shftd_);
          repl_new_fvals_.conservativeResize(repl_new_pt_abs_.cols());
          dbg_->prntModelData("chooseAndReplacePoint-04");

          setIsReplNeeded(true);

          for (int ii = 0; ii < repl_new_pt_abs_.cols(); ii++) {
            Case *new_case = new Case(base_case_);
            new_case->SetRealVarValues(repl_new_pt_abs_.col(ii));
            addReplCase(new_case);

            repl_pt_case_uuid_.push_back(new_case->id());
          }
          return true; //TODO Probably not necessary

        } else if (areReplPtsCd()) {

          dbg_->prntModelData("chooseAndReplacePoint-05");
          f_succeeded.conservativeResize(nfp_new_pt_abs_.cols(), 1);
          f_succeeded.fill(false);

          for (int ii = 0; ii < repl_new_pt_abs_.cols(); ii++) {

            int nvars = (int) repl_cases_[0]->GetRealVarVector().size();
            auto c = repl_cases_hash_[repl_pt_case_uuid_[ii]];
            repl_new_pt_abs_.col(ii) = c->GetRealVarVector();
            repl_new_fvals_(ii) = fmult_ * c->objf_value();

            std::map<string, string> stateMap = c->GetState();
            f_succeeded(ii) = true;
          }

          dbg_->prntModelData("chooseAndReplacePoint-06");
          setAreReplPtsCd(false);
          repl_pt_case_uuid_.clear();
        }
        if (f_succeeded.all()) {
          break;
        }
      }
    }

    if (repl_pt_found_ && f_succeeded.all()) {
      setIsReplNeeded(false);

      if (!repl_cases_.size() == 0) {
        clearReplCasesList();
      }

      //!<Normalize polynomial value>
      pivot_polys_[pos] = normalizePoly(pos, repl_new_pt_shftd_);
      dbg_->prntPivotPolys("chooseAndReplacePoint-03");
      dbg_->prntModelData("chooseAndReplacePoint-07");

      //!<Orthogonalize polynomials on present block (all)>
      int block_beginning;
      int block_end;
      if (pos <= dim) {
        block_beginning = 1;
        block_end = min(points_num-1, dim);
      } else {
        block_beginning = dim + 1;
        block_end = points_num-1;
      }

      //!<Re-orthogonalize>
      pivot_polys_[pos] = orthogzToOthrPolys(pos, block_end);

      //!<Orthogonalize block>
      orthogzBlock(repl_new_pt_shftd_, pos, block_beginning, block_end);

      //!<Update model and recompute polynomials>
      //TODO: check maximum cache size
      int nc_cp = (int)cached_pts_.cols();
      cached_pts_.conservativeResize(dim, nc_cp+1);
      cached_pts_.col(nc_cp) = pts_abs_.col(pos);

      int nc_fv = (int)cached_fvals_.cols();
      cached_fvals_.conservativeResize(nc_fv+1);
      cached_fvals_.col(nc_fv) = fvals_.col(pos);

      dbg_->prntModelData("chooseAndReplacePoint-08");
      pts_abs_.col(pos) = repl_new_pt_abs_;

      fvals_.col(pos) = repl_new_fvals_;

      pivot_values_((int)pos) *= new_pivot_value;

      modeling_polys_.clear();
      dbg_->prntPolys("chooseAndReplacePoint", modeling_polys_[0]);
      dbg_->prntModelData("chooseAndReplacePoint-09");

      success = true;
    }
  }
  return success;
}

// _________________________________________________________
// IMPROVEMODELNFP
bool TRFrame::improveModelNfp() {

  auto rel_pivot_threshold = settings_->parameters().tr_piv_thld;
  auto tol_radius = settings_->parameters().tr_rad_tol;
  auto radius_factor = settings_->parameters().tr_rad_fac;

  auto radius = radius_;
  auto pivot_threshold = rel_pivot_threshold*std::fmin(1, radius);
  auto dim = pts_shftd_.rows();
  auto p_ini = pts_shftd_.cols();
  auto shift_center = pts_abs_.col(0);

  auto tr_center = tr_center_;
  bool exit_flag = true;

  Eigen::Matrix<bool, 1, Dynamic> f_succeeded;
  f_succeeded.conservativeResize(1);
  f_succeeded.fill(false);

  int poly_i;
  double new_pivot_value = 0.0;

  dbg_->prntPivotPolys("improveModelNfp-00");
  dbg_->prntModelData("improveModelNfp-00");

  auto pivot_polynomials = pivot_polys_;
  auto polynomials_num = pivot_polynomials.size();

  if (lb_.rows() <= 0) {
    lb_.conservativeResize(dim);
    lb_.fill(-std::numeric_limits<double>::infinity());
  }

  if (ub_.rows() <= 0) {
    ub_.conservativeResize(dim);
    ub_.fill(std::numeric_limits<double>::infinity());
  }

  auto bl_shifted = lb_ - shift_center;
  auto bu_shifted = ub_ - shift_center;

  // Test if the model is already FL but old
  // Distance measured in inf norm
  auto tr_center_pt = pts_shftd_.col(tr_center);

  if ((tr_center_pt.lpNorm<Infinity>() > radius_factor * radius) && (p_ini > dim+1)) {
    exit_flag = false; // Needs to rebuild
  } else {
    // The model is not old>
    if (p_ini < polynomials_num) {
      int block_end = (int)p_ini;
      int block_begining = 0;
      // The model is not complete
      if (p_ini < dim+1) {
        // The model is not yet fully linear
        block_begining = 1;
      } else {
        // We can add a point to the quadratic block
        block_begining = (int)dim + 1;
      }
      int next_position = (int)p_ini;
      double radius_used = radius;

      // Possibly try smaller radius
      for (int attempts = 1; attempts<=3; attempts++) {

        // Iterate through available polynomials
        for (poly_i = next_position; poly_i<=block_end; poly_i++) {

          if(!areImprPtsCd()) {

            nfp_poly_ = orthogzToOthrPolys(poly_i, p_ini);
            std::tie(nfp_new_pts_shftd_, nfp_new_pivots_, nfp_pt_found_) =
              ptNew(nfp_poly_, tr_center_pt, radius_used, bl_shifted, bu_shifted, pivot_threshold);

          } else if(areImprPtsCd()) {
            // Using previously computed points/pivot-points
          }

          if (nfp_pt_found_) {
            for (int found_i = 0; found_i < nfp_new_pts_shftd_.cols(); found_i++) {

              if(!areImprPtsCd()) {

                nfp_new_pt_shftd_ = nfp_new_pts_shftd_.col(found_i);
                auto new_pivot_value = nfp_new_pivots_(found_i);

                nfp_new_pt_abs_ = unshiftPt(nfp_new_pt_shftd_);
                nfp_new_fvals_.conservativeResize(nfp_new_pt_abs_.cols());

                setIsImprNeeded(true);

                for (int ii = 0; ii < nfp_new_pt_abs_.cols(); ii++) {
                  Case *new_case = new Case(base_case_);
                  new_case->SetRealVarValues(nfp_new_pt_abs_.col(ii));
                  addImprCase(new_case);

                  pt_case_uuid_.push_back(new_case->id());
                }

                return 5;  //TODO Probably not necessary

              } else if (areImprPtsCd()) {
                f_succeeded.conservativeResize(nfp_new_pt_abs_.cols(),1);
                f_succeeded.fill(false);
                for (int ii = 0; ii < nfp_new_pt_abs_.cols(); ii++) {

                  auto c = impr_cases_hash_[pt_case_uuid_[ii]];
                  nfp_new_pt_abs_.col(ii) = c->GetRealVarVector();

                  if (settings_->mode() == Settings::Optimizer::OptimizerMode::Maximize) {
                    nfp_new_fvals_(ii) = -c->objf_value();
                  } else {
                    nfp_new_fvals_(ii) = c->objf_value();
                  }

                  std::map <string, string> stateMap = c->GetState();
                  // stateMap["EvalSt"] gives "pending" after evaluations of
                  // cases in tests...
                  // TODO Update status of cases when called from tests
                  // if(stateMap["EvalSt"] == "OKAY") {
                  f_succeeded(ii) = true;
                  // }
                }
                setAreImprPtsCd(false);
                pt_case_uuid_.clear();
              }
              if (f_succeeded.all()) {
                break;
              }
            }
            if (f_succeeded.all()) {
              break;
              // Stop trying pivot polynomials for poly_i
            }
          }
          // Attempt another polynomial if did not break
          if (poly_i < block_end) {
            nfp_new_pts_shftd_.conservativeResize(0,0);
            nfp_new_pivots_.conservativeResize(0);
            nfp_new_pt_shftd_.conservativeResize(0);
            nfp_new_pt_abs_.conservativeResize(0);
            nfp_pt_found_ = false;
            pt_case_uuid_.clear();
          }
        }

        if (nfp_pt_found_ && f_succeeded.all()) {
          break; // (for attempts)
        } else if (nfp_pt_found_) {
          // Reduce radius if it did not break
          radius_used = 0.5 * radius_used;
          if (radius_used < tol_radius) {
            break;
          }
        } else {
          break;
        }
      }

      if (nfp_pt_found_ && f_succeeded.all()) {
        setIsImprNeeded(false);

        if (!impr_cases_.size() == 0) {
          clearImprCasesList();
        }

        dbg_->prntPivotPolys("update-this-polynomial-00");
        // Update this polynomial in the set
        pivot_polys_[poly_i] = nfp_poly_;
        dbg_->prntPivotPolys("update-this-polynomial-01+swap-polynomials-00");

        // Swap polynomials
        std::swap(pivot_polys_[poly_i], pivot_polys_[next_position]);
        dbg_->prntPivotPolys("swaped-polynomials-01");

        // Add point
        int nr, nc;
        nr = (int)pts_shftd_.rows();
        nc = (int)pts_shftd_.cols();
        pts_shftd_.conservativeResize(nr, nc+1);
        pts_shftd_.col(nc) = nfp_new_pt_shftd_;

        // Normalize polynomial value
        pivot_polys_[next_position] = normalizePoly(next_position,
                                                    nfp_new_pt_shftd_);

        // Re-orthogonalize
        pivot_polys_[next_position] = orthogzToOthrPolys(next_position, p_ini);

        // Orthogonalize polynomials on present block (deffering subsequent ones)
        orthogzBlock(nfp_new_pt_shftd_, next_position, block_begining, p_ini);

        // Update model and recompute polynomials
        nr = (int)pts_abs_.rows();
        nc = (int)pts_abs_.cols();
        pts_abs_.conservativeResize(nr, nc+1);
        pts_abs_.col(nc) = nfp_new_pt_abs_;

        nr = (int)fvals_.rows();
        nc = (int)fvals_.cols();
        fvals_.conservativeResize(nr, nc+1);
        fvals_.col(nc) = nfp_new_fvals_;

        dbg_->prntPivotPolys("new_pivot_value-use-00");
        pivot_values_(next_position) = new_pivot_value;
        dbg_->prntPivotPolys("new_pivot_value-use-01");

        exit_flag = true;
      } else {
        exit_flag = false;
      }
    } else {
      // The model is already complete
      exit_flag = false;
    }
  }

  dbg_->prntPivotPolys("improveModelNfp-01");
  dbg_->prntModelData("improveModelNfp-01");
  return exit_flag;
}

// _________________________________________________________
// ENSUREIMPROVEMENT
int TRFrame::ensureImpr() {
  bool model_complete = isModComplete();
  bool model_fl = isLambdaPoised();
  bool model_old = isModOld();
  int exit_flag = 0;
  bool success = false;

  dbg_->prntModelData("EnsureImpr-00");

  if (!model_complete && (!model_old || !model_fl)) {
    //!<Calculate a new point to add>
    success = improveModelNfp(); //!<improve model>

    if ((success) || (!success && impr_cases_.size() > 0)) {
      exit_flag = 1;
    }
  } else if ((model_complete) && (!model_old)){
    success = chooseNReplacePt(); //!<replace point>
    if ((success) || (!success && repl_cases_.size() > 0))  {
      exit_flag = 2;
    }
  }
  dbg_->prntModelData("EnsureImpr-01");

  if ((!success) && (impr_cases_.size() == 0) && (repl_cases_.size() == 0)) {
    bool model_changed = rebuildModel();

    if (!model_changed) {
      if (!model_complete) {
        success = improveModelNfp();

      } else {
        //!<Replace point>
        success = chooseNReplacePt();
      }
    } else {
      success = true;
    }
    if (model_old) {
      exit_flag = 3;
    } else {
      exit_flag = 4;
    }
  }
  dbg_->prntModelData("EnsureImpr-02");
  return exit_flag;
}