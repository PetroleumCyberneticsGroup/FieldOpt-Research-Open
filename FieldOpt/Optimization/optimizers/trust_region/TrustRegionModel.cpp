/*********************************************************************
 Created by thiagols on 27.11.18
 Copyright (C) 2018 Thiago Lima Silva<thiagolims@gmail.com>
 Modified 2018-2019 Mathias Bellout <mathias.bellout@ntnu.no>

 This file is part of the FieldOpt project.

 FieldOpt is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published
 by the Free Software Foundation, either version 3 of the License,
 or (at your option) any later version.

 FieldOpt is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FieldOpt. If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#include "TrustRegionModel.h"
#include <Settings/optimizer.h>

#include "Utilities/printer.hpp"
#include "Utilities/stringhelpers.hpp"
#include <Utilities/verbosity.h>

#include <numeric>
#include <functional>
#include <limits>

using namespace Eigen;
using std::vector;
using std::swap;
using std::tuple;
using std::tie;
using std::bind;
using std::fmax;
using std::fmin;

namespace Optimization {
namespace Optimizers {

bool compare(
    const int &lhs,
    const int &rhs,
    const double *distances) {
  return distances[lhs] <= distances[rhs];
}

TrustRegionModel::TrustRegionModel(
    VectorXd& lb,
    VectorXd& ub,
    Case *base_case,
    Settings::Optimizer *settings) {

    init_points_computed_ = false;
    impr_points_computed_ = false;
    repl_points_computed_ = false;

    is_initialized_ = false;
    needs_improvement_ = false;
    needs_replacement_ = false;
    model_changed_ = false;

    lb_ = lb;
    ub_ = ub;
    base_case_ = base_case;

    settings_ = settings;
    radius_ = settings_->parameters().tr_initial_radius;
    tr_center_ = 0;
    cache_max_ = (int)3*pow(dim_,2);

    pivot_values_.resize(0);
    cached_fvalues_.resize(0);

    index_vector_.resize(0);
    distances_.resize(0);

    cached_points_.resize(0,0);
    points_shifted_.resize(0,0);

    // rebuildModel();
    // moveToBestPoint();
    // computePolynomialModels();

    // if (points_abs_.cols() < 2) {
    // ensureImprovement();
    // }

    SNOPTSolver_ = new SNOPTSolver();

}

void TrustRegionModel::moveToBestPoint() {
  auto best_i = findBestPoint();
  if (best_i != tr_center_) {
    tr_center_ = best_i;
  }

// Here should rebuild polynomials!!!
}


void TrustRegionModel::criticalityStep() {

    // factor b/e radius & criticality measure
    double mu = settings_->parameters().tr_criticality_mu;

    // factor used to reduce radius
    double omega = settings_->parameters().tr_criticality_omega;

    // Ensure the final radius reduction is not drastic
    double beta = settings_->parameters().tr_criticality_beta;

    // Tolerance of TR algorithm
    double tol_radius = settings_->parameters().tr_tol_radius;
    double tol_f = settings_->parameters().tr_tol_f;

    double init_radius_ = settings_->parameters().tr_initial_radius;

    while (!isLambdaPoised() || isOld()) {

        ensureImprovement();
        computePolynomialModels();
        bool model_changed = rebuildModel();

        if (!model_changed) {
            Printer::ext_info("[criticalityStep] Model did not change.",
                              "Optimization", "Trust Region Model");
            break;
        }
    }

    // Get gradient (g) at TR center
    Polynomial center_polynomial = modeling_polynomials_[0];
    // Note: assuming center polynomial is at position "0"

    double c;
    VectorXd g(center_polynomial.dimension);
    MatrixXd H(center_polynomial.dimension, center_polynomial.dimension);

    std::tie(c, g, H) = coefficientsToMatrices(
            center_polynomial.dimension, center_polynomial.coefficients);

    // Project gradient
    auto x_center = points_abs_.col(tr_center_);
    auto g_proj = ub_.cwiseMin(lb_.cwiseMax(x_center - g)) - x_center;
    auto crit_measure = mu * g_proj;

    while (radius_ > crit_measure.minCoeff()) {
        radius_ *= omega;

        while (!isLambdaPoised() || isOld()) {

            ensureImprovement();
            computePolynomialModels();
            bool model_changed = rebuildModel();

            if (!model_changed) {
                Printer::ext_info("[criticalityStep] Model did not change.",
                                  "Optimization", "Trust Region Model");
                break;
            }

        }

        if ((radius_ < tol_radius
        || beta * crit_measure.norm() < tol_f)
        && radius_ < 100 * tol_radius) {
            break;

            // Note: condition structured as: ((A || C) && B)
            // whereas in CG it is: (A || C && B) (CLion complains)

            // CG comment: Better break. Not the end of algorithm,
            // but satisfies stopping condition for outer algorithm
            // anyway...

        }
    }

    // The final radius is increased to avoid a drastic reduction
    radius_ = std::min(
            std::max(radius_, beta * crit_measure.norm()),
            init_radius_);
}

VectorXd TrustRegionModel::measureCriticality() {
  auto p = getModelingPolynomials()[0];
  double c;
  VectorXd g(p.dimension);
  MatrixXd H(p.dimension, p.dimension);

  //!<g is just the gradient, measured on the tr_center>
  tie(c, g, H) = coefficientsToMatrices(p.dimension, p.coefficients);

  //!<Projected gradient>
  VectorXd tr_center_pt = points_shifted_.col(tr_center_);
  VectorXd xmax = lb_.cwiseMax(tr_center_pt - g);
  VectorXd xmin = ub_.cwiseMin(xmax);
  return xmin - tr_center_pt;
}

double TrustRegionModel::checkInterpolation() {

  double tol1, tol2;
  if (radius_ < 1e-3) {
    tol1 = 100*std::numeric_limits<double>::epsilon();
  } else {
    tol1 = 10*std::numeric_limits<double>::epsilon();
  }
  tol2 = 10*sqrt(std::numeric_limits<double>::epsilon());

  // Remove shift center from all points
  int points_num = (int)points_abs_.cols();
  MatrixXd hh = points_abs_;

  for (int ii=points_num; ii==0; ii--) {
    hh.col(ii) = hh.col(ii) - points_abs_.col(tr_center_);
  }

  double cval, cdiff, conda;

  int max_diff = -1;
  for (int kk = 0; kk < dim_; kk++) {
    for (int ii = 0; ii < points_num; ii++) {

      cval = evaluatePolynomial(pivot_polynomials_[kk], hh.col(ii));
      cdiff = std::abs(fvalues_(ii) - cval);

      if (cdiff > max_diff) {
        max_diff = cdiff;
      }

      // This should be written using only Eigen methods...
      conda = std::max(
          tol1 * fvalues_.array().abs().maxCoeff(), tol2);

      if (std::abs(cdiff) > conda) {
        Printer::ext_warn("cmg:tr_interpolation_error",
                          "Interpolation error");

      }

    }
  }

  // Comment: from MATLAB code: fvalues_ should be a matrix...
  // difference = abs(model.fvalues(k, m) - this_value);
  // and
  // max(tol_1*max(abs(model.fvalues(k, :))), tol_2)
}

void TrustRegionModel::submitTempInitCases() {
  initialization_cases_ = temp_init_cases_;

  int nvars = initialization_cases_[0]->GetRealVarVector().size();
  points_abs_.setZero(nvars, initialization_cases_.size());
  fvalues_.setZero(initialization_cases_.size());

  int ii = 0;
  for (Case *c : initialization_cases_) {
    points_abs_.col(ii) = c->GetRealVarVector();
    fvalues_(ii) = c->objective_function_value();
    ii++;
  }
}

void TrustRegionModel::submitTempImprCases() {
    improvement_cases_ = temp_impr_cases_;

    for (Case *c : improvement_cases_) {
        improvement_cases_hash_.insert(c->id(), c);
    }
}

void TrustRegionModel::submitTempReplCases() {
  replacement_cases_ = temp_repl_cases_;

  for (Case *c : replacement_cases_) {
    replacement_cases_hash_.insert(c->id(), c);
  }
}

bool TrustRegionModel::rebuildModel() {
  double pivot_threshold = settings_->parameters().tr_pivot_threshold * std::fmin(1, radius_);

  //!<All points we know>
  all_points_.resize(points_abs_.rows(), points_abs_.cols() + cached_points_.cols());
  if (cached_points_.size() == 0) {
    all_points_ = points_abs_;
  } else {
    all_points_ << points_abs_;
    all_points_ << cached_points_;
  }

  //!<All function values we know>
  all_fvalues_.resize(fvalues_.size() + cached_fvalues_.size());

  if (cached_fvalues_.size() == 0) {
    all_fvalues_ = fvalues_;
  } else {
    all_fvalues_ << fvalues_;
    all_fvalues_ << cached_fvalues_;
  }

  int dim = all_points_.rows();
  int n_points = all_points_.cols();
  if (tr_center_ != 0) {
    //!<Center will be first>
    all_points_.col(0).swap(all_points_.col(tr_center_));
    auto center_fvalue = all_points_(tr_center_);
    all_fvalues_(tr_center_) = all_fvalues_(0);
    all_fvalues_(0) = center_fvalue;
  }
  //!<Calculate distances from tr center>
  points_shifted_.resize(dim, n_points);
  distances_.resize(n_points);
  points_shifted_.setZero(dim, n_points);
  distances_.setZero(n_points);

  //!<Shift all points to TR center>
  for (int i = 1; i < n_points; i++) {
    //!<Compute distances>
    points_shifted_.col(i) << all_points_.col(i) - all_points_.col(0);

    //<!distances in infinity or 2-norm>
    distances_(i) = points_shifted_.col(i).lpNorm<Infinity>();
  }

  //!<Reorder points based on their distances to the tr center>
  index_vector_.setLinSpaced(distances_.size(), 0, distances_.size() - 1);
  std::sort(index_vector_.data(),
          index_vector_.data() + index_vector_.size(),
          std::bind(compare,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  distances_.data()));

  //!<All points>
  sortMatrixByIndex(points_shifted_, index_vector_);
  sortMatrixByIndex(all_points_, index_vector_);

  //!<All fvalues>
  sortVectorByIndex(distances_, index_vector_);
  sortVectorByIndex(all_fvalues_, index_vector_);

  nfpBasis(dim);//!<build nfp polynomial basis>

  int polynomials_num = pivot_polynomials_.size();
  pivot_values_.resize(polynomials_num);
  pivot_values_.setZero(polynomials_num);

  //!<Constant term>
  int last_pt_included = 0;
  int poly_i = 1;
  pivot_values_(0) = 1;

  //!<Gaussian elimination (using previous points)>
  for (int iter = 1; iter < polynomials_num; iter++) {

    pivot_polynomials_[poly_i] =
            orthogonalizeToOtherPolynomials(poly_i, last_pt_included);

    double max_layer;
    double farthest_point = (double)distances_(distances_.size()-1);
    double distance_farthest_point = (farthest_point/radius_);

    int block_beginning;
    int block_end;

    if (poly_i <= dim) {
      block_beginning = 1;
      block_end = dim;
      //!<linear block -- we allow more points (*2)

      max_layer = std::min(2 * settings_->parameters().tr_radius_factor, distance_farthest_point);
      if (iter > dim) {
        //!< We already tested all linear terms.
        //!< We do not have points to build a FL model.
        //!< How did this happen??? see Comment [1]>
        break;

      }
    } else { //!<Quadratic block -- being more carefull>
      max_layer = min(settings_->parameters().tr_radius_factor, distance_farthest_point);
      block_beginning = dim + 1;
      block_end = polynomials_num - 1;
    }

    max_layer = std::fmax(1, max_layer);

    VectorXd all_layers;
    all_layers.setLinSpaced(ceil(max_layer), 1, max_layer);

    double max_absval = 0;
    double pt_max = 0;
    for (int i = 0; i < all_layers.size(); i++) {
      auto layer = all_layers(i);
      double dist_max = layer * radius_;
      for (int n = last_pt_included + 1; n < n_points; n++) {
        if (distances_(n) > dist_max) {
          break; //!<for n>
        }

        auto val = evaluatePolynomial(
                pivot_polynomials_[poly_i], points_shifted_.col(n));
        val = val / dist_max; //!<minor adjustment>
        if (abs(max_absval) < abs(val)) {
          max_absval = val;
          pt_max = n;
        }
        if (abs(max_absval) > pivot_threshold) {
          break; //!<for(layer)>
        }
      }
    }

    if (abs(max_absval) > pivot_threshold) {
      //!<Points accepted>
      int pt_next = last_pt_included + 1;
      if (pt_next != pt_max) {
        points_shifted_.col(pt_next).swap(points_shifted_.col(pt_max));
        all_points_.col(pt_next).swap(all_points_.col(pt_max));
        std::swap(all_fvalues_[pt_next], all_fvalues_[pt_max]);
        std::swap(distances_[pt_next], distances_[pt_max]);
      }

      pivot_values_(pt_next) = max_absval;

      //!<Normalize polynomial value>
      auto point = points_shifted_.col(pt_next);
      pivot_polynomials_[poly_i] = normalizePolynomial(poly_i, point);

      //!<Re-orthogonalize (just to make sure it still assumes 0 in previous points).
      //!< Unnecessary in infinite precision>
      pivot_polynomials_[poly_i] = orthogonalizeToOtherPolynomials(poly_i, last_pt_included);

      //!<Orthogonalize polynomials on present block (deferring subsequent ones)>
      orthogonalizeBlock(
              points_shifted_.col(poly_i), poly_i, block_beginning, poly_i);

      last_pt_included = pt_next;
      poly_i++;
    } else {
      //!<These points don't render good pivot value for this>
      //!<specific polynomial. Exchange some polynomials>
      //!<and try to advance moving this polynomial>
      //!<to the end of the block>

      shiftPolynomialToEndBlock(poly_i, block_end);

            //!<Comment [1]: If we are on the linear block,>
            //!<this means we won't be able to build a Fully Linear model>
    }

    tr_center_ = 0;
    points_abs_ = all_points_.leftCols(last_pt_included + 1);
    points_shifted_ = points_shifted_.leftCols(last_pt_included + 1);
    fvalues_ = all_fvalues_.head(last_pt_included + 1);

    auto cache_size = std::fmin(n_points - last_pt_included - 1, 3 * pow(dim, 2));
    modeling_polynomials_.clear();

    //!<Points not included>
    if (cache_size > 0) {
      cached_points_ = all_points_.middleCols(last_pt_included, cache_size);
      cached_fvalues_ = fvalues_.segment(last_pt_included, cache_size);
    } else {
      cached_points_.resize(0, 0);
      cached_fvalues_.resize(0);
    }
  } // end for-loop: (iter < polynomials_num)

  //!<Clean auxiliary objects>
  all_points_.resize(0, 0);
  all_fvalues_.resize(0);

  return last_pt_included < n_points; //!<model has changed>
}

bool TrustRegionModel::improveModelNfp() {

  auto rel_pivot_threshold = settings_->parameters().tr_pivot_threshold;
  auto tol_radius = settings_->parameters().tr_tol_radius;
  auto radius_factor = settings_->parameters().tr_radius_factor;

  auto radius = radius_;
  auto pivot_threshold = rel_pivot_threshold*std::fmin(1, radius);
  auto dim = points_shifted_.rows();
  auto p_ini = points_shifted_.cols();
  auto shift_center = points_abs_.col(0);

  auto tr_center = tr_center_;
  bool exit_flag = true;

  Eigen::Matrix<bool, 1, Dynamic> f_succeeded;
  f_succeeded.resize(1);
  f_succeeded.fill(false);

  int poly_i;
  double new_pivot_value;

  auto pivot_polynomials = pivot_polynomials_;
  auto polynomials_num = pivot_polynomials.size();

  if (lb_.rows() <= 0) {
    lb_.resize(dim);
    lb_.fill(-Eigen::Infinity);
  }

  if (ub_.rows() <= 0) {
    ub_.resize(dim);
    ub_.fill(Eigen::Infinity);
  }

  auto bl_shifted = lb_ - shift_center;
  auto bu_shifted = ub_ - shift_center;
  // tol_shift = 10*eps(max(1, max(abs(shift_center))));
  // Use: tol_shift = 10*std::numeric_limits<double>::epsilon();

  //TODO: test this function
  auto unshift_point = [shift_center, bl_shifted, bu_shifted](Eigen::VectorXd x) {
    Eigen::VectorXd shifted_x = x + shift_center;
    shifted_x = shifted_x.cwiseMin(bl_shifted+shift_center);
    shifted_x = shifted_x.cwiseMax(bu_shifted+shift_center);
    return shifted_x;
  };


  //!<Test if the model is already FL but old
  //!<Distance measured in inf norm>
  auto tr_center_pt = points_shifted_.col(tr_center);

  if ((tr_center_pt.lpNorm<Infinity>() > radius_factor*radius) && (p_ini > dim+1)) {
    exit_flag = false; //!<Needs to rebuild>
  } else {
    //!<The model is not old>
    if (p_ini < polynomials_num) {
      int block_end = (int)p_ini;
      int block_begining = 0;
      //!<The model is not complete>
      if (p_ini < dim+1) {
        //!<The model is not yet fully linear>
        block_begining = 1;
      } else {
        //!<We can add a point to the quadratic block>
        block_begining = (int)dim+1;
      }
      int next_position = (int)p_ini;
      double radius_used = radius;

      //!<Possibly try smaller radius>
      for (int attempts = 1; attempts<=3; attempts++) {

        //!<Iterate through available polynomials>
        for (poly_i = next_position; poly_i<=block_end; poly_i++) {

          if(!areImprovementPointsComputed()) {

              nfp_polynomial_ = orthogonalizeToOtherPolynomials(poly_i, p_ini);
              std::tie(nfp_new_points_shifted_, nfp_new_pivots_, nfp_point_found_) =
                  pointNew(nfp_polynomial_, tr_center_pt, radius_used, bl_shifted, bu_shifted, pivot_threshold);

          } else if(areImprovementPointsComputed()) {
              // Using previously computed points/pivot-points
          }

          if (nfp_point_found_) {
            for (int found_i = 0; found_i < nfp_new_points_shifted_.cols(); found_i++) {

              if(!areImprovementPointsComputed()) {

                nfp_new_point_shifted_ = nfp_new_points_shifted_.col(found_i);
                auto new_pivot_value = nfp_new_pivots_(found_i);

                nfp_new_point_abs_ = unshift_point(nfp_new_point_shifted_);
                nfp_new_fvalues_.resize(nfp_new_point_abs_.cols());

                setIsImprovementNeeded(true);

                for (int ii = 0; ii < nfp_new_point_abs_.cols(); ii++) {
                    Case *new_case = new Case(base_case_);
                    new_case->SetRealVarValues(nfp_new_point_abs_.col(ii));
                    addImprovementCase(new_case);

                    pt_case_uuid_.push_back(new_case->id());
                }

                return 5;  //TODO Probably not necessary

              } else if (areImprovementPointsComputed()) {
                  f_succeeded.resize(nfp_new_point_abs_.cols(),1);
                  f_succeeded.fill(false);
                  for (int ii = 0; ii < nfp_new_point_abs_.cols(); ii++) {

                    int nvars = (int)improvement_cases_[0]->GetRealVarVector().size();
                    auto c = improvement_cases_hash_[pt_case_uuid_[ii]];
                    nfp_new_point_abs_.col(ii) = c->GetRealVarVector();
                    nfp_new_fvalues_(ii) = c->objective_function_value();

                    std::map <string, string> stateMap = c->GetState();
                    // stateMap["EvalSt"] gives "pending" after evaluations of
                    // cases in tests...
                    // TODO Update status of cases when called from tests
                    // if(stateMap["EvalSt"] == "OKAY") {
                      f_succeeded(ii) = true;
                    // }
                  }
                  setAreImprovementPointsComputed(false);
              }
              if (f_succeeded.all()) {
                break;
              }
            }  // end-if: (found_i < new_points_shifted.cols())
            if (f_succeeded.all()) {
              break;  //!<Stop trying pivot polynomials for poly_i>
            }
          } // end-if (point_found)
          //!<Attempt another polynomial if did not break>
          if (poly_i < block_end) { // Keep the last values
              nfp_new_points_shifted_.resize(0,0);
              nfp_new_pivots_.resize(0);
              nfp_new_point_shifted_.resize(0);
              nfp_new_point_abs_.resize(0);
              nfp_point_found_ = false;
              pt_case_uuid_.clear();
          }

        } // end-for-loop: (poly_i <= block_end)

        if (nfp_point_found_ && f_succeeded.all()) {
          break; //!<(for attempts)>
        } else if (nfp_point_found_) {
          //!<Reduce radius if it did not break>
          radius_used = 0.5*radius_used;
          if (radius_used < tol_radius) {
            break;
          }
        } else {
          break;
        }
      }  // end-for-loop: (attempts<=3)

      if (nfp_point_found_ && f_succeeded.all()) {
        setIsImprovementNeeded(false);

        //!<Update this polynomial in the set>
        pivot_polynomials_[poly_i] = nfp_polynomial_;
        //!<Swap polynomials>
        std::swap(pivot_polynomials_[poly_i], pivot_polynomials_[next_position]);

        //!<Add point>
        int nr, nc;
        nr = (int)points_shifted_.rows();
        nc = (int)points_shifted_.cols();
        points_shifted_.conservativeResize(nr, nc+1);
        points_shifted_.col(nc) = nfp_new_point_shifted_;

        //!<Normalize polynomial value>
        pivot_polynomials_[next_position] = normalizePolynomial(next_position, nfp_new_point_shifted_);

        //!<Re-orthogonalize>
        pivot_polynomials_[next_position] = orthogonalizeToOtherPolynomials(next_position, p_ini);

        //!<Orthogonalize polynomials on present block (deffering subsequent ones)>
        orthogonalizeBlock(nfp_new_point_shifted_, next_position, block_begining, p_ini);

        //!<Update model and recompute polynomials>
        nr = (int)points_abs_.rows();
        nc = (int)points_abs_.cols();
        points_abs_.conservativeResize(nr, nc+1);
        points_abs_.col(nc) = nfp_new_point_abs_;

        nr = (int)fvalues_.rows();
        nc = (int)fvalues_.cols();
        fvalues_.conservativeResize(nr, nc+1);
        fvalues_.col(nc) = nfp_new_fvalues_;

        pivot_values_(next_position) = new_pivot_value;
        exit_flag = true;
      }
    } else {
      //!<The model is already complete>
      exit_flag = false;
    }
  }
  return exit_flag;
}

int TrustRegionModel::ensureImprovement() {

  bool model_complete = isComplete();
  bool model_fl = isLambdaPoised();
  bool model_old = isOld();
  int exit_flag = 4;
  bool success = false;

  if (!model_complete && (!model_old || !model_fl)) {
    //!<Calculate a new point to add>
    success = improveModelNfp(); //!<improve model>

    if (success) {
      exit_flag = 1;
    }
  } else if ((model_complete) && (!model_old)){
    //!<Replace some point with a new one that improves geometry>
    success = chooseAndReplacePoint(); //!<replace point>
    if (success) {
      exit_flag = 2;
    }
  }
  if (!success) {
    bool model_changed = rebuildModel();
    if (!model_changed) {
      if (!model_complete) {
        //!<Improve model>
        success = improveModelNfp();

      } else {
        //!<Replace point>
        success = chooseAndReplacePoint();
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
  return exit_flag;
}

bool TrustRegionModel::isLambdaPoised() {

  int dim = (int)points_abs_.rows();
  int points_num = (int)points_abs_.cols();

  double pivot_threshold = settings_->parameters().tr_pivot_threshold;
  bool result = false;

  if (settings_->parameters().tr_basis.compare("dummy")) {
    result = true;
  } else if(points_num >= dim+1) {

    //!<Fully linear, already>
    result = true;
    //!<but lets double check>
//        if (pivot_values_.lpNorm<Infinity>() > settings_->parameters().tr_pivot_threshold) {
//            Printer::ext_warn("Low pivot values.", "Optimization", "TrustRegionOptimization");
//        }
  } else {
    result = false;
  }
  return result;
}

int TrustRegionModel::changeTrCenter(
    VectorXd new_point,
    double fvalue) {

  int exit_flag = 4;
  return exit_flag;
}

std::tuple<VectorXd, double> TrustRegionModel::solveTrSubproblem() {
  double radius = radius_;
  auto obj_pol = modeling_polynomials_[0];
  auto x_tr_center = points_abs_.col(tr_center_);
  obj_pol = shiftPolynomial(obj_pol, -x_tr_center); //!<Shift to origin>

  Eigen::VectorXd trial_point(getXDim(), 1);
  double trial_fval;
  int exit_flag;


  std::tie(trial_point, trial_fval, exit_flag) = minimizeTr(obj_pol, x_tr_center, radius,lb_, ub_);
  double current_fval = fvalues_(tr_center_);
  double trial_decrease = current_fval - trial_fval;

  if (current_fval <= trial_fval) {
    Printer::ext_info("current_fval <= trial_fval", "Optimization", "TrustRegionModel");
  }

  //  double tol_interp = max(1e-08, eps(max(double(1), fvalues_.maxCoeff()))*1e03);
  double tol_interp = 1e-08;
  int n_points = points_abs_.cols();

  for (int k=0; k<n_points; k++) {
    auto val = evaluatePolynomial(obj_pol, points_abs_.col(k));
    auto error_interp = abs(val - fvalues_(k));
    if (error_interp > tol_interp) {
      Printer::ext_info("error_interp <= tol_interp", "Optimization", "TrustRegionModel");
    }
  }
  return make_tuple(trial_point, trial_decrease);
}

int TrustRegionModel::tryToAddPoint(VectorXd new_point, double new_fvalue) {
  int exit_flag = 0;
  bool point_added = false;
  double relative_pivot_threshold = settings_->parameters().tr_pivot_threshold;

  if (!isComplete()) {
    //!<Add this point>
    exit_flag = addPoint(new_point, new_fvalue, relative_pivot_threshold);
    if (exit_flag >0){
      point_added = true;
    }
  }

  if (!point_added) {
    //!<Save information about this new point, just not to lose>
    int nc = cached_points_.cols();
    cached_points_.resize(cached_points_.rows(), nc+1);
    cached_points_.col(nc) = new_point;

    int nc_f = cached_fvalues_.cols();
    cached_fvalues_.resize(nc_f+1);
    cached_fvalues_(nc_f) = new_fvalue;

    //!<Either add a geometry improving point or rebuild model>
    exit_flag = ensureImprovement();

  } else {
    exit_flag = 1;
  }
  return exit_flag;
}

int TrustRegionModel::addPoint(VectorXd new_point, double new_fvalue, double relative_pivot_threshold) {

  int exit_flag;
  double pivot_threshold = min(double(1), radius_)*relative_pivot_threshold;
  int dim = points_abs_.rows();
  int last_p = points_abs_.cols()-1;

  VectorXd new_point_shifted(dim);
  VectorXd shift_center(dim);
  int pivot_value = 0;
  bool success = false;

  int polynomials_num = pivot_polynomials_.size();

  if (last_p >= 0) {
    shift_center = points_abs_.col(0);
    new_point_shifted = new_point - shift_center;

    int nc = points_shifted_.cols();
    points_shifted_.conservativeResize(points_shifted_.rows(), nc+1);
    points_shifted_.col(nc) = new_point_shifted;
  } else {
    Printer::ext_info("cmg:no_point",
                      "Tr model had no point. Should this ever happen ?", "addPoint");
  }

  int next_position = last_p+1;
  if (next_position == 0) {
    //!<Should be rebuilding model!>
    tr_center_ = 0;
    //pivot_polynomials_ = bandPrioritizingBasis(dim); //TODO: implement a band prioritizing basis
    exit_flag = 1;
    Printer::ext_info("cmg:no_point",
                      "TR model had no point. This should never happen", "addPoint");
  } else if (next_position < polynomials_num) {
    int block_beginning =0;
    int block_end =0;
    if (next_position <= dim)  {
      //!<Add to linear block>
      block_beginning =1;
      block_end = dim;
    } else {
      //!<Add to quadratic block>
      block_beginning = dim+1;
      block_end = polynomials_num-1;
    }

    tie(pivot_value, success) = choosePivotPolynomial(next_position, block_end, pivot_threshold);

    if (success) {
      //!<Normalize polynomial value>
      pivot_polynomials_[next_position] = normalizePolynomial(next_position, new_point_shifted);

      //!<Re-orthogonalize>
      pivot_polynomials_[next_position] = orthogonalizeToOtherPolynomials(next_position, next_position-1);

      //!<Orthogonalize polynomials on present block (deferring subsequent ones)>
      orthogonalizeBlock(new_point_shifted, next_position, block_beginning, next_position-1);

      exit_flag = 1;
    } else {
      exit_flag = 0;
    }
  } else {
    //!<Model is full. Should remove another point to add this one>
    exit_flag = 0;
  }
  if (exit_flag >0) {
    points_shifted_.col(next_position) = new_point_shifted;
    fvalues_.resize(fvalues_.cols()+1);
    fvalues_(next_position) = new_fvalue;

    int nr = (int)points_abs_.rows();
    int nc = (int)points_abs_.cols();
    points_abs_.conservativeResize(nr, nc+1);
    points_abs_.col(next_position) = new_point;

    modeling_polynomials_.clear();
    pivot_values_.resize(pivot_values_.cols()+1);
    pivot_values_(next_position) = pivot_value;
  }
  return exit_flag;
}

tuple<double, bool> TrustRegionModel::choosePivotPolynomial(int initial_i, int final_i, double tol) {
  int last_point = initial_i - 1;
  auto incumbent_point = points_shifted_.col(initial_i);
  bool success = false;
  double pivot_value = 0;

  for (int k=initial_i; k<=final_i; k++) {
    auto polynomial = orthogonalizeToOtherPolynomials(k, last_point);
    auto val = evaluatePolynomial(polynomial, incumbent_point);
    if (abs(val) > tol) {
      //!<Accept polynomial>
      success = true;

      //!<Swap polynomials>
      pivot_polynomials_[k] = pivot_polynomials_[initial_i];
      pivot_polynomials_[initial_i] = polynomial;
      pivot_value = val;
      break;
    } else {
      //!<We won't even save the orthogonalization>
      success = false;
    }
  }
  return make_tuple(pivot_value, success);
}


void TrustRegionModel::computePolynomialModels() {

  int dim = (int)points_abs_.rows();
  int points_num = (int)points_abs_.cols();
  int functions_num = 1; //!<Code currently supports only 1 function>
  int linear_terms = dim+1;
  int full_q_terms = (dim+1)*(dim+2)/2;
  std::vector<Polynomial> polynomials(functions_num);

  if ((linear_terms < points_num) && (points_num < full_q_terms)) {
    //!<Compute quadratic model>
    polynomials = computeQuadraticMNPolynomials();
  }

  //!<Ideally we should check if the model is a badly conditioned system
  if ((points_num <= linear_terms) || (points_num == full_q_terms)) {
    //!<Compute model with incomplete (complete) basis>
    auto l_alpha = nfpFiniteDifferences(points_num);
    for (int k=functions_num-1; k>=0; k--) {
      polynomials[k] = combinePolynomials(points_num, l_alpha);
      polynomials[k] = shiftPolynomial(polynomials[k], points_shifted_.col(tr_center_));
    }
  }
  modeling_polynomials_ = polynomials;
}

void TrustRegionModel::shiftPolynomialToEndBlock(
    int point_index,
    int block_end) {

  std::swap(pivot_polynomials_[point_index], pivot_polynomials_[block_end]);
  for (int pi=block_end-1; pi>point_index; pi--) {
        std::swap(pivot_polynomials_[pi], pivot_polynomials_[point_index]);
  }
}

void TrustRegionModel::sortVectorByIndex(
    VectorXd &vec,
    const VectorXd &ind) {
  VectorXd vec_ord(vec.size());
  for (int i=0; i < vec.size(); i++) {
    int index = int(ind(i));
    vec_ord(i) = vec(index);
  }

  for(int i=0; i<vec.size(); i++) {
    vec(i) = vec_ord(i);
  }
  vec_ord.resize(0);
}

void TrustRegionModel::sortVectorByIndex(
    RowVectorXd &vec,
    const VectorXd &ind) {
  VectorXd vec_ord(vec.size());
  for (int i=0; i < vec.size(); i++) {
    int index = int(ind(i));
    vec_ord(i) = vec(index);
  }

  for(int i=0; i<vec.size(); i++) {
    vec(i) = vec_ord(i);
  }
  vec_ord.resize(0);
}

void TrustRegionModel::sortMatrixByIndex(
    Matrix<double,Dynamic,Dynamic> &points,
    const VectorXd &ind) {

  Matrix<double,Dynamic,Dynamic> points_ord(points.rows(), points.cols());
  for (int i=0; i<points.cols(); i++) {
    int index = int(ind(i));
    points_ord.col(i) << points.col(index);
  }

  for (int i=0; i<points.cols(); i++) {
    points.col(i) << points_ord.col(i);
  }
  points_ord.resize(0,0);
}

void TrustRegionModel::nfpBasis(int dim) {
  //!<Number of terms>
  int poly_num = (dim+1)*(dim+2)/2;
  int linear_size = dim+1;

  //!<Calculating basis of polynomials>
  pivot_polynomials_.resize(poly_num);
  pivot_polynomials_[poly_num-1].dimension = dim;
  pivot_polynomials_[poly_num-1].coefficients.resize(poly_num);
  pivot_polynomials_[poly_num-1].coefficients.setZero(poly_num);

  for (int i=0; i<linear_size;i++) {
    pivot_polynomials_[i].dimension = dim;
    pivot_polynomials_[i].coefficients.resize(poly_num);
    pivot_polynomials_[i].coefficients.setZero(poly_num);
    pivot_polynomials_[i].coefficients(i) = 1;
  }

  //!<Quadratic entries>
  int c0 = 0;
  int m = 0;
  int n = 0;
  VectorXd g0(dim);
  g0.setZero(dim);

  for (int poly_i=linear_size; poly_i<poly_num; poly_i++) {
    MatrixXd H(dim,dim);
    H.setZero(dim,dim);
    if (m == n) {
      H(m,n) = 2;
    } else {
      H(m,n) = 1;
      H(n,m) = 1;
    }

    pivot_polynomials_[poly_i] = matricesToPolynomial(c0, g0, H);
    if (n < dim-1) {
      n++;
    } else {
      m++;
      n = m;
    }
  }
}

Polynomial TrustRegionModel::matricesToPolynomial(
    double c0,
    const VectorXd &g0,
    const MatrixXd &H) {

  int dim = (int)g0.size();
  int n_terms = (dim+1)*(dim+2)/2;
  VectorXd coefficients(n_terms);
  coefficients.setZero(n_terms);

  //!<Zero order>
  coefficients(0) = c0;

  //!<First order>
  int ind_coefficients = dim;
  coefficients.segment(1,ind_coefficients) = g0;

  //!<Second order>
  for (int k=0; k<dim; k++) {
    for (int m=0; m <= k; m++) {
      ind_coefficients = ind_coefficients + 1;
      coefficients(ind_coefficients) = H(k, m);

      if (H(m, k) != H(k, m)) {
        Printer::ext_info("H not symmetrical: ",
                          "Optimization", "Trust Region");
      }
    }
  }

  Polynomial p;
  p.dimension = dim;
  p.coefficients = coefficients;
  return p;
}

std::tuple<double, VectorXd, MatrixXd> TrustRegionModel::coefficientsToMatrices(
    int dimension,
    VectorXd coefficients) {
  if (coefficients.size() == 0) {
    coefficients.resize(dimension);
    coefficients.setZero(dimension);
  }

  int n_terms = (dimension+1)*(dimension+2)/2;
  if (coefficients.size() != n_terms) {
    Printer::ext_warn("Wrong dimension of coefficients.",
                      "Optimization", "TrustRegionModel");
    throw std::runtime_error(
        "Failed to convert polynomial coefficients to matrices.");
  }

  //!< Constant term>
  double c = coefficients(0);

  //!< Order one term>
  int idx_coefficients = dimension + 1;
  auto g = coefficients.segment(1,idx_coefficients-1);

  //!< Second order term>
  MatrixXd H(dimension, dimension);
  H.setZero(dimension,dimension);

  for (int k=0; k<dimension; k++) {
    for (int m=0; m < k; m++) {
      idx_coefficients++;
      H(k,m) = coefficients(idx_coefficients);
      H(m,k) = H(k,m);
    }
  }

  return std::make_tuple(c, g, H);
}


Polynomial TrustRegionModel::normalizePolynomial(
    int poly_i,
    VectorXd point) {
  auto polynomial = pivot_polynomials_[poly_i];
  auto val = evaluatePolynomial(polynomial, point);
  for (int i=0; i<3; i++) {
    polynomial = multiplyPolynomial(polynomial, (double) 1/val);
    val = evaluatePolynomial(polynomial, point);
    if ((val - 1) == 0) {
      break;
    }
  }
  return polynomial;
}

Polynomial TrustRegionModel::orthogonalizeToOtherPolynomials(
    int poly_i,
    int last_pt) {
  auto polynomial = pivot_polynomials_[poly_i];
  for (int n=0; n<=last_pt; n++) {
    if (n != poly_i) {
      polynomial = zeroAtPoint(polynomial,
                               pivot_polynomials_[n],
                               points_shifted_.col(n));
    }
  }
  return polynomial;
}

void TrustRegionModel::orthogonalizeBlock(
    VectorXd point,
    int np,
    int block_beginning,
    int block_end) {
  for (int p = block_beginning; p<=block_end; p++) {
    if (p != np) {
      pivot_polynomials_[p] = zeroAtPoint(pivot_polynomials_[p], pivot_polynomials_[np], point);
    }
  }
}


Polynomial TrustRegionModel::zeroAtPoint(
    Polynomial p1,
    Polynomial p2,
    VectorXd x) {
  Polynomial p;
  p = p1;

  auto px = evaluatePolynomial(p, x);
  auto p2x = evaluatePolynomial(p2, x);
  int iter = 1;
  while (px != 0) {
    p = addPolynomial(p, multiplyPolynomial(p2, -px/p2x));
    px = evaluatePolynomial(p, x);

    if (iter >= 2)
      break;
    iter++;
  }
  return p;
}

double TrustRegionModel::evaluatePolynomial(
    Polynomial p1,
    VectorXd x) {
  double c;
  double terms[3];
  VectorXd g(p1.dimension);
  MatrixXd H(p1.dimension, p1.dimension);
  std::tie(c, g, H) = coefficientsToMatrices(p1.dimension, p1.coefficients);

  terms[0] = c;
  terms[1] = g.transpose()*x;
  terms[2] = (x.transpose()*H*x);
  terms[2] *= 0.5;

  return terms[0] + terms[1] + terms[2];

}

Polynomial TrustRegionModel::addPolynomial(
    Polynomial p1,
    Polynomial p2) {
  if (p1.dimension != p2.dimension) {
    Printer::ext_warn(
        "Summation of polynomials "
        "with different dimensions.",
        "Optimization", "TrustRegionModel");
    throw std::runtime_error("Failed to compute add polynomials. "
                             "They have different dimensions.");
  }

  Polynomial p;
  p.dimension = p1.dimension;
  p.coefficients.resize(p1.coefficients.size());
  p.coefficients = p1.coefficients + p2.coefficients;

  if (p.coefficients.size() == 0) {
    Printer::ext_warn("Empty resulting polynomial.",
                      "Optimization", "TrustRegionModel");
    throw std::runtime_error("Failed to compute add polynomials. "
                             "The resulting polynomial is empty.");
  }

  return p;
}

Polynomial TrustRegionModel::multiplyPolynomial(
    Polynomial p1,
    double factor) {
  Polynomial p = p1;
  p.coefficients = p1.coefficients.array()*factor;
  return p;
}

int TrustRegionModel::findBestPoint() {

  int dim = (int)points_abs_.rows();
  int n_points = (int)points_abs_.cols();
  int best_i = 0;
  auto min_f = std::numeric_limits<double>::infinity();


  for (int k=0; k<n_points; k++) {
    VectorXd point = points_abs_.col(k);
    if (((point - lb_).minCoeff() >= 0) &&  ((ub_ - point).minCoeff() > 0)) {
      auto val = fvalues_(k);
      if (val < min_f) {
        min_f = val;
        best_i = k;
      }
    }
  }
  return best_i;
}

std::vector<Polynomial>
TrustRegionModel::computeQuadraticMNPolynomials() {

  int dim = (int)points_abs_.rows();
  int points_num = (int)points_abs_.cols();
  int functions_num = (int)fvalues_.rows();
  std::vector<Polynomial> polynomials(functions_num);

  MatrixXd points_shifted = MatrixXd::Zero(dim, points_num-1);
  RowVectorXd fvalues_diff = RowVectorXd::Zero(functions_num, points_num-1);

  int m2 = 0;
  for (int m = tr_center_; (m<points_num) && (m != tr_center_); m++) {
    points_shifted.col(m2) = points_abs_.col(m) - points_abs_.col(tr_center_);
    fvalues_diff(m2) = fvalues_(m) - fvalues_(tr_center_);
    m2++;
  }

  MatrixXd M = points_shifted.transpose()*points_shifted;
  M += 0.5*(MatrixXd)(M.array().square());

  //!<Solve symmetric system>
  //!< //TODO: choose best algorithm for solving
  //!< linear system of equations automatically.

  // Works w/ 3rd-party-Eigen BUT crashes with system-Eigen
  // PartialPivLU< Ref<MatrixXd> > lu(M);

  // Works w/ system-Eigen AND 3rd-party-Eigen
  PartialPivLU<MatrixXd> lu(M);

  lu.compute(M); //!<Update LU matrix>
  auto mult_mn = lu.solve(fvalues_diff.transpose());

  //TODO: raise a warning if the system is badly conditioned
  // using the resulting conditioning number (tol=1e4*eps(1))

  for (int n=0; n<functions_num; n++) {
    VectorXd g = VectorXd::Zero(dim);
    MatrixXd H = MatrixXd::Zero(dim,dim);
    for (int m=0; m<points_num-1; m++) {
      g += mult_mn(m)*points_shifted.col(m);
      H += mult_mn(m)*(points_shifted.col(m)*points_shifted.col(m).transpose());
    }
    auto c = fvalues_(tr_center_);
    polynomials[n] = matricesToPolynomial(c, g, H);
  }
  return polynomials;
}

RowVectorXd TrustRegionModel::nfpFiniteDifferences(int points_num) {
  //!<Change so we can interpolate more functions at the same time>
  int dim = (int)points_shifted_.cols();
  RowVectorXd l_alpha = fvalues_;
    std::vector<Polynomial> polynomials = std::vector<Polynomial>(pivot_polynomials_.begin(), pivot_polynomials_.begin() + points_num);

  //!<Remove constant polynomial>
  for (int m=1; m<points_num; m++) {
    auto val = evaluatePolynomial(polynomials[0], points_shifted_.col(m));
    l_alpha(m) = l_alpha(m) - l_alpha(0)*val;
  }

  //!<Remove terms corresponding to degree 1 polynomials>
  for (int m=dim+1; m <points_num; m++) {
    for (int n=1; n<dim+1; n++) {
      auto val = evaluatePolynomial(polynomials[n], points_shifted_.col(m));
      l_alpha(m) = l_alpha(m) - l_alpha(n)*val;
    }
  }
  polynomials.clear();
  return l_alpha;
}

Polynomial TrustRegionModel::combinePolynomials(
    int points_num,
    RowVectorXd coefficients) {

  auto polynomials = std::vector<Polynomial>(
      pivot_polynomials_.begin(),
      pivot_polynomials_.begin() + points_num);

  int terms = (int)polynomials.size();

  if ((terms == 0) || (coefficients.size() != terms)) {
    Printer::ext_warn("Polynomial and coefficients have different sizes.",
                      "Optimization", "TrustRegionModel");
    throw std::runtime_error(
        "Failed to combine polynomials. Polynomial "
        "and coefficients have different dimensions.");
  }

  auto p = multiplyPolynomial(polynomials[0], coefficients(0));
  for (int k = 1; k < terms; k++) {
    p = addPolynomial(p, multiplyPolynomial(polynomials[k], coefficients[k]));
  }
  return p;
}

Polynomial TrustRegionModel::shiftPolynomial(Polynomial polynomial, VectorXd s) {

  double terms[3];
  double c;
  VectorXd g(polynomial.dimension);
  MatrixXd H(polynomial.dimension, polynomial.dimension);

  std::tie(c, g, H) = coefficientsToMatrices(polynomial.dimension, polynomial.coefficients);

  double c_mod = c + (double)(g.transpose()*s) + 0.5*(double)(s.transpose()*H*s);

  VectorXd g_mod(g.size());
  g_mod = g + H*s;

  return matricesToPolynomial(c_mod,g_mod, H);
}

bool TrustRegionModel::isComplete() {

  int dim = (int)points_abs_.rows();
  int points_num = (int)points_abs_.cols();
  int max_terms = ((dim+1)*(dim+2))/2;
  if (points_num > max_terms) {
    Printer::ext_warn("Too many points in the Trust Region model.",
                      "Optimization", "TrustRegionModel");
  }
  return (points_num >= max_terms );;
}

bool TrustRegionModel::isOld() {
  VectorXd distance(points_abs_.rows());
  distance = points_abs_.col(0) - points_abs_.col(tr_center_);
  return (distance.lpNorm<Infinity>() > settings_->parameters().tr_radius_factor);
}

bool TrustRegionModel::chooseAndReplacePoint() {
  double pivot_threshold = settings_->parameters().tr_exchange_threshold;
  double radius = radius_;

  int dim = points_shifted_.rows();
  int points_num = points_shifted_.cols();
  auto points_shifted = points_shifted_;

  int tr_center = tr_center_;
  auto shift_center = points_abs_.col(0);
  auto tr_center_x = points_shifted_.col(tr_center);

  auto pivot_values = pivot_values_;
  auto pivot_polynomials = pivot_polynomials_;

  int linear_terms = dim+1;

  auto bl_shifted = lb_ - shift_center;
  auto bu_shifted = ub_ - shift_center;
  bool success = false;

  Eigen::Matrix<bool, 1, Dynamic> f_succeeded;
  f_succeeded.resize(1);
  f_succeeded.fill(false);

//  unshift_point = @(x) max(min(x + shift_center, bu), bl);
//  tol_shift = 10*eps(max(1, norm(shift_center, inf)));

  //TODO: test this function
  auto unshift_point = [shift_center, bl_shifted, bu_shifted](Eigen::VectorXd x) {
    Eigen::VectorXd shifted_x = x + shift_center;
    shifted_x = shifted_x.cwiseMin(bl_shifted+shift_center);
    shifted_x = shifted_x.cwiseMax(bu_shifted+shift_center);
    return shifted_x;
  };

  //!<Reorder points based on their pivot_values>
  piv_order_.setLinSpaced(pivot_values_.size(), 0, pivot_values_.size() - 1);
  std::sort(piv_order_.data(),
            piv_order_.data() + piv_order_.size(),
            std::bind(compare,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      pivot_values_.data()));

  int polynomials_num = pivot_polynomials_.size();
  //!<Could iterate through all pivots, but will try just dealing with the worst one>
  auto pos = piv_order_(0);
  if ((pos == 0) || (pos == tr_center) || (pos <= linear_terms  && points_num > linear_terms)) {
    //!<Better to just rebuild model>
    success = false;
  } else {
    double new_pivot_value = 1;
    auto current_pivot_value = pivot_values_(pos);

    if (!areReplacementPointsComputed()) {
      std::tie(repl_new_point_shifted_, repl_new_pivots_, repl_point_found_) = pointNew(pivot_polynomials_[pos], tr_center_x, radius_, bl_shifted, bu_shifted, pivot_threshold);
    }

    if (repl_point_found_) {
      for (int found_i = 0; found_i < repl_new_point_shifted_.cols(); found_i++) {

        if (!areReplacementPointsComputed()) {

          repl_new_point_shifted_ = repl_new_point_shifted_.col(found_i);
          auto new_pivot_value = repl_new_pivots_(found_i);

          repl_new_point_abs_ = unshift_point(repl_new_point_shifted_);
          repl_new_fvalues_.resize(repl_new_point_abs_.cols());

          setIsReplacementNeeded(true);

          for (int ii = 0; ii < repl_new_point_abs_.cols(); ii++) {
            Case *new_case = new Case(base_case_);
            new_case->SetRealVarValues(repl_new_point_abs_.col(ii));
            addReplacementCase(new_case);

            repl_pt_case_uuid_.push_back(new_case->id());
          }
          return 5; //TODO Probably not necessary

        } else if (areReplacementPointsComputed()) {
          f_succeeded.resize(nfp_new_point_abs_.cols(), 1);
          f_succeeded.fill(false);

          for (int ii = 0; ii < repl_new_point_abs_.cols(); ii++) {

            int nvars = (int) replacement_cases_[0]->GetRealVarVector().size();
            auto c = replacement_cases_hash_[repl_pt_case_uuid_[ii]];
            repl_new_point_abs_.col(ii) = c->GetRealVarVector();
            repl_new_fvalues_(ii) = c->objective_function_value();

            std::map<string, string> stateMap = c->GetState();
            f_succeeded(ii) = true;

          }
          setAreReplacementPointsComputed(false);
        }
        if (f_succeeded.all()) {
          break;
        }
      }
    }

    if (repl_point_found_ && f_succeeded.all()) {
      setIsReplacementNeeded(false);

      //!<Normalize polynomial value>
      pivot_polynomials_[pos] = normalizePolynomial(pos, repl_new_point_shifted_);

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
      pivot_polynomials_[pos] = orthogonalizeToOtherPolynomials(pos, block_end);

      //!<Orthogonalize block>
      orthogonalizeBlock(repl_new_point_shifted_, pos, block_beginning, block_end);

      //!<Update model and recompute polynomials>
      cached_points_ = points_abs_.leftCols(points_abs_.rows());
      cached_fvalues_ = fvalues_.head(fvalues_.size());

      int nr = (int)points_abs_.rows();
      int nc = (int)points_abs_.cols();

      points_abs_.conservativeResize(nr, nc+1);
      points_abs_.col(nc) = repl_new_point_abs_;

      nr = (int)fvalues_.rows();
      nc = (int)fvalues_.cols();
      fvalues_.conservativeResize(nr, nc+1);
      fvalues_.col(nc) = repl_new_fvalues_;

      pivot_values_((int)pos) *= new_pivot_value;

      modeling_polynomials_.clear();
      success = true;
    }
  }
}


tuple<Eigen::MatrixXd, Eigen::RowVectorXd, bool>
TrustRegionModel::pointNew(Polynomial polynomial,
                           Eigen::VectorXd tr_center_point,
                           double radius_used,
                           Eigen::VectorXd bl,
                           Eigen::VectorXd bu,
                           double pivot_threshold) {

    cout << "Find new pt -> Calling minimize_tr() twice" << endl;

    Eigen::VectorXd new_point_min(getXDim(), 1);
    Eigen::VectorXd new_point_max(getXDim(), 1);
    Eigen::MatrixXd new_points(getXDim(), 2);

    double pivot_min, pivot_max;
    Eigen::VectorXd new_pivot_values(new_points.cols(), 1);

    int exitflag_min, exitflag_max;
    bool point_found;

    tie(new_point_min, pivot_min, exitflag_min) = minimizeTr(polynomial,
                                                             tr_center_point,
                                                             radius_used,
                                                             bl, bu);

    auto polynomial_max = multiplyPolynomial(polynomial, -1);

    tie(new_point_max, pivot_max, exitflag_max) = minimizeTr(polynomial_max,
                                                             tr_center_point,
                                                             radius_used,
                                                             bl, bu);

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
            new_points.resize(dim_,1);
            new_points.col(0) = new_point_min;

            new_pivot_values.resize(new_points.cols(), 1);
            new_pivot_values << pivot_min;
        }
        point_found = true;

    } else if ((exitflag_max >= 0)
               && (abs(pivot_max) >= pivot_threshold)) {
        new_points.resize(dim_,1);
        new_points.col(0) = new_point_max;

        new_pivot_values.resize(new_points.cols(), 1);
        new_pivot_values << pivot_max;
        point_found = true;

    } else {
        point_found = false;
        new_points.resize(dim_, 0);

        new_pivot_values.resize(new_points.cols(), 1);
        new_pivot_values << 0;
    }

    cout << "Found " + to_string(new_points.cols()) +
            " new pts/pivot-pts" << endl;

    return make_tuple(new_points, new_pivot_values, point_found);
}

tuple<Eigen::VectorXd, double, int>
TrustRegionModel::minimizeTr(Polynomial polynomial,
                             Eigen::VectorXd x_tr_center,
                             double radius,
                             Eigen::VectorXd bl,
                             Eigen::VectorXd bu) {

    bool dbg = false;

    Eigen::VectorXd x(dim_, 1);
    Eigen::VectorXd bl_tr(dim_, 1), bu_tr(dim_, 1);
    Eigen::VectorXd bl_mod(dim_, 1), bu_mod(dim_, 1);

    double fval;
    double tol_norm, tol_arg, tol_tr, tol_eps;
    int exitflag = -1;

    // TR tol
    // CG: tol_arg = max(1, norm(x_tr_center, inf));
    if (x_tr_center.lpNorm<Infinity>() > 1.0) {
        Printer::ext_warn("X > 1 in EPS(X).",
                          "minimizeTr",
                          "TrustRegionModel");
    }
    tol_tr = 10*std::numeric_limits<double>::epsilon();

    if (dbg) { // Can be deleted
        cout << "Debug input polyn/x_tr_center/radius" << endl;
        auto pc = polynomial.coefficients;
        cout << "polynomial.coefficients = \n" << pc << endl;
        cout << "x_tr_center = \n" << x_tr_center << endl;
        cout << "radius = \n" << radius << endl ;
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

    if (dbg) { // Can be deleted
        cout << "Debug bounds_mod/tol_tr" << endl;
        cout << "bl_mod = \n" << bl_mod << endl;
        cout << "bu_mod = \n" << bu_mod << endl;
        cout << "tol_tr = \n" << tol_tr << endl;
    }

    double c;
    Eigen::VectorXd g(polynomial.dimension);
    Eigen::MatrixXd H(polynomial.dimension,
                      polynomial.dimension);

    tie(c, g, H) =
            coefficientsToMatrices(polynomial.dimension,
                                   polynomial.coefficients);

    if (dbg) { // Can be deleted
        cout << "Debug polynomial matrices c, g, H:" << endl;
        cout << "c = \n" << c << endl ;
        cout << "g = \n" << g << endl ;
        cout << "H = \n" << H << endl;
    }

    cout << "Minimize TR using SNOPT";

    // Set prob name (incl. spec) to be solved
    Optimizer::EmbeddedProblem prob;
    prob.setProbName("TRMod_A");
    // prob.setProbName("Rosenbrock"); // dbg

    // Set prob dims
    prob.setNunVars(getXDim());
    prob.setNumLinConst(0);
    prob.setNunNnlConst(0);

    // Set init point, bounds on vars
    prob.setXInit(x_tr_center);
    prob.setXUb(bu_mod);
    prob.setXLb(bl_mod);

    VectorXd FUb(1 + prob.getNunNnlConst());
    VectorXd FLb(1 + prob.getNunNnlConst());

    // Bounds on objective of tr prob
    FUb(0) = 1e20;
    FLb(0) = -1e20;

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

    if (prob.getSNOPTExitCode() == 1) {
        exitflag = 1;
    }

    if (false) {
        cout << "Debug SNOPT xsol/fsol/exitflag:";
        cout << "xsol = \n" << prob.getXSol() << endl;
        cout << "fsol = \n" << prob.getFSol() << endl;
        cout << "exit = \n" << exitflag << endl ;
    }

    return make_tuple(x, fval, exitflag);
}

}  // namespace Optimizers
}  // namespace Optimization
