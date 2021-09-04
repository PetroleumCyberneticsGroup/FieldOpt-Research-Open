/***********************************************************
Created by bellout on 21.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
TRFrame built from TrustRegionModel.cpp
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

#include "TRFrame.h"

namespace Optimization {
namespace Optimizers {

bool comp(const int &lhs, const int &rhs,
             const double *distances) {
  return distances[lhs] < distances[rhs];
}

bool compAbs(const int &lhs, const int &rhs,
                const double *distances) {
  return abs(distances[lhs]) < abs(distances[rhs]);
}

using OptMode = Settings::Optimizer::OptimizerMode;

// _________________________________________________________
// TRMOD CONSTRUCTOR
TRFrame::TRFrame(VectorXd &lb, VectorXd &ub,
                 Case *base_case, Settings::Optimizer *settings) {
  settings_ = settings;
  setSettings(settings);
  dbg_ = new TRDebug(tr_prob_name_, this);
  base_case_ = base_case;
  SNOPTSolver_ = new SNOPTSolver(); // subproblem solver

  base_case_->set_objf_value(fmult_ * base_case_->objf_value());
  setXDim(base_case_->GetRealVarVector().size());

  lb_ = lb;
  ub_ = ub;
  tr_center_ = 0; // index
  // cache_max_ = (int)3*pow(dim_,2); // not used

  ti_ = current_time();
  t0_ = ti_;
}

// ---------------------------------------------------------
// TR 1BUILD
#include "TRBuild.cpp"

// ---------------------------------------------------------
// TR 2SOLVE
#include "TRSolve.cpp"

// ---------------------------------------------------------
// TR 3CHECK
#include "TRCheck.cpp"
// TRFrame::findBestPt
// TRFrame::getCritGrad
// TRFrame::setSettings

// ---------------------------------------------------------
// TR 4POLY
#include "TRPolys.cpp"
// TRFrame::choosePivPoly
// TRFrame::normalizePoly
// TRFrame::orthogzToOthrPolys
// TRFrame::orthogzBlock
// TRFrame::zeroAtPt
// TRFrame::evaluatePoly
// TRFrame::addPoly
// TRFrame::multiplyPoly
// TRFrame::combinePolys
// TRFrame::shiftPoly
// TRFrame::coeffsToMatrices
// TRFrame::matricesToPoly

// _________________________________________________________
// SETSETTINGS
void TRFrame::setSettings(Settings::Optimizer *settings) {
  tr_init_rad_   = settings_->parameters().tr_init_rad;
  tr_tol_f_      = settings_->parameters().tr_tol_f;
  tr_eps_c_      = settings_->parameters().tr_eps_c;
  tr_eta_0_      = settings_->parameters().tr_eta_0;
  tr_eta_1_      = settings_->parameters().tr_eta_1;
  tr_piv_thld_   = settings_->parameters().tr_piv_thld;
  tr_add_thld_   = settings_->parameters().tr_add_thld;
  tr_exch_thld_  = settings_->parameters().tr_xch_thld;
  tr_rad_max_    = settings_->parameters().tr_rad_max;
  tr_rad_fac_    = settings_->parameters().tr_rad_fac;
  tr_tol_rad_    = settings_->parameters().tr_rad_tol;
  tr_gamma_inc_  = settings_->parameters().tr_gamma_inc;
  tr_gamma_dec_  = settings_->parameters().tr_gamma_dec;
  tr_crit_mu_    = settings_->parameters().tr_crit_mu;
  tr_crit_omega_ = settings_->parameters().tr_crit_omega;
  tr_crit_beta_  = settings_->parameters().tr_crit_beta;
  tr_lower_bnd_  = settings_->parameters().tr_lower_bnd;
  tr_upper_bnd_  = settings_->parameters().tr_upper_bnd;
  tr_iter_max_   = settings_->parameters().tr_iter_max;
  tr_num_init_x_ = settings_->parameters().tr_num_init_x;
  tr_basis_      = settings_->parameters().tr_basis;
  tr_init_smpln_ = settings_->parameters().tr_init_smpln;
  tr_rng_seed_   = settings_->parameters().rng_seed;
  tr_prob_name_  = settings_->parameters().tr_prob_name;
  vp_ = settings_->verbParams();

  fmult_ = 1.0;
  if (settings_->mode() == OptMode::Maximize) {
    fmult_ = -1.0;
  }

  radius_ = tr_init_rad_;
  pivot_values_.conservativeResize(0);
  pivot_order_.conservativeResize(0);
  cached_fvals_.conservativeResize(0);
  cached_pts_.conservativeResize(0,0);
  pts_shftd_.conservativeResize(0,0);
  pts_shftd_temp_.conservativeResize(0,0);
  index_vector_.conservativeResize(0);
  distances_.conservativeResize(0);
}

// _________________________________________________________
// SUBMITTEMPINITCASES
void TRFrame::submitTempInitCases() {
  init_cases_ = init_cases_temp_;

  int nvars = (int)init_cases_[0]->GetRealVarVector().size();
  pts_abs_.setZero(nvars, init_cases_.size());
  fvals_.setZero(init_cases_.size());

  int ii = 0;
  for (Case *c : init_cases_) {
    pts_abs_.col(ii) = c->GetRealVarVector();

    if (settings_->mode() == Settings::Optimizer::OptimizerMode::Maximize) {
      fvals_(ii) = -c->objf_value();
    } else {
      fvals_(ii) = c->objf_value();
    }

    ii++;
  }
}

// _________________________________________________________
// SUBMITTEMPIMPRCASES
void TRFrame::submitTempImprCases() {
  impr_cases_ = impr_cases_temp_;
  for (Case *c : impr_cases_) {
    impr_cases_hash_.insert(c->id(), c);
  }
}

// _________________________________________________________
// SUBMITTEMPREPLCASES
void TRFrame::submitTempReplCases() {
  repl_cases_ = repl_cases_temp_;
  for (Case *c : repl_cases_) {
    repl_cases_hash_.insert(c->id(), c);
  }
}

// _________________________________________________________
// CLEARIMPRCASESLIST
void TRFrame::clearImprCasesList() {
  impr_cases_.clear();
  impr_cases_temp_.clear();
}

// _________________________________________________________
// CLEARREPLCASESLIST
void TRFrame::clearReplCasesList() {
  repl_cases_.clear();
  repl_cases_temp_.clear();
}

}
}