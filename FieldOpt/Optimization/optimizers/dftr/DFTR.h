/***********************************************************
Created by bellout on 21.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
DFTR built from TrustRegionOptimization.h
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

#ifndef FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_DFTR_H_
#define FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_DFTR_H_

#include "Optimization/optimizer.h"
#include "Optimization/optimizers/dftr/TRFrame.h"
#include "Optimization/optimizers/dftr/TRDebug.h"

#include <Eigen/Core>
#include <iostream>
#include <string>

namespace Optimization {
namespace Optimizers {

using namespace Eigen;
using Model::Properties::VarPropContainer;
using Optimization::Constraints::ConstraintHandler;

using TC = Optimization::Optimizer::TerminationCondition;

using Utilities::FileHandling::CreateDir;
using Utilities::FileHandling::FileExists;
using Utilities::FileHandling::DeleteFile;
using Utilities::FileHandling::WriteLineToFile;

class DFTR : public Optimizer {

 public:
  DFTR(Settings::Optimizer *settings,
       Case *base_case,
       VarPropContainer *variables,
       Reservoir::Grid::Grid *grid,
       Logger *logger,
       CaseHandler *case_handler = nullptr,
       ConstraintHandler *constraint_handler = nullptr
  );

  TC IsFinished() override;
  TRFrame *getTRMod() { return trm_; };

 protected:
  void iterate() override;
  void handleEvaluatedCase(Case *c) override;

 private:
  TRFrame *trm_;

  // DFTR management
  void updateRadius();
  bool testCriticality();

  // FO-DFTR interface management
  bool submitTempCases();
  bool handleTempCases(Case *c);
  bool ensureImprPostProc();

  // DFTR constructor aux
  void setSettings(Settings::Optimizer *settings);
  void setLowerUpperBounds();
  void computeInitPts();
  void projToBnds(VectorXd *point);
  void createLogFile();
  void printIteration(double fval_current);

  // DFTR dbg
  bool F = false;
  bool T = true;
  double D = 0.0;
  VectorXd V = VectorXd::Zero(0);

  // TR core properties
  double infd_ = -std::numeric_limits<double>::infinity();
  double tr_init_rad_;
  double tr_tol_f_, tr_eps_c_, tr_eta_0_, tr_eta_1_; // Tols
  double tr_piv_thld_, tr_add_thld_, tr_xch_thld_; // Thesholds
  double tr_rad_max_, tr_rad_fac_, tr_rad_tol_; // Radii
  double tr_gamma_inc_, tr_gamma_dec_; // Gamma factors
  double tr_crit_mu_, tr_crit_omega_, tr_crit_beta_; // Criticality
  double tr_lower_bnd_, tr_upper_bnd_;
  VectorXd lb_, ub_; // Bounds [scalars from settings, lb/ub vectors]
  int tr_iter_max_, tr_num_init_x_, tr_rng_seed_; // Iter max + seed
  string tr_basis_, tr_init_smpln_, tr_prob_name_;

  int tr_lim_inf_ = 0;

  // TR management
  double rho_;         // Agreement factor>
  double ared_;        // Actual reduction>
  double sum_rho_;     // Cumulative rho>
  double sum_rho_sqr_; // Cumulative squared rho
  double delay_reduc_; // Delay reduction>
  double fmult_;       // Maximize <-> minimize multiplier

  // int n_initial_points_;

  // TR iteration
  VectorXd trial_point_;
  VectorXd trial_step_;
  double fval_trial_, pred_red_, cr_rad_ini_;

  int mchange_ = 0;
  bool iter_modl_fl_ = false;
  bool impr_modl_nx_ = false;

  critExecStat cr_stat_ = FAILED;
  string getCrStat() {
    if (cr_stat_ == ONGOING) { return "ONGOING";
    } else if (cr_stat_ == SUCCESS) { return "SUCCESS";
    } else if (cr_stat_ == FAILED) { return "FAILED"; }
  }

  // FO constructs
  Settings::Optimizer *settings_;
  Model::Properties::VarPropContainer *variables_;
  Case *base_case_;

  // Dbg constructs
  string md_ = "Optimization/optimizers/dftr";
  string cl_ = "DFTR";
  string im_ = "", wm_ = "", em_ = "";
  Settings::VerbParams vp_;
  QString dftr_log_;

  class ConfSmry : public Loggable {
   public:
    ConfSmry(DFTR *opt) { opt_ = opt; }
    LogTarget GetLogTarget() override {
      return LOG_SUMMARY;
    };

    QUuid GetId() override {
      return QUuid(); // Null UUID
    };
   private:
    DFTR *opt_;
  };

};

}
}

#endif //FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_DFTR_H_
