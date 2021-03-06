/***********************************************************
Created by bellout on 21.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
DFTR built from TrustRegionOptimization.cpp
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

#include "Utilities/printer.hpp"
#include "Utilities/stringhelpers.hpp"
#include "Utilities/random.hpp"

#include "DFTR.h"

namespace Optimization {
namespace Optimizers {

using OptMode = Settings::Optimizer::OptimizerMode;
// using CS = TRFrame::critExecStat; // remove later

DFTR::DFTR(Settings::Optimizer *settings,
           Case *base_case,
           VarPropContainer *variables,
           Reservoir::Grid::Grid *grid,
           Logger *logger,
           CaseHandler *case_handler,
           ConstraintHandler *constraint_handler)
  : Optimizer(settings, base_case, variables, grid,
              logger, case_handler, constraint_handler) {
  setSettings(settings);
  variables_ = variables;
  base_case_ = base_case;
  setLowerUpperBounds();
  trm_ = new TRFrame(lb_, ub_, base_case_, settings_);
  computeInitPts();
  createLogFile();
}

// __________________________________________________________
// ISFINISHED
TC DFTR::IsFinished() {
  TC tc = NOT_FINISHED;

  if (trm_->isModInit() && !trm_->getModPolys().empty()) {
    if (trm_->measureCrit().norm() < tr_tol_f_) {
      tc = DFTR_CRIT_NORM_1ST_ORD_LT_TOLF;

    } else if (tr_lim_inf_ > std::ceil(tr_iter_max_*.20)) {
      tc = DFTR_MAX_NUM_RHO_INF_MET;
    }
  }
  if (trm_->getRad() < tr_rad_tol_) {
    tc = MIN_STEP_LENGTH_REACHED;
  }
  if (iteration_ == tr_iter_max_) {
    tc = MAX_ITERS_REACHED;
  }
  if (tc != NOT_FINISHED) {
  }
  return tc;
}

// _________________________________________________________
// UPDATERADIUS
void DFTR::updateRadius() {
  // From time to time a step may end a bit outside the TR
  double step_sz = min(trm_->getRad(), trial_step_.lpNorm<Infinity>());
  double inc_fac = 0.0;

  auto A = (rho_ == infd_);
  auto B = (mchange_ == 4);
  trm_->dbg_->prntUpdateRad(1, A, B, rho_, tr_eta_1_, step_sz, inc_fac);

  // Radius update>
  if (rho_ > tr_eta_1_) {
    inc_fac = max(1.0, tr_gamma_inc_ * (step_sz / trm_->getRad()));
    trm_->setRad(min(inc_fac * trm_->getRad(), tr_rad_max_));
    trm_->dbg_->prntUpdateRad(2, F, F, trm_->getRad(), step_sz, inc_fac);

  } else if (iter_modl_fl_ && (A || B || (cr_stat_ == SUCCESS))) {
    // A good model should have provided a better point.
    // Reduce radius, since error is related to it:
    //        | f(x) - m(x) | < K*(radius)^2
    //
    // rho == -inf -> too short step size
    // mchange_ == 4 -> Couldn't add point, had to rebuild model

    if (trm_->getRad() <= 2 * tr_rad_tol_ / tr_gamma_dec_) {
      delay_reduc_ = delay_reduc_ + 1;
    } else {
      delay_reduc_ = 0;
    }

    if (delay_reduc_ >= 3 || delay_reduc_ == 0 || (cr_stat_ == SUCCESS)) {
      trm_->setRad(tr_gamma_dec_ * trm_->getRad());
      delay_reduc_ = 0;
    }

    trm_->dbg_->prntUpdateRad(3, iter_modl_fl_, F, trm_->getRad());
  } else {
    trm_->dbg_->prntUpdateRad(4, F, F, trm_->getRad());
  }
  iteration_++;
}

// ---------------------------------------------------------
// TESTCRITICALITY
bool DFTR::testCriticality() {
  //! -> Measure criticality
  auto cr_mod = trm_->measureCrit();
  trm_->dbg_->prntTstCrit(1, getCrStat(), cr_mod, trm_->getCurrPt(),
                           trm_->getCurrFval(), trm_->getRad());

  //! -> Crit step (possibly close to x*)
  if (cr_mod.norm() <= tr_eps_c_) {
    trm_->dbg_->prntTstCrit(2);

    if (cr_stat_ != ONGOING) { cr_rad_ini_ = trm_->getRad(); }

    cr_stat_ = trm_->critStep(cr_rad_ini_);
    if (cr_stat_ == ONGOING) {
      ensureImprPostProc();
      return true;
    }

    trm_->dbg_->prntTstCrit(3, getCrStat(), cr_mod);
    if (cr_mod.norm() < tr_tol_f_) { return true; } // check if not obsolete

  } else {
    cr_stat_ = FAILED;
  }
  return false;
}

// ---------------------------------------------------------
// ITERATE
void DFTR::iterate() {
  if (iteration_ == 0) {
    if (submitTempCases()) { return; }
  }
  assert(trm_->isModInit());

  // dbg?
  // auto x_curr = trm_->getCurrPt(); // VectorXd
  // auto f_curr = trm_->getCurrFval(); // double
  double err_mod = 0.0;

  if (impr_modl_nx_) {
    trm_->dbg_->prntIterBool(1);
    mchange_ = trm_->ensureImpr();
    impr_modl_nx_ = ensureImprPostProc();
  }

  // logic
  // repl pts comp/needed, impr pts comp/needed
  auto IC = (trm_->areImprPtsCd() && trm_->isImprNeeded());
  auto RC = (trm_->areReplPtsCd() && trm_->isReplNeeded());

  auto NIC = !trm_->areImprPtsCd() && !trm_->isImprNeeded();
  auto NRC = !trm_->areReplPtsCd() && !trm_->isReplNeeded();
  trm_->dbg_->prntIterBool(2, IC, RC, NIC, NRC);

  if (NIC && NRC) { // 1ST PART ITER LOGIC
    if (enable_logging_) { logger_->AddEntry(this); }

    // CODE ->
    if ((trm_->getRad() < tr_rad_tol_) || (iteration_ >= tr_iter_max_)) {
      trm_->dbg_->prntIterBool(3);
      return; // end of the algorithm

    } else {
      // Move among points that are part of the model
      // if (trm_->isLambdaPoised()) {
      trm_->dbg_->prntIterBool(4);
      trm_->moveToBestPt();
      trm_->computePolyMods();
      err_mod = trm_->checkInterp();
      // }

      printIteration(trm_->getCurrFval());
      if (testCriticality()) { return; };

      //! -> Check LambdaPoised
      trm_->dbg_->prntIterBool(5);
      iter_modl_fl_ = trm_->isLambdaPoised();

      //! -> Compute step
      trm_->dbg_->prntIterBool(6);
      tie(trial_point_, pred_red_) = trm_->solveTrSubprob();
      trial_step_ = trial_point_ - trm_->getCurrPt();

      trm_->dbg_->prntIterBool(7); // to be filled

      //! --
      if ((pred_red_ < tr_rad_tol_ * 1e-2)
        || (pred_red_ < tr_rad_tol_ * abs(trm_->getCurrFval()) && (trial_step_.norm()) < tr_rad_tol_)
        || (pred_red_ < tr_tol_f_ * abs(trm_->getCurrFval()) * 1e-3)) {

        trm_->dbg_->prntIterBool(8);

        rho_ = trm_->infd_;
        tr_lim_inf_++;
        mchange_ = trm_->ensureImpr();
        impr_modl_nx_ = ensureImprPostProc();

      } else {
        trm_->dbg_->prntIterBool(9); //! -> Evaluate objective at trial point>
        Case *new_case = new Case(base_case_);
        new_case->SetRealVarValues(trial_point_);
        case_handler_->AddNewCase(new_case);
        return;
      }

    } // <- CODE

  } else { // 2ND PART ITER LOGIC
    auto impr_cs = trm_->getImprCases();
    auto repl_cs = trm_->getReplCases();

    if (IC || RC) {
      trm_->dbg_->prntIterBool(10);

      if (impr_cs.empty() == 0){
        for (Case* c: impr_cs) {
          if (isImprovement(c)) { updateTentativeBestCase(c); }
          if (enable_logging_) { logger_->AddEntry(this); }
        }
      }

      if (repl_cs.empty() == 0){
        for (Case* c: repl_cs) {
          if (isImprovement(c)) { updateTentativeBestCase(c); }
          if (enable_logging_) { logger_->AddEntry(this); }
        }
      }

      // CODE -> Continue with model improvement
      mchange_ = trm_->ensureImpr();
      if ((mchange_ == 1) || (mchange_ == 2)) {
        iteration_--;
      }
      trm_->dbg_->prntIterBool(11);
      updateRadius();
      return;
      // <- CODE

    } else {
      trm_->dbg_->prntIterBool(12);
      if (trm_->isImprNeeded() && impr_cs.empty() == 0) {
        case_handler_->AddNewCases(impr_cs);
        return;
      }

      if (trm_->isReplNeeded() && repl_cs.empty() == 0) {
        case_handler_->AddNewCases(repl_cs);
        return;
      }
    }

  } // if (!IC && !RC)
  trm_->dbg_->prntIterBool(13);
}

// _________________________________________________________
// SUBMITTEMPCASES
bool DFTR::submitTempCases() {
  auto PN = (!trm_->areInitPtsCd() && !trm_->isModInit());
  auto PY = (trm_->areInitPtsCd() && !trm_->isModInit());

  // both repl, impr points computed
  auto IRC = (trm_->areImprPtsCd() && trm_->areReplPtsCd());
  // repl pts comp/needed, impr pts comp/needed
  auto IC = (trm_->areImprPtsCd() && trm_->isImprNeeded());
  auto RC = (trm_->areReplPtsCd() && trm_->isReplNeeded());
  trm_->dbg_->prntTempCases(1, PN, PY, IRC, IC, RC);

  // iteration_ == 0
  if (PN) {
    auto init_cases = trm_->getInitCases();
    case_handler_->AddNewCases(init_cases);
    trm_->dbg_->prntTempCases(2);
    return true;
  }

  // iteration_ == 0
  if (PY && !IRC) {
    bool is_model_changed = trm_->rebuildModel();
    trm_->setHasModChg(is_model_changed);
    trm_->moveToBestPt();
    trm_->computePolyMods();
    trm_->dbg_->prntTempCases(3);

    if(trm_->getNPts() < 2) { // trm_ has only 1 point
      trm_->dbg_->prntTempCases(4);
      mchange_ = trm_->ensureImpr();

      auto impr_cs = trm_->getImprCases();
      auto repl_cs = trm_->getReplCases();
      // Might be 0 if point_found in improveModelNfp is false.
      // We also need to handle exit_flag=3 or exit_flag=4

      if (trm_->isImprNeeded() && impr_cs.empty() == 0) {
        case_handler_->AddNewCases(impr_cs);
        return true;
      }

      if (trm_->isReplNeeded() && repl_cs.empty() ==0) {
        case_handler_->AddNewCases(repl_cs);
        return true;
      }

    } else { // trm_ has two points
      trm_->setIsModInit(true);
      iteration_++;
      trm_->dbg_->prntTempCases(5);
    }
  }

  // iteration_ == 0
  if (PY && (IC || RC)) {
    mchange_ = trm_->ensureImpr();
    if(!trm_->isImprNeeded() || !trm_->isReplNeeded()) {
      trm_->setIsModInit(true);
      iteration_++;
      trm_->dbg_->prntTempCases(6);
    }
    if (enable_logging_) { logger_->AddEntry(this); }
  }
  trm_->dbg_->prntTempCases(7);
  return false;
}

// _________________________________________________________
// HANDLETEMPCASES
bool DFTR::handleTempCases(Case *c) {
  auto A = case_handler_->QueuedCases().empty();
  auto B = case_handler_->CasesBeingEvaluated().empty();

  auto PN = !trm_->areInitPtsCd() && !trm_->isModInit();
  auto PY = trm_->areInitPtsCd() && !trm_->isModInit();

  // replacement, improvement points needed
  auto IRN = (trm_->isImprNeeded() || trm_->isReplNeeded());

  if (PN && iteration_ == 0) {

    trm_->addTempInitCase(c); // Collect init case
    if (A && B) { // All init cases evaluated
      trm_->submitTempInitCases();
      trm_->setAreInitPtsCd(true);
    }
    return true;

  } else if ((PY && iteration_ == 0) || (IRN && iteration_ > 0)) {

    if (trm_->isImprNeeded()) { // If impr needed
      trm_->addTempImprCase(c); // Collect impr case
      if (A && B) { // All impr cases evaluated
        trm_->submitTempImprCases();
        trm_->setAreImprPtsCd(true);
      }
      return true;
    }

    if (trm_->isReplNeeded()) { // If repl needed
      trm_->addTempReplCase(c); // Collect repl case
      if (A && B) { // all repl cases evaluated
        trm_->submitTempReplCases();
        trm_->setAreReplPtsCd(true);
      }
      return true;
    }
  }
  return false;
}

// _________________________________________________________
// HANDLEEVALUATEDCASE
void DFTR::handleEvaluatedCase(Case *c) {
  if (isImprovement(c)) { updateTentativeBestCase(c); }

  if (!handleTempCases(c)) {
    fval_trial_ = fmult_ * c->objf_value();
    ared_ = trm_->getCurrFval() - fval_trial_; // Actual reduction
    rho_ = ared_ / (pred_red_); // Agreement factor

    if (rho_ > tr_eta_1_) { // Acceptance of the trial point
      if (isImprovement(c)) { // Successful iteration
        updateTentativeBestCase(c);
      }

      // Including new point as TR center
      mchange_ = trm_->changeTrCenter(trial_point_, fval_trial_);
    } else {
      auto point_added = trm_->tryToAddPt(trial_point_, fval_trial_);
      if (!point_added) { impr_modl_nx_ = true; }
    }
    sum_rho_ += rho_;
    sum_rho_sqr_ += pow(rho_, 2);

    if (!impr_modl_nx_) { updateRadius(); }
  }
}

// _________________________________________________________
// ENSUREIMPRPOSTPROC
bool DFTR::ensureImprPostProc(){
  auto impr_cases = trm_->getImprCases(); // improve model
  auto repl_cases = trm_->getReplCases(); // replace point

  if (trm_->isImprNeeded() && impr_cases.empty() == 0) {
    case_handler_->AddNewCases(impr_cases);
    trm_->dbg_->prntEnsImprPostProc(1);
    return true;
  }

  if (trm_->isReplNeeded() && repl_cases.empty() == 0) {
    case_handler_->AddNewCases(repl_cases);
    trm_->dbg_->prntEnsImprPostProc(2);
    return true;
  }

  auto E = (impr_cases.empty() && repl_cases.empty());
  auto F = (cr_stat_ == ONGOING);
  trm_->dbg_->prntEnsImprPostProc(3, 0, 0, 0, F,
    impr_cases.empty(), repl_cases.empty());

  if (F && E) {
    if (vp_.vOPT > 1) {
      wm_ = "[@ensureImprPostProc] Crit step did not generate a case.";
      ext_warn(wm_, md_, cl_, vp_.lnw);
    }
  }

  //! Only update radius if not within criticality execution
  if ((impr_cases.empty()) && (repl_cases.empty()) && !F) {
    trm_->dbg_->prntEnsImprPostProc(4);
    updateRadius();
  }
  return false;
}

// _________________________________________________________
// SETSETTINGS
void DFTR::setSettings(Settings::Optimizer *settings) {
  settings_ = settings;
  tr_init_rad_   = settings_->parameters().tr_init_rad;
  tr_tol_f_      = settings_->parameters().tr_tol_f;
  tr_eps_c_      = settings_->parameters().tr_eps_c;
  tr_eta_0_      = settings_->parameters().tr_eta_0;
  tr_eta_1_      = settings_->parameters().tr_eta_1;
  tr_piv_thld_   = settings_->parameters().tr_piv_thld;
  tr_add_thld_   = settings_->parameters().tr_add_thld;
  tr_xch_thld_   = settings_->parameters().tr_xch_thld;
  tr_rad_max_    = settings_->parameters().tr_rad_max;
  tr_rad_fac_    = settings_->parameters().tr_rad_fac;
  tr_rad_tol_    = settings_->parameters().tr_rad_tol;
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

  // Log base case
  if (enable_logging_) { logger_->AddEntry(this); }

  // TR management
  rho_ = 0;
  sum_rho_ = 0;
  sum_rho_sqr_ = 0;
  delay_reduc_ = 0;
}

// _________________________________________________________
// SETLOWERUPPERBNDS
void DFTR::setLowerUpperBounds() {
  lb_.conservativeResize(base_case_->GetRealVarIdVector().size());
  ub_.conservativeResize(base_case_->GetRealVarIdVector().size());

  if (constraint_handler_ != nullptr) { // All actual cases

    if (constraint_handler_->HasBoundaryConstraints()) {
      // Use constraint_handler lower/upper bounds
      lb_ = constraint_handler_->GetLowerBounds(
        base_case_->GetRealVarIdVector());
      ub_ = constraint_handler_->GetUpperBounds(
        base_case_->GetRealVarIdVector());

    } else {
      // Use lower/upper bounds specified in settings
      if (settings_->parameters().tr_lower_bnd > 0
        && settings_->parameters().tr_upper_bnd > 0) {
        lb_.fill(settings_->parameters().tr_lower_bnd);
        ub_.fill(settings_->parameters().tr_upper_bnd);

      } else {
        wm_ = "Lower/upper bounds for DF-TR algorithm not specified.";
        ext_warn(wm_, md_, cl_, vp_.lnw);
        throw std::runtime_error(wm_);
      }
    }

  } else { // constraint_handler_ == nullptr in unit tests
    lb_.fill(settings_->parameters().tr_lower_bnd);
    ub_.fill(settings_->parameters().tr_upper_bnd);
  }
}

// _________________________________________________________
// COMPUTEINITPTS
void DFTR::computeInitPts() {
  trm_->dbg_->prntProgInit(1, lb_, ub_);

  trm_->addInitCase(base_case_);
  trm_->setXDim(variables_->ContinuousVariableSize());
  int n_cont_vars = trm_->getXDim();
  auto init_pt = base_case_->GetRealVarVector();

  //! Find another point since only one initial guess provided
  if (settings_->parameters().tr_num_init_x == -1) {
    if (settings_->parameters().tr_init_smpln == "Random") {

      //! Find (random) 2nd initial point
      auto rng = get_random_generator(tr_rng_seed_);
      VectorXd scnd_pt(n_cont_vars, 1);
      scnd_pt.setZero(n_cont_vars);
      for (int i = 0; i < n_cont_vars; ++i) {
        scnd_pt(i) = random_double(rng, 0, 1);
      }
      trm_->dbg_->prntProgInit(3, init_pt, scnd_pt);

      //! Second point must not be too close
      while (scnd_pt.lpNorm<Infinity>() < tr_piv_thld_) {
        scnd_pt *= 2;
        trm_->dbg_->prntProgInit(2, scnd_pt);
      }
      scnd_pt = scnd_pt.array() - .5;
      auto scnd_bs = scnd_pt;
      scnd_pt *= tr_init_rad_; //! scaling measure
      scnd_pt = init_pt + scnd_pt;

      if(lb_.size() > 0 && ub_.size() > 0) { projToBnds(&scnd_pt); }
      trm_->dbg_->prntProgInit(3, init_pt, scnd_pt, tr_init_rad_);

      //! Append case corresponding to 2nd init point
      Case *scnd_case = new Case(base_case_);
      scnd_case->SetRealVarValues(scnd_pt);

      // Hack: Comment this line to test improveModelNfp()
      // -> currently, this leads to a crash
      // This will be done later in the test code
      trm_->addInitCase(scnd_case);

    } else {
      //TODO: implement other sampling method such as Uniform.
    }

  } else {
    //TODO: get other initial points from the parameters
  }

  trm_->setInitRad(tr_init_rad_);
}

// _________________________________________________________
// PROJECTTOBOUNDS
void DFTR::projToBnds(VectorXd *point) {
  if (ub_.size() > 0) { *point = ub_.cwiseMin(*point); }
  if (lb_.size() > 0) { *point = lb_.cwiseMax(*point); }
}

// _________________________________________________________
// CREATELOGFILE
void DFTR::createLogFile() {
  QString out_dir = logger_->getOutputDir() + "/tr-dfo-out";
  if (out_dir.length() > 0) { CreateDir(out_dir, vp_); }
  QString fn = "/dftr_" + QString::fromStdString(tr_prob_name_) + ".csv";
  dftr_log_ = out_dir + fn;

  // Delete existing logs if --force flag is on
  if (FileExists(dftr_log_, vp_)) {
    DeleteFile(dftr_log_, vp_);
  }

  WriteLineToFile(fn, dftr_log_);
  const QString header = "iter        fval         rho      radius  pts";
  WriteLineToFile(header, dftr_log_);
}

// _________________________________________________________
// PRINTITERATION
void DFTR::printIteration(double fval_current) {
  stringstream ss;
  ss << setw(4) << right << iteration_ << setprecision(3)
     << setw(12) << scientific << right << fval_current
     << setw(12) << scientific << right << rho_
     << setw(12) << scientific << right << trm_->getRad()
     << setw(5) << right << trm_->getNPts();
  string str = ss.str();
  WriteLineToFile(QString::fromStdString(str), dftr_log_);
}

}
}

