/***********************************************************
Copyright (C) 2019
Einar J.M. Baumann <einar.baumann@gmail.com>
Created by einar on Jan. 14. 2019

Modified 2021 Mathias Bellout
<chakibbb.pcg@gmail.com>

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

#include "Optimization/optimizers/SPSA.h"
#include "Utilities/stringhelpers.hpp"
#include <cmath>

namespace Optimization {
namespace Optimizers {


SPSA::SPSA(Settings::Optimizer *settings,
     Case *base_case,
     Model::Properties::VarPropContainer *variables,
     Reservoir::Grid::Grid *grid,
     Logger *logger,
     CaseHandler *case_handler,
     Constraints::ConstraintHandler *constraint_handler)
    : Optimizer(settings, base_case, variables, grid, logger, case_handler, constraint_handler) {

  D_ = base_case->GetRealVarVector().size();
  gen_ = get_random_generator(settings->parameters().rng_seed);
  //gen_ = get_random_generator(1);

  auto params = settings->parameters();
  max_iterations_ = params.spsa_max_iterations;
  init_step_magnitude_ = params.spsa_init_step_magnitude;
  alpha_ = params.spsa_alpha;
  gamma_ = params.spsa_gamma;
  A_ = params.spsa_A;
  a_ = params.spsa_a;
  c_ = params.spsa_c;
  assert(a_ != 0.0 || init_step_magnitude_ != 0.0);

  estimate_ = new Case(base_case);
  perturbations_evaluated_ = false;
  perturbations_valid_ = false;

  if (constraint_handler_ != nullptr) { // All actual cases
    ub_ = constraint_handler_->GetUpperBounds(base_case->GetRealVarIdVector());
    lb_ = constraint_handler_->GetLowerBounds(base_case->GetRealVarIdVector());
  } else { // constraint_handler_ == nullptr in unit tests
    ub_.resize(base_case->GetRealVarIdVector().size());
    lb_.resize(base_case->GetRealVarIdVector().size());
    ub_.fill(settings->parameters().upper_bound);
    lb_.fill(settings->parameters().lower_bound);
  }
}

Optimization::Optimizer::TerminationCondition SPSA::IsFinished() {
  if (case_handler_->CasesBeingEvaluated().size() > 0 || case_handler_->QueuedCases().size() > 0) {
    return NOT_FINISHED;
  }
  if (iteration_ >= max_iterations_) {
    im_ = "Max iterations reached. Terminating";
    ext_info(im_, "Optimization", "SPSA");
    logger_->AddEntry(this);
    return MAX_ITERS_REACHED;
  }
  else {
    return NOT_FINISHED;
  }
}

void SPSA::handleEvaluatedCase(Case *c) {
  if (isImprovement(c)) {
    updateTentativeBestCase(c);
    if (vp_.vOPT >= 2) {
      im_ = "Found new tentative best case in iteration ";
      im_ += num2str(iteration_) + ": " + num2str(tentative_best_case_->objf_value());
      Printer::ext_info(im_, "Optimization", "SPSA");
    }
  }
  bool first_done = perturbations_.first->state.eval == Case::CaseState::E_DONE || perturbations_.first->state.eval == Case::CaseState::E_BOOKKEEPED;
  bool second_done = perturbations_.second->state.eval == Case::CaseState::E_DONE || perturbations_.second->state.eval == Case::CaseState::E_BOOKKEEPED;

  if (first_done && second_done) {
    perturbations_evaluated_ = true;
    perturbation_evaluation_failed_ = false;
    if (vp_.vOPT >= 3) {
      im_ = "Both perturbations evaluated in iteration " + num2str(iteration_);
      ext_info(im_, "Optimization", "SPSA");
    }
    updateGradient();
    if (iteration_ == 1 && init_step_magnitude_ != 0.0) {
      if (vp_.vOPT >= 3) {
        im_ = "First iteration done. Computing a.";
        ext_info(im_, "Optimization", "SPSA"); }
      compute_a();
    }
    update_a_k();
    updateEstimate();
  }
  else if ( ! first_done && ! second_done ) {
      std::stringstream ss;
      ss << "Non-successfully evaluated case returned (codes: " 
          << perturbations_.first->state.eval << ", " << perturbations_.second->state.eval
          << "). Attempting to recover.";
      Printer::ext_warn(ss.str(), "Optimization", "SPSA");
      perturbation_evaluation_failed_ = true;
  }
}

void SPSA::iterate() {
  logger_->AddEntry(this);

  if (iteration_ >= max_iterations_) {
    ext_info("Reached max iterations.", "Optimization", "SPSA");
    return;
  }

  if ( ! perturbations_evaluated_ && iteration_ != 0) {
    im_ = "Iterate called before perturbations were evaluated.";
    ext_warn(im_, "Optimization", "SPSA");
    return;
  }

  iteration_++;
  if (vp_.vOPT >= 3) {
    im_ = "Starting iteration " + Printer::num2str(iteration_);
    ext_info(im_, "Optimization", "SPSA");
  }

  update_c_k();

  int max_attempts = 1*D_;
  int attempt = 1;
  perturbations_valid_ = false;
  auto first = new Case(estimate_);
  auto second = new Case(estimate_);
  while (!perturbations_valid_ && attempt <= max_attempts) {
    updateSPVector();
    createPerturbations(first, second);
    attempt++;
  }
  perturbations_ = std::pair<Case *, Case *>(first, second);

  if (!perturbations_valid_) {
    Printer::ext_warn("Unable to generate valid pair of perturbations after " + Printer::num2str(max_attempts) + " attempts.", "Optimization", "SPSA");
    if ((iteration_ - tentative_best_case_iteration_) > 0.1*max_iterations_) {
      wm_ = "No new best case found in the last (0.1*max_iterations) iterations.";
      wm_ += " Setting estimate to tentative best case.";
      ext_warn(wm_, "Optimization", "SPSA");

      estimate_->SetRealVarValues(tentative_best_case_->GetRealVarVector());
      tentative_best_case_iteration_ = iteration_;
    }
    Printer::ext_warn("Skipping iteration." + Printer::num2str(iteration_), "Optimization", "SPSA");
    iterate();
  }
  else {
    case_handler_->AddNewCase(perturbations_.first);
    case_handler_->AddNewCase(perturbations_.second);
    perturbations_evaluated_ = false;
  }
}

void SPSA::update_a_k()
{
  a_k_ = a_ / pow((A_ + iteration_), alpha_);
  if (vp_.vOPT >= 3) {
    im_ = "Updated gain sequence a_k: " + Printer::num2str(a_k_);
    ext_info(im_, "Optimization", "SPSA");
  }
}

void SPSA::update_c_k()
{
  c_k_ = c_ / pow(iteration_, gamma_);
  if (vp_.vOPT >= 1) {
    Printer::ext_info("Updated gain sequence c_k: " + Printer::num2str(c_k_), "Optimization", "SPSA");
  }
}

void SPSA::updateSPVector()
{
  if (vp_.vOPT >= 4) {
    ext_info("Updating SP vector.", "Optimization", "SPSA");
  }
  delta_k_ = random_symmetric_bernoulli_eigen(gen_, D_);
}

void SPSA::createPerturbations(Case *first, Case *second)
{
  if (vp_.vOPT >= 5) {
    im_ = "Perturbation: +- " + eigenvec_to_str(c_k_ * delta_k_);
    ext_info(im_, "Optimization", "SPSA");
  }
  first->SetRealVarValues( estimate_->GetRealVarVector() + c_k_ * delta_k_ );
  second->SetRealVarValues(estimate_->GetRealVarVector() - c_k_ * delta_k_ );

  bool first_ok;
  bool second_ok;
  if (constraint_handler_ != nullptr) { // All actual cases
    first_ok = constraint_handler_->CaseSatisfiesConstraints(first);
    second_ok = constraint_handler_->CaseSatisfiesConstraints(second);
  }

  if ( first_ok && second_ok ) {
    perturbations_valid_ = true;
  }
  else if ( !first_ok && second_ok ) {
    if (vp_.vOPT >= 1) {
      im_ = "Positive perturbation violates constraints.";
      ext_info(im_, "Optimization", "SPSA");
    }
    first->SetRealVarValues( estimate_->GetRealVarVector() );
    second->SetRealVarValues(estimate_->GetRealVarVector() - 2.0*(c_k_ * delta_k_));

    if (constraint_handler_ != nullptr) { // All actual cases
      second_ok = constraint_handler_->CaseSatisfiesConstraints(second);
    }

    if (second_ok) {
      if (vp_.vOPT >= 1) {
        ext_info("Using estimate.", "Optimization", "SPSA");
      }
      perturbations_valid_ = true;
    } else {
      if (vp_.vOPT >= 1) {
        ext_info("Unable to correct.", "Optimization", "SPSA");
      }
      perturbations_valid_ = false;
    }
  }
  else if ( first_ok && !second_ok ) {
    if (vp_.vOPT >= 1) {
      im_ = "Negative perturbation violates constraints.";
      ext_info(im_, "Optimization", "SPSA");
    }
    second->SetRealVarValues( estimate_->GetRealVarVector() );
    first->SetRealVarValues( estimate_->GetRealVarVector() + 2.0*(c_k_ * delta_k_) );

    if (constraint_handler_ != nullptr) { // All actual cases
      first_ok = constraint_handler_->CaseSatisfiesConstraints(first);
    }

    if (first_ok) {
      if (vp_.vOPT >= 1) {
        ext_info("Using estimate.", "Optimization", "SPSA");
      }
      perturbations_valid_ = true;
    } else {
      if (vp_.vOPT >= 1) {
        ext_info("Unable to correct.", "Optimization", "SPSA");
      }
      perturbations_valid_ = false;
    }
  }
  else {
    if (vp_.vOPT >= 1) {
      im_ = "Both perturbations violate constraints.";
      ext_info(im_, "Optimization", "SPSA");
    }
    for (int i=0; i < D_; ++i) {
      double first_val = first->GetRealVarVector()[i];
      if (first_val > ub_[i] || first_val < lb_[i]) {
          delta_k_[i] = -1 * delta_k_[i];
      }
    }
    first->SetRealVarValues( estimate_->GetRealVarVector() + 2.0*(c_k_ * delta_k_) );
    second->SetRealVarValues( estimate_->GetRealVarVector() );

    if (constraint_handler_ != nullptr) { // All actual cases
      first_ok = constraint_handler_->CaseSatisfiesConstraints(first);
      second_ok = constraint_handler_->CaseSatisfiesConstraints(second);
    }

    if (first_ok && second_ok) {
      if (vp_.vOPT >= 1) Printer::ext_info("Managed to correct perturbations.", "Optimization", "SPSA");
      perturbations_valid_ = true;
    }
    else {
      if (vp_.vOPT >= 1) Printer::ext_info("Unable to correct perturbations.", "Optimization", "SPSA");
      perturbations_valid_ = false;
    }

  }
}

void SPSA::updateGradient() {
  if (vp_.vOPT >= 3) {
    ext_info("Updating gradient.", "Optimization", "SPSA");
  }
  double yplus = perturbations_.first->objf_value();
  double yminus = perturbations_.second->objf_value();
  if (mode_ == Settings::Optimizer::OptimizerMode::Maximize) {
    yplus = yplus * -1;
    yminus = yminus * -1;
  }
  g_k_.resize(D_);
  for (int i = 0; i < D_; ++i) {
    g_k_[i] = (yplus - yminus) / (2 * c_k_ * delta_k_[i]);
  }
  if (vp_.vOPT >= 4) {
    im_ = "Updated gradient vector: " + eigenvec_to_str(g_k_);
    ext_info(im_, "Optimization", "SPSA");
  }
}

void SPSA::updateEstimate() {
  if (vp_.vOPT >= 3) Printer::ext_info("Updating estimate.", "Optimization", "SPSA");

  if ((iteration_ - tentative_best_case_iteration_) > 0.1*max_iterations_) {
    wm_ = "No new best case found in the last (0.2*max_iterations) iterations. ";
    wm_ += "Setting estimate to tentative best case.";
    ext_warn(wm_, "Optimization", "SPSA");
    estimate_->SetRealVarValues(tentative_best_case_->GetRealVarVector());
    tentative_best_case_iteration_ = iteration_;
  } else {
    if (vp_.vOPT >= 4) {
      im_ = "Taking step: " + eigenvec_to_str(-1.0 * a_k_ * g_k_);
      ext_info(im_, "Optimization", "SPSA");
    }
    estimate_->SetRealVarValues(estimate_->GetRealVarVector() - a_k_ * g_k_);

    if (constraint_handler_ != nullptr) { // All actual cases
      constraint_handler_->SnapCaseToConstraints(estimate_);
    }
  }
}

void SPSA::compute_a()
{
  a_ = init_step_magnitude_ / abs(g_k_.mean()) * pow(A_ + 1.0, alpha_);
  if (a_ < 0.0) a_ = 0.0;
  if (a_ > 1.0) a_ = 1.0;
  if (vp_.vOPT >= 2) {
    Printer::ext_info("Computed parameter a: " + Printer::num2str(a_), "Optimization", "SPSA");
  }
}

}
}
