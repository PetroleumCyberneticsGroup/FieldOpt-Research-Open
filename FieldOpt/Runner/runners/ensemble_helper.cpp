/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
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

#include <iostream>
#include "Utilities/math.hpp"
#include "Utilities/random.hpp"
#include "ensemble_helper.h"

#include "Utilities/verbosity.h"
#include "Utilities/printer.hpp"

namespace Runner {

EnsembleHelper::EnsembleHelper() {
  current_case_ = 0;
  rzn_queue_ = std::vector<std::string>();
  rzn_busy_ = std::vector<std::string>();
}

EnsembleHelper::EnsembleHelper(const Settings::Ensemble &ensemble,
                               int rng_seed) {
  ensemble_ = ensemble;
  current_case_ = 0;
  rzn_queue_ = std::vector<std::string>();
  rzn_busy_ = std::vector<std::string>();
  rng_ = get_random_generator(rng_seed*3);
  for (std::string alias : ensemble.GetAliases()) {
    assigned_workers_[alias] = std::vector<int>();
  }

  // \todo  Use SetNSelect() to set n_select_
  n_select_ = ensemble.NSelect();

  // dbg -> enables rlz subset (random iter needs proper debug)
  // n_select_ = ensemble.GetAliases().size() - 5;
  // cout << "ensemble.GetAliases().size(): " << ensemble.GetAliases().size() << endl;
  // cout << "n_select_:                    " << n_select_ << endl;

  if (n_select_ > ensemble.GetAliases().size()) {
    n_select_ = ensemble.GetAliases().size();
    if ( VERB_RUN >= 3 ) {
      cout << "Override: n_select_ now =: # of rlzs in *.ens "
      << "file (=" << n_select_ << ")" << endl;
    }
  }
  assert(n_select_ <= ensemble.GetAliases().size());
  assert(n_select_ > 0);
}

void EnsembleHelper::SetActiveCase(Optimization::Case *c) {
  if (!IsCaseDone()) {
    std::cerr
        << "ERROR: Unable to set new active case "
           "before the previous case is done." << std::endl;
    throw std::runtime_error("Error in EnsembleHelper.");
  }

  current_case_ = c;
  selectRealizations();
  eval_start_time_ = std::chrono::high_resolution_clock::now();
}

bool EnsembleHelper::IsCaseDone() const {
  return rzn_queue_.empty() && rzn_busy_.empty();
}

bool EnsembleHelper::IsCaseAvailableForEval() const {
  return !rzn_queue_.empty();
}

void EnsembleHelper::GetCaseForEval(Optimization::Case &c) {

  std::string next_alias = rzn_queue_.back();
  rzn_queue_.pop_back();
  rzn_busy_.push_back(next_alias);

  c.SetEnsembleRealization(
      QString::fromStdString(next_alias));
}

Optimization::Case *EnsembleHelper::GetCaseForEval() {
  if (!IsCaseAvailableForEval()) {
    std::cerr
        << "ERROR: No more cases are available for evaluation."
        << std::endl;
    throw std::runtime_error("Error in EnsembleHelper.");
  }
  auto case_copy = new Optimization::Case(current_case_);
  std::string next_alias = rzn_queue_.back();
  rzn_queue_.pop_back();
  rzn_busy_.push_back(next_alias);

  // TODO: This is a duplicate case that
  //  wont get deleted, i.e. a MEMORY LEAK.
  case_copy->SetEnsembleRealization(
      QString::fromStdString(next_alias));

  return case_copy;
}

void EnsembleHelper::SubmitEvaluatedRealization(Optimization::Case *c) {

  long alias_pos =
      distance(rzn_busy_.begin(),
               find(rzn_busy_.begin(), rzn_busy_.end(),
                    c->GetEnsembleRlz().toStdString()));

  if (alias_pos < 0 || alias_pos >= rzn_busy_.size()) {
    std::cerr
        << "ERROR: Unable to find alias in list of busy realizations."
        << std::endl;
    throw std::runtime_error("Error in EnsembleHelper.");
  }

  if (c->state.eval == Optimization::Case::CaseState::E_DONE) {

    current_case_->SetRealizationOfv(
      c->GetEnsembleRlz(),
        c->objf_value());

  } else {
    std::cout << "WARNING: Ensemble realization case "
              << c->GetEnsembleRlz().toStdString()
              << " was not successfully evaluated. It will "
                 "not be further considered."
              << std::endl;
  }
  rzn_busy_.erase(rzn_busy_.begin() + alias_pos);
}

Optimization::Case *EnsembleHelper::GetEvaluatedCase() {
  if (!IsCaseDone()) {
    std::cerr
        << "ERROR: Unable to get case before all "
           "selected realizations have been evaluated."
        << std::endl;
    throw std::runtime_error("Error in EnsembleHelper.");
  }
  rzn_queue_ = std::vector<std::string>();
  rzn_busy_ = std::vector<std::string>();

  current_case_->set_objf_value(
    current_case_->GetEnsembleAverageOfv());

  if ( VERB_OPT >= 3 ) {
    cout << "submitting evaluated ensemble case; objf = "
         << current_case_->objf_value() << endl;
  }

  // -------------------------------------------------------
  auto eval_end_time = std::chrono::high_resolution_clock::now();
  auto time_diff =
      std::chrono::duration_cast<std::chrono::milliseconds>
          (eval_end_time - eval_start_time_);

  // -------------------------------------------------------
  current_case_->SetSimTime(time_diff.count() / 1000);
  current_case_->state.eval =
      Optimization::Case::CaseState::EvalStatus::E_DONE;

  return current_case_;
}

void EnsembleHelper::selectRealizations() {
  auto all_aliases = ensemble_.GetAliases();

  if (n_select_ == all_aliases.size()) {
    if (VERB_RUN >=2) {
      Printer::ext_info("Selecting all realizations",
                        "Runner", "EnsembleHelper");
    }
    for (auto alias : all_aliases) {
      rzn_queue_.push_back(alias);
    }

  } else {
    // \todo debug: make sure this function does what it is supposed
    auto indices = unique_random_integers(rng_, 0,
                                          all_aliases.size() - 1, n_select_);
    if (VERB_RUN >=2) {
      // /todo make proper printing-of-vector in printer.hpp
      stringstream ss; for (int ii=0; ii < indices.size(); ii++) {
        ss << indices[ii] << " ";
      }
      Printer::ext_info("Selecting subset of realizations: " + ss.str(),
                        "Runner", "EnsembleHelper");
    }
    for (auto idx : indices) {
      rzn_queue_.push_back(all_aliases[idx]);
    }
  }
}

Settings::Ensemble::Realization EnsembleHelper::GetRlz(
    const std::string &alias) const {
  return ensemble_.GetRealization(alias);
}

Settings::Ensemble::Realization EnsembleHelper::GetBaseRealization() const {

  auto base_alias = ensemble_.GetAliases()[0];
  if ( VERB_RUN >= 2 ) {
    std::stringstream ss;
    ss << "Ensemble base rlz: " << base_alias << endl;
    cout << ss.str();
  }
  return ensemble_.GetRealization(base_alias);
}

int EnsembleHelper::NBusyCases() const {
  return rzn_busy_.size();
}

int EnsembleHelper::NQueuedCases() const {
  return rzn_queue_.size();
}

std::string EnsembleHelper::GetStateString() const {
  std::stringstream str;
  str << "EnsembleHelper: ";
  if (current_case_ == 0) {
    str << "Case not set.";

  } else {
    str << "Current ensemble case done: "
        << (IsCaseDone() ? "Yes" : "No") << endl
        << "N. Queued Cases: " << NQueuedCases() << string(9, ' ')
        << "N. Busy Cases:   " << NBusyCases();
  }
  return str.str();
}

bool EnsembleHelper::HasAssignedWorkers(const std::string &alias) const {
  return assigned_workers_.at(alias).size() > 0;
}

int EnsembleHelper::GetAssignedWorker(const std::string &alias,
    std::vector<int> free_workers) {

  for (int rank : free_workers) {
    // Check if realization has been assigned
    // to one of the free workers
    if (isAssignedToWorker(alias, rank)) {
      return rank; // If it has, return that rank
    }
  }

  // If not, assign one of the free workers to the realization
  return AssignNewWorker(alias, free_workers);
}

int EnsembleHelper::AssignNewWorker(const std::string &alias,
    std::vector<int> free_workers) {

  assert(free_workers.size() > 0);
  auto loads = workerLoads(free_workers);

  int least_loaded_worker = loads.front().first;
  assigned_workers_[alias].push_back(least_loaded_worker);

  return least_loaded_worker;
}

bool EnsembleHelper::isAssignedToWorker(const std::string &alias,
    const int rank) const {

  auto workers = assigned_workers_.at(alias);
  if (std::find(workers.begin(), workers.end(), rank) != workers.end()) {
    return true;

  } else {
    return false;
  }
}

int EnsembleHelper::nRealizationsAssignedToWorker(const int &rank) const {

  int count = 0;
  for (std::pair<std::string, std::vector<int> > assigments : assigned_workers_) {
    if (isAssignedToWorker(assigments.first, rank)) {
      count++;
    }
  }
  return count;
}

std::vector< pair<int, int> >
EnsembleHelper::workerLoads(std::vector<int> free_workers) const {

  std::vector< pair<int, int> > loads;
  for (int rank : free_workers) {
    loads.push_back(
        std::pair<int, int>(rank,
                            nRealizationsAssignedToWorker(rank))
    );
  }

  std::sort(loads.begin(), loads.end(),
            []( pair<int, int> l1, pair<int, int> l2)
            { return l1 < l2; } );

  return loads;
}

} // End ensemble_helper
