/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>
Created by einar on 3/27/17.

Modified 2017-2021 Mathias Bellout
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

#include "Utilities/math.hpp"
#include "Utilities/random.hpp"
#include "RGARDD.h"

namespace Optimization {
namespace Optimizers {

RGARDD::
RGARDD(Settings::Optimizer *settings,
       Optimization::Case *base_case,
       Model::Properties::VarPropContainer *variables,
       Reservoir::Grid::Grid *grid,
       Logger *logger,
       CaseHandler *case_handler,
       Constraints::ConstraintHandler *constraint_handler)
  : GeneticAlgorithm(settings, base_case, variables, grid,
                     logger, case_handler, constraint_handler) {
  vp_ = settings->verbParams();

  // Need an even number of chromosomes
  assert(population_size_ % 2 == 0);

  if (settings->parameters().discard_parameter < 0) {
    discard_parameter_ = 1.0/population_size_;
  } else {
    discard_parameter_ = settings->parameters().discard_parameter;
  }

  stagnation_limit_ = settings->parameters().stagnation_limit;
  mating_pool_ = population_;

  if (enable_logging_) {
    logger_->AddEntry(this);
    logger_->AddEntry(new ConfigurationSummary(this));
  }
}

void RGARDD::iterate() {
  if (!case_handler_->QueuedCases().empty()
    || !case_handler_->CasesBeingEvaluated().empty()) {
    wm_ = "Iteration requested while evaluation queue is not empty. Skipping call.";
    ext_warn(wm_, md_, cl_, vp_.lnw);
    return;
  }

  if (vp_.vOPT >= 2) { print_state("Iterating with RGARDD"); }
  if (enable_logging_) { logger_->AddEntry(this); }

  population_ = sortPopulation(population_);

  if (is_stagnant()) {
    if (vp_.vOPT >= 1) {
      im_ = "The population has stagnated in generation";
      im_ += num2str(iteration_) + ". Repopulating";
      ext_info(im_, md_, cl_, vp_.lnw);
    }
    repopulate();
    iteration_++;
    return;
  }

  mating_pool_ = selection(population_);
  double r;
  for (int i = 0; i < population_size_ / 2; ++i) {
    r = random_double(gen_);
    Chromosome p1 = mating_pool_[i]; // Parent 1
    Chromosome p2 = mating_pool_[population_size_ / 2 + i]; // Parent 2
    vector<Chromosome> offspring;

    if (r > p_crossover_ && p1.rea_vars != p2.rea_vars) {
      offspring = crossover(vector<Chromosome>{p1, p2});
    } else {
      offspring = mutate(vector<Chromosome>{p1, p2});
    }
    offspring[0].createNewCase();
    offspring[1].createNewCase();

    case_handler_->AddNewCase(offspring[0].case_pointer);
    case_handler_->AddNewCase(offspring[1].case_pointer);
    mating_pool_[i] = offspring[0];
    mating_pool_[population_size_ / 2 + i] = offspring[1];
  }
  iteration_++;
}

void RGARDD::handleEvaluatedCase(Case *c) {
  int index = -1;
  for (int i = 0; i < mating_pool_.size(); ++i) {
    if (mating_pool_[i].case_pointer == c) {
      index = i;
      break;
    }
  }

  if (index == -1) {
    if (isImprovement(c)) {
      wm_ = "Unable to handle case which would have been an improvement.";
    } else {
      wm_ = "Unable to handle case (would not have been an improvement).";
    }
    ext_warn(wm_, md_, cl_, vp_.lnw);
    return;
  }

  if (isBetter(c, population_[index].case_pointer)) {
    population_[index] = mating_pool_[index];
    if (isImprovement(c)) {
      updateTentativeBestCase(c);
      if (enable_logging_) {
        logger_->AddEntry(this);
      }
      if (vp_.vOPT >= 3) {
        im_ = "New best in generation " + num2str(iteration_);
        im_ += ": " + num2str(GetTentativeBestCase()->objf_value());
        ext_info(im_, "RGARDD", "Optimization", vp_.lnw);
      }
    }
  }
}

vector<GeneticAlgorithm::Chromosome>
RGARDD::selection(vector<Chromosome> population) {
  auto mating_pool = population_;
  int n_repl = floor(population_size_ * discard_parameter_);
  for (int i = population_size_-n_repl; i < population_size_; ++i) {
    mating_pool[i] = population_[i-population_size_+n_repl];
  }
  return sortPopulation(mating_pool);
}

vector<GeneticAlgorithm::Chromosome>
RGARDD::crossover(vector<Chromosome> mating_pool) {
  assert(mating_pool.size() == 2);

  auto p1 = mating_pool[0];
  auto p2 = mating_pool[1];
  auto r = random_doubles(gen_, 0, 1, n_vars_);
  Eigen::VectorXd dirs(n_vars_);

  for (int i = 0; i < n_vars_; ++i) {
    if (r[i] < 0.5) {
      dirs(i) = 0;
    } else {
      dirs(i) = p1.rea_vars(i) - p2.rea_vars(i);
    }
  }

  Chromosome o1 = Chromosome(p1.case_pointer);
  Chromosome o2 = Chromosome(p2.case_pointer);

  double s = abs(p1.ofv() - p2.ofv()) /
    (population_[0].ofv() - population_[population_size_-1].ofv());
  o1.rea_vars = p1.rea_vars + s * dirs;
  o2.rea_vars = p2.rea_vars + s * dirs;
  snap_to_bounds(o1);
  snap_to_bounds(o2);
  return vector<Chromosome>{o1, o2};
}

vector<GeneticAlgorithm::Chromosome>
RGARDD::mutate(vector<Chromosome> mating_pool) {
  assert(mating_pool.size() == 2);
  double s = pow(1.0 - (iteration_*1.0/max_generations_), decay_rate_);
  auto p1 = mating_pool[0];
  auto p2 = mating_pool[1];
  auto o1 = Chromosome(p1);
  auto o2 = Chromosome(p2);

  Eigen::VectorXd
    dir = random_doubles_eigen(gen_,
                               -mutation_strength_,
                               mutation_strength_,
                               n_vars_);

  o1.rea_vars = p1.rea_vars + s * dir.cwiseProduct(upper_bound_ - lower_bound_);
  o2.rea_vars = p2.rea_vars + s * dir.cwiseProduct(upper_bound_ - lower_bound_);

  // Snap to bound  constraints
  snap_to_bounds(o1);
  snap_to_bounds(o2);
  return vector<Chromosome> {o1, o2};
}

void RGARDD::snap_to_bounds(Chromosome &chrom) {
  for (int i = 0; i < chrom.rea_vars.size(); ++i) {
    if (chrom.rea_vars(i) < lower_bound_(i)) {
      chrom.rea_vars(i) = lower_bound_(i);
    } else if (chrom.rea_vars(i) > upper_bound_(i)) {
      chrom.rea_vars(i) = upper_bound_(i);
    }
  }
}

bool RGARDD::is_stagnant() {
  // Using the sums of the variable values in each chromosome
  vector<double> list_of_sums;
  for (auto chrom : population_) {
    list_of_sums.push_back(chrom.rea_vars.sum());
  }
  double stdev = calc_standard_deviation(list_of_sums);
  return stdev <= stagnation_limit_;
}

void RGARDD::repopulate() {
  for (int i = 0; i < population_size_ - 1; ++i) {
    auto new_case = generateRandomCase();
    population_[i+1] = Chromosome(new_case);
    case_handler_->AddNewCase(new_case);
  }
  mating_pool_ = population_;
}

Loggable::LogTarget RGARDD::ConfigurationSummary::GetLogTarget() {
  return LOG_SUMMARY;
}

map<string, string> RGARDD::ConfigurationSummary::GetState() {
  map<string, string> statemap;
  statemap["Name"] = "Genetic Algorithm (RGARDD)";
  statemap["Mode"] = opt_->mode_ == Settings::Optimizer::OptimizerMode::Maximize ? "Maximize" : "Minimize";
  statemap["Max Generations"] = boost::lexical_cast<string>(opt_->max_generations_);
  statemap["Max Evaluations"] = boost::lexical_cast<string>(opt_->max_evaluations_);
  statemap["Population Size"] = boost::lexical_cast<string>(opt_->population_size_);
  statemap["Discard Parameter"] = boost::lexical_cast<string>(opt_->discard_parameter_);
  statemap["Crossover Probability"] = boost::lexical_cast<string>(opt_->p_crossover_);
  statemap["Decay Rate"] = boost::lexical_cast<string>(opt_->decay_rate_);
  statemap["Mutation Strength"] = boost::lexical_cast<string>(opt_->mutation_strength_);

  string constraints_used;
  if (opt_->constraint_handler_ != nullptr) { // All actual cases
    for (auto cons : opt_->constraint_handler_->constraints()) {
      constraints_used += cons->name() + " ";
    }
  }
  statemap["Constraints"] = constraints_used;

  return statemap;
}

QUuid RGARDD::ConfigurationSummary::GetId() {
  return {}; // Null UUID
}

map<string, vector<double>> RGARDD::ConfigurationSummary::GetValues() {
  map<string, vector<double>> valmap;
  return valmap;
}

void RGARDD::print_state(string header) {
  std::stringstream ss;
  ss << header << "|";
  ss << "Iteration: " << iteration_ << "|" << "Current best case: ";
  ss << tentative_best_case_->id().toString().toStdString() << "|";
  ss << "              OFV: " << tentative_best_case_->objf_value() << "|";
  ext_info(ss.str(), md_, cl_, vp_.lnw);
}
}
}
