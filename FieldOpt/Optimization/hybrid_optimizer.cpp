/******************************************************************************
   Copyright (C) 2015-2018 Einar J.M. Baumann <einar.baumann@gmail.com>

   This file is part of the FieldOpt project.

   FieldOpt is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   FieldOpt is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FieldOpt.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#include <Optimization/optimizers/compass_search.h>
#include <Optimization/optimizers/APPS.h>
#include <Optimization/optimizers/RGARDD.h>
#include <Optimization/optimizers/bayesian_optimization/EGO.h>
#include <Utilities/printer.hpp>
#include <Utilities/verbosity.h>
#include "hybrid_optimizer.h"

namespace Optimization {

HybridOptimizer::HybridOptimizer(Settings::Optimizer *settings,
         Case *base_case,
         Model::Properties::VariablePropertyContainer *variables,
         Reservoir::Grid::Grid *grid,
         Logger *logger
)
    : Optimizer(settings, base_case, variables, grid, logger)
{
    variables_ = variables;
    grid_ = grid;
    iteration_ = 0;

    assert(settings->HybridComponents().size() == 2);
    primary_settings_ = new Settings::Optimizer(settings->HybridComponents()[0]);
    secondary_settings_ = new Settings::Optimizer(settings->HybridComponents()[1]);

    initializeComponent(0);
    active_component_ = 0;

    component_improvement_found_ = true;

    max_hybrid_iterations_ = settings->parameters().hybrid_max_iterations;
    if (settings->parameters().hybrid_switch_mode == "OnConvergence") {
        switch_mode_ = HybridSwitchMode::ON_CONVERGENCE;
    }
    else {
        throw std::runtime_error("Hybrid optimizer switch mode not recognized.");
    }
    if (settings->parameters().hybrid_termination_condition == "NoImprovement") {
        hybrid_termination_condition_ = HybridTerminationCondition::NO_IMPROVEMENT;
    }
    else {
        throw std::runtime_error("Hybrid optimizer termination condition not recognized.");
    }
}

Optimizer::TerminationCondition HybridOptimizer::IsFinished() {
    if (case_handler_->CasesBeingEvaluated().size() > 0) {
        return TerminationCondition::NOT_FINISHED;
    }
    else if (hybrid_termination_condition_ == HybridTerminationCondition::NO_IMPROVEMENT
             && component_improvement_found_ == false && iteration_ > 1) {
        Printer::ext_info("No improvement found in previous component run. Terminating.", "Optimization", "HybridOptimizer");
        if (enable_logging_) {
            logger_->AddEntry(this);
        }
        return TerminationCondition::MINIMUM_STEP_LENGTH_REACHED;
    }
    else if (iteration_ < max_hybrid_iterations_) {
        return TerminationCondition::NOT_FINISHED;
    }
    else {
        if (active_component_ == 0) {
            return TerminationCondition::NOT_FINISHED;
        }
        else {
            if (secondary_->IsFinished() != NOT_FINISHED) {
                if (enable_logging_) {
                    logger_->AddEntry(this);
                }
            }
            return secondary_->IsFinished();
        }
    }
}
void HybridOptimizer::handleEvaluatedCase(Case *c) {
    if (active_component_ == 0) {
        if (isImprovement(c)) {
            if (VERB_OPT >= 1) {
                std::stringstream ss;
                ss << "Found better case. Passing to primary." << "|";
                ss << " ID: " << c->id().toString().toStdString() << "|";
                ss << "OFV: " << c->objective_function_value();
                Printer::ext_info(ss.str(), "Optimization", "HybridOptimizer");
            }
            updateTentativeBestCase(c);
        }
        primary_->handleEvaluatedCase(c);
    }
    else {
        if (isImprovement(c)) {
            if (VERB_OPT >= 1) {
                std::stringstream ss;
                ss << "Found better case. Passing to secondary." << "|";
                ss << " ID: " << c->id().toString().toStdString() << "|";
                ss << "OFV: " << c->objective_function_value();
                Printer::ext_info(ss.str(), "Optimization", "HybridOptimizer");
            }
            updateTentativeBestCase(c);
        }
        secondary_->handleEvaluatedCase(c);
    }

}
void HybridOptimizer::iterate() {
    if (enable_logging_) {
        logger_->AddEntry(this);
    }
    if (active_component_ == 0) { // Primary is active.
        if (primary_->IsFinished() == TerminationCondition::NOT_FINISHED) { // Primary is not finished.
            if (VERB_OPT >= 1) { Printer::ext_info("Iterating with primary.", "Optimization", "HybridOptimizer"); }
            primary_->iterate();
        }
        else { // Primary is finished. Switch to secondary.
            if (VERB_OPT >= 1) { Printer::ext_info("Primary component converged. Switching to secondary.", "Optimization", "HybridOptimizer"); }
            iteration_++;
            primary_best_case_ = tentative_best_case_;
            initializeComponent(1);
            active_component_ = 1;
            if (case_handler_->QueuedCases().size() == 0) { // Iterate if the constructor does not generate cases
                secondary_->iterate();
            }

            if (iteration_ > 1 && isBetter(primary_best_case_, secondary_best_case_) == false) {
                // If primary didn't find an improvement over previous secondary
                component_improvement_found_ = false;
            }
        }
    }
    else { // Secondary is active
        if (secondary_->IsFinished() == TerminationCondition::NOT_FINISHED) { // Secondary is not finished.
            if (VERB_OPT >= 1) { Printer::ext_info("Iterating with secondary.", "Optimization", "HybridOptimizer"); }
            secondary_->iterate();
        }
        else { // Secondary is finished. Switch back to primary.
            if (VERB_OPT >= 1) { Printer::ext_info("Secondary component converged. Switching to primary.", "Optimization", "HybridOptimizer"); }
            iteration_++;
            secondary_best_case_ = tentative_best_case_;
            initializeComponent(0);
            active_component_ = 0;
            if (case_handler_->QueuedCases().size() == 0) { // Iterate if the constructor does not generate cases
                primary_->iterate();
            }

            if (isBetter(secondary_best_case_, primary_best_case_) == false) {
                // If secondary didn't find an improvement over previous primary
                component_improvement_found_ = false;
            }
        }
    }
}
void HybridOptimizer::initializeComponent(int component=0) {
    assert(component == 0 || component == 1);

    Settings::Optimizer *opt_settings;
    std::string compstr = "";
    if (component == 0) {
        opt_settings = primary_settings_;
        compstr = "primary";
    }
    else {
        opt_settings = secondary_settings_;
        compstr = "secondary";
    }
    opt_settings->SetRngSeed(opt_settings->parameters().rng_seed + iteration_ * 7);
    opt_settings->set_mode(mode_);

    Optimizer *opt;

    switch (opt_settings->type()) {
        case Settings::Optimizer::OptimizerType::Compass:
            Printer::ext_info("Using Compass Search as " + compstr + " component in hybrid algorithm.", "Optimization", "HybridOptimizer");
            opt = new Optimizers::CompassSearch(
                opt_settings, tentative_best_case_, variables_, grid_, logger_,
                case_handler_, constraint_handler_
            );
            break;
        case Settings::Optimizer::OptimizerType::APPS:
            Printer::ext_info("Using APPS as " + compstr + " component in hybrid algorithm.", "Optimization", "HybridOptimizer");
            opt = new Optimizers::APPS(
                opt_settings, tentative_best_case_, variables_, grid_, logger_,
                case_handler_, constraint_handler_
            );
            break;
        case Settings::Optimizer::OptimizerType::GeneticAlgorithm:
            Printer::ext_info("Using Genetic Algorithm as " + compstr + " component in hybrid algorithm.|RNG Seed: "
                               + Printer::num2str(opt_settings->parameters().rng_seed), "Optimization", "HybridOptimizer");
            opt = new Optimizers::RGARDD(
                opt_settings, tentative_best_case_, variables_, grid_, logger_,
                case_handler_, constraint_handler_
            );
            break;
        case Settings::Optimizer::OptimizerType::EGO:
            Printer::ext_info("Using EGO as " + compstr + " component in hybrid algorithm.", "Optimization", "HybridOptimizer");
            opt = new Optimizers::BayesianOptimization::EGO(
                opt_settings, tentative_best_case_, variables_, grid_, logger_,
                case_handler_, constraint_handler_
            );
            break;
        default:
            throw std::runtime_error("Unable to initialize hybrid optimizer: algorithm not recognized.");
    }
    opt->DisableLogging();
    if (component == 0)
        primary_ = opt;
    else
        secondary_ = opt;
}
}
