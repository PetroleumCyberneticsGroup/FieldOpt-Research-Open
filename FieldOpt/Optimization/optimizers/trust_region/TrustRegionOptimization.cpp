/******************************************************************************
   Created by thiagols on 26.11.18
   Copyright (C) 2018 Thiago Lima Silva<thiagolims@gmail.com>

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
#include <Utilities/verbosity.h>
#include "Utilities/printer.hpp"
#include "Utilities/stringhelpers.hpp"
#include "TrustRegionOptimization.h"

namespace Optimization {
namespace Optimizers {

TrustRegionOptimization::TrustRegionOptimization(Settings::Optimizer *settings,
         Case *base_case,
         Model::Properties::VariablePropertyContainer *variables,
         Reservoir::Grid::Grid *grid,
         Logger *logger,
         CaseHandler *case_handler,
         Constraints::ConstraintHandler *constraint_handler
) : Optimizer(settings, base_case, variables, grid, logger, case_handler, constraint_handler) {
    if (constraint_handler_->HasBoundaryConstraints()) {
        lb_ = constraint_handler_->GetLowerBounds(base_case->GetRealVarIdVector());
        ub_ = constraint_handler_->GetUpperBounds(base_case->GetRealVarIdVector());
    }
    else {
        lb_.resize(base_case->GetRealVarIdVector().size());
        ub_.resize(base_case->GetRealVarIdVector().size());
        lb_.fill(settings->parameters().lower_bound);
        ub_.fill(settings->parameters().upper_bound);
    }

    settings_ = settings;
    variables_ = variables;
    base_case_ = base_case;

    if (enable_logging_) { // Log base case
        logger_->AddEntry(this);
    }

    TrustRegionOptimization::computeInitialPoints();

    if (enable_logging_) {
        logger_->AddEntry(new ConfigurationSummary(this));
    }
}
Optimization::Optimizer::TerminationCondition TrustRegionOptimization::IsFinished() {
    TerminationCondition tc = NOT_FINISHED;
    if (case_handler_->CasesBeingEvaluated().size() > 0)
        return tc;
    if (evaluated_cases_ > max_evaluations_)
        tc = MAX_EVALS_REACHED;
    if (tc != NOT_FINISHED) {
        if (enable_logging_) {
            logger_->AddEntry(this);
        }
    }
    return tc;
}
void TrustRegionOptimization::handleEvaluatedCase(Case *c) {
    if (!tr_model_ ) {
        if (initial_points_.cols()  < 1) {
            Printer::ext_warn("Insufficient number of points to create the Trust Region model.", "Optimization", "TrustRegionOptimization");
            throw std::runtime_error("Failed to initialize Trust Region Optimizer.");
        } else if (initial_points_.cols()  == 1 ) { //!<just evaluated the second point, so we are ready to build the quadratic model for the TR
            initial_points_.col(1) = c->GetRealVarVector();
            initial_fvalues_(1) = c->objective_function_value();
        }
        tr_model_ = new TrustRegionModel(initial_points_, initial_fvalues_, settings_);  //!<creates the initial trust region
    }
    //TODO: improve trust region with the new point

    if (isImprovement(c)) {
        updateTentativeBestCase(c);
        Printer::ext_info("Found new tentative best case: " + Printer::num2str(c->objective_function_value()), "Optimization", "EGO");
    }
}
void TrustRegionOptimization::iterate() {
    // TODO: implement the main iteration for the TR optimization algorithm

    if (enable_logging_) {
        logger_->AddEntry(this);
    }
}

void TrustRegionOptimization::computeInitialPoints() {
    int n_cont_vars = variables_->ContinousVariableSize();
    if (settings_->parameters().tr_init_guesses == -1) { //!<only one initial guess was provided, thus the algorithm finds another point.
        if (settings_->parameters().tr_init_sampling_method == "Random") {
            auto rng = get_random_generator(settings_->parameters().rng_seed);

            VectorXd pos = VectorXd::Zero(n_cont_vars);
            for (int i = 0; i < n_cont_vars; ++i) {
                pos(i) = random_double(rng, lb_(i), ub_(i));
            }
            Case * init_case = new Case(base_case_);
            init_case->SetRealVarValues(pos);
            case_handler_->AddNewCase(init_case);
            initial_points_.setZero(n_cont_vars, 2); //!<only two points to start the quadratic model for the trust region>

            initial_fvalues_.setZero(2);
            initial_fvalues_(0) = base_case_->objective_function_value() ;
        } else {
            //TODO: implement other sampling method such as Uniform.
        }
    } else {
        //TODO: get other initial points from the parameters
    }
}

Loggable::LogTarget TrustRegionOptimization::ConfigurationSummary::GetLogTarget() {
    return LOG_SUMMARY;
}
QUuid TrustRegionOptimization::ConfigurationSummary::GetId() {
    return QUuid(); // Null UUID
}

}
}