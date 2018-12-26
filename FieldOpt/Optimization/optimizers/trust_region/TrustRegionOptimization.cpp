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

TrustRegionOptimization::TrustRegionOptimization(
        Settings::Optimizer *settings,
        Case *base_case,
        Model::Properties::VariablePropertyContainer *variables,
        Reservoir::Grid::Grid *grid,
        Logger *logger,
        CaseHandler *case_handler,
        Constraints::ConstraintHandler *constraint_handler
) : Optimizer(settings, base_case, variables, grid, logger, case_handler, constraint_handler) {

    settings_ = settings;
    variables_ = variables;
    base_case_ = base_case;

    setLowerUpperBounds();

    if (enable_logging_) { // Log base case
        logger_->AddEntry(this);
    }

    n_initial_points_ = 0;
    computeInitialPoints();

    if (enable_logging_) {
        logger_->AddEntry(new ConfigurationSummary(this));
    }
}

void TrustRegionOptimization::setLowerUpperBounds() {
    
  if (constraint_handler_->HasBoundaryConstraints()) {

    // Use constraint_handler lower/upper bounds
    lb_ = constraint_handler_->GetLowerBounds(
        base_case_->GetRealVarIdVector());
    ub_ = constraint_handler_->GetUpperBounds(
        base_case_->GetRealVarIdVector());

  } else {

    // Use lower/upper bounds specified in tr-params
    if (settings_->parameters().tr_lower_bound
        && settings_->parameters().tr_upper_bound) {

      lb_.resize(base_case_->GetRealVarIdVector().size());
      ub_.resize(base_case_->GetRealVarIdVector().size());
      lb_.fill(settings_->parameters().tr_lower_bound);
      ub_.fill(settings_->parameters().tr_upper_bound);

    } else {

      Printer::ext_warn(
          "Lower/upper bounds for DF-TR algorithm not specified.",
          "Optimization", "TrustRegionOptimization");
      throw std::runtime_error(
          "Lower/upper bounds for DF-TR algorithm not specified.");

    }
  }
}

// ===================================================================
Optimization::Optimizer::TerminationCondition
TrustRegionOptimization::IsFinished() {
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
    IOFormat frmt(3, 0, " ", "\n", "             [", "]");

    if ((n_initial_points_  < 1) && (case_handler_->CasesBeingEvaluated().size() == 0)) {
        Printer::ext_warn("Insufficient number of points to create the Trust Region model.", "Optimization", "TrustRegionOptimization");
        throw std::runtime_error("Failed to initialize Trust Region Optimizer.");
    } else if (n_initial_points_ == 1 ) { //!<just evaluated the second point, so we are ready to build the quadratic model for the TR
        initial_points_.col(1) = c->GetRealVarVector();
        initial_fvalues_(1) = c->objective_function_value();
        n_initial_points_++;

        //!<dbg>
        //cout << "[          ] initial_points_" << endl;
        //cout << initial_points_.format(frmt) << endl;

        //cout << "[          ] initial_fvalues_" << endl;
        //cout << initial_fvalues_.format(frmt) << endl;

        tr_model_ = new TrustRegionModel(initial_points_, initial_fvalues_,
                                         lb_, ub_, settings_);  //!<creates the initial trust region

        //TODO: improve trust region with the new point

        if (isImprovement(c)) {
            updateTentativeBestCase(c);
            Printer::ext_info("Found new tentative best case: " + Printer::num2str(c->objective_function_value()),
                              "Optimization", "Trust Region");
        }
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
    auto initial_point = base_case_->GetRealVarVector();
    //!<dbg>
    //IOFormat frmt(3, 0, " ", "\n", "             [", "]");

    //!<Find another point since only one initial guess provided
    if (settings_->parameters().tr_init_guesses == -1) {
        if (settings_->parameters().tr_init_sampling_method == "Random") {

            //!<Find (random) 2nd initial point>
            auto rng = get_random_generator(settings_->parameters().rng_seed);

            VectorXd second_point(n_cont_vars, 1);
            second_point.setZero(n_cont_vars);
            for (int i = 0; i < n_cont_vars; ++i) {
                second_point(i) = random_double(rng, 0, 1);
            }

            //!<dbg>
            //cout << "[          ] second_point 1" << endl;
            //cout << second_point.format(frmt) << endl;

            //cout << "[          ] initial_point 1" << endl;
            //cout << initial_point.format(frmt) << endl;

            //cout << "lp_infinity:" << endl;
            //cout << second_point.lpNorm<Infinity>() << endl;
            //cout << "tr_pivot_threshold" << endl;
            //cout << settings_->parameters().tr_pivot_threshold << endl;


            while (second_point.lpNorm<Infinity>()  < settings_->parameters().tr_pivot_threshold) { //!<Second point must not be too close>
                second_point << 2*second_point.array();
            }


            second_point << (second_point.array() - 0.5);
            second_point *= settings_->parameters().tr_initial_radius;
            second_point << initial_point.array() + second_point.array();
            if(lb_.size()>0 && ub_.size() >0) {
                projectToBounds(&second_point);
            }

            //!<dbg>
            //cout << "[          ] second_point" << endl;
            //cout << second_point.format(frmt) << endl;


            //!<Compute case corresponding to 2nd init point>
            Case *second_case = new Case(base_case_);
            second_case->SetRealVarValues(second_point);
            case_handler_->AddNewCase(second_case);

            //!<Establish initial point matrix (2 cols since only
            //!< two points are needed to build quad model for trust region)>
            initial_points_.setZero(n_cont_vars, 2);
            initial_points_.col(0) = initial_point;

            //!<Establish initial feval matrix>
            initial_fvalues_.setZero(2);
            initial_fvalues_(0) = base_case_->objective_function_value();
            n_initial_points_++;

            //!<dbg>
            //cout << "[          ] initial_points_" << endl;
            //cout << initial_points_.format(frmt) << endl;

            //cout << "[          ] initial_fvalues_" << endl;
            //cout << initial_fvalues_.format(frmt) << endl;

        } else {
            //TODO: implement other sampling method such as Uniform.
        }

    } else {
        //TODO: get other initial points from the parameters
    }
}

void TrustRegionOptimization::projectToBounds(VectorXd *point) {
    if (ub_.size() > 0) {
        *point = ub_.cwiseMin(*point);
    }
    if (lb_.size() > 0) {
        *point = lb_.cwiseMax(*point);
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
