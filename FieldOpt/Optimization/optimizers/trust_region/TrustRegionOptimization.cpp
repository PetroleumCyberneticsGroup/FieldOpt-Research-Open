/*********************************************************************
 Created by thiagols on 26.11.18
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

    n_initial_points_ = 1;

    // Construct shell of TRModel, does not initialize
    // model, i.e., is_model_initialized_ = false
    tr_model_ = new TrustRegionModel(lb_, ub_, base_case_, settings_);

    // [1] Find 2nd point; make and add cases corresponding
    // to 1st & 2nd points to initialization_cases_ list
    // [2] Fields initial_points_ and initial_fvalues_ are first
    // established here using fields from base_case (thus, these
    // need to be removed from constructor above, where they play
    // no part)
    computeInitialPoints();

    // Log configuration
    if (enable_logging_) {
        logger_->AddEntry(new ConfigurationSummary(this));
    }
}

void TrustRegionOptimization::iterate() {

    if (iteration_ == 0) {

        if (!tr_model_->areInitPointsComputed()
            && !tr_model_->isInitialized()) {
            cout << "init_cases -> case_handler_ (return)" << endl;
            auto init_cases = tr_model_->getInitializationCases();
            case_handler_->AddNewCases(init_cases);
            cout << "# of cases in case_handler: " << init_cases.size() << endl;
            return;
        }

        if (tr_model_->areInitPointsComputed()
            && !tr_model_->areImprovementPointsComputed()
            && !tr_model_->areReplacementPointsComputed()
            && !tr_model_->isInitialized()) {
            cout << "Initializing TRModel" << endl;

            tr_model_->setModelChanged(tr_model_->rebuildModel());
            tr_model_->moveToBestPoint();
            tr_model_->computePolynomialModels();

            if (tr_model_->hasOnlyOnePoint()) {
                cout << "Improve TRModel" << endl;

                int exit_flag = tr_model_->ensureImprovement();
                auto improvement_cases = tr_model_->getImprovementCases(); //<improve model>
                auto replacement_cases = tr_model_->getReplacementCases(); //<replace point>
                // Might be 0 if point_found in improveModelNfp is false. We also need to handle exit_flag=3 or exit_flag=4

                if (tr_model_->isImprovementNeeded()
                    && !improvement_cases.size() == 0) {
                    cout << "impr_cases -> case_handler_ (return)" << endl;
                    case_handler_->AddNewCases(improvement_cases);
                    return;
                }

                if (tr_model_->isReplacementNeeded() && !replacement_cases.size() ==0) {
                    cout << "repl_cases -> case_handler_ (return)" << endl;
                    case_handler_->AddNewCases(replacement_cases);
                    return;
                }
            } else {
                cout << "Initialized TRModel [1]" << endl;
                tr_model_->setIsInitialized(true);
            }
        } // end-if: Initializing TRModel

        if (tr_model_->areInitPointsComputed()
            && tr_model_->isImprovementNeeded()
            && tr_model_->areImprovementPointsComputed()
            && !tr_model_->isInitialized()) {

          cout << "Continue model improvement process (nfp model)" << endl;
          int exit_flag = tr_model_->ensureImprovement();

          if (!tr_model_->isImprovementNeeded()) {
              cout << "Improvement process successful (nfp model)" << endl;
              cout << "Initialized TRModel [2]" << endl;
              tr_model_->setIsInitialized(true);
          }
        }  // end-if: finish model improvement process

        if (tr_model_->areInitPointsComputed()
            && tr_model_->isReplacementNeeded()
            && tr_model_->areReplacementPointsComputed()
            && !tr_model_->isInitialized()) {

          cout << "Continue model improvement process (replacement of points)" << endl;
          int exit_flag = tr_model_->ensureImprovement();

          if (!tr_model_->isReplacementNeeded()) {
            cout << "Improvement process successful (replacement of points)" << endl;
            cout << "Initialized TRModel [2]" << endl;
            tr_model_->setIsInitialized(true);
          }
        }  // end-if: finish model improvement process
    }  // end-if: (iteration_ == 0)

    if(tr_model_->isInitialized()) {

        cout << "TRMmodel fully initialized" << endl;
        cout << "Starting first iteration" << endl;
        // termination crit.: check for # of iterations
        // termination crit.: check for radius size

        // check if model is lambda-poised
        // check criticality step

        // compute step: solve_tr_subproblem
        iteration_++;
        return;
    }  // end-if: (tr_model_->isInitialized)

    // Save for later
    //    Printer::ext_warn(
    //            "Insufficient # of points to create the Trust Region model.",
    //            "Optimization", "TrustRegionOptimization");
    //    throw std::runtime_error(
    //            "Failed to initialize Trust Region Optimizer.");

}

void TrustRegionOptimization::handleEvaluatedCase(Case *c) {

    if (iteration_ == 0) {

        // Collect init points
        if (!tr_model_->areInitPointsComputed()
           && !tr_model_->isInitialized()) {

            // Collect initialization case
            tr_model_->addTempInitCase(c);

            // All initialization cases have been evaluated
            if (case_handler_->QueuedCases().size() == 0
                && case_handler_->CasesBeingEvaluated().size() == 0) {

                    cout << "Return evaluated init_cases_ (sz:"
                         << tr_model_->getSizeInitCases() << ")" << endl;
                    tr_model_->submitTempInitCases();
                    tr_model_->setAreInitPointsComputed(true);
            }
            return;

        } else if (tr_model_->areInitPointsComputed()
                   && !tr_model_->isInitialized()) {

            // Collect initialization case
            tr_model_->addTempImprCase(c);

            // All improvement cases have been evaluated
            if (case_handler_->QueuedCases().size() == 0
                && case_handler_->CasesBeingEvaluated().size() == 0) {

                cout << "Return evaluated impr_cases_ (sz:"
                     << tr_model_->getSizeImprCases() << ")" << endl;

              cout << "Return evaluated repl_cases_ (sz:"
                   << tr_model_->getSizeReplCases() << ")" << endl;
                tr_model_->submitTempImprCases();
                tr_model_->submitTempReplCases();
                tr_model_->setAreImprPointsComputed(true);
            }
            return;
        }

    } else { // Initialization is done; handling "normal" case
        // do normal stuff
    }
}

void TrustRegionOptimization::computeInitialPoints() {

    tr_model_->setXDim(variables_->ContinousVariableSize());
    int n_cont_vars = tr_model_->getXDim();

    auto initial_point = base_case_->GetRealVarVector();

    cout << "Adding base case to Case List" << endl;
    tr_model_->addInitializationCase(base_case_);

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

            //!<Second point must not be too close>
            while (second_point.lpNorm<Infinity>()
                    < settings_->parameters().tr_pivot_threshold) {
                second_point << 2*second_point.array();
            }

            second_point << (second_point.array() - 0.5);
            second_point *= settings_->parameters().tr_initial_radius;
            second_point << initial_point.array() + second_point.array();
            if(lb_.size()>0 && ub_.size() >0) {
                projectToBounds(&second_point);
            }

            //!<Append case corresponding to 2nd init point>
            Case *second_case = new Case(base_case_);
            second_case->SetRealVarValues(second_point);

            // Hack: Comment this line to test improveModelNfp()
            // -> currently, this leads to a crash
            // This will be done later in the test code
            tr_model_->addInitializationCase(second_case);
            cout << "Adding 2nd point to Case List" << endl;

            n_initial_points_++;

        } else {
            //TODO: implement other sampling method such as Uniform.
        }

    } else {
        //TODO: get other initial points from the parameters
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

Optimization::Optimizer::TerminationCondition
TrustRegionOptimization::IsFinished() {

    TerminationCondition tc = NOT_FINISHED;
    if (case_handler_->CasesBeingEvaluated().size() > 0) {
        return tc;
    }

    if (evaluated_cases_ > max_evaluations_) {
          tc = MAX_EVALS_REACHED;
    }

    if (tc != NOT_FINISHED) {
        if (enable_logging_) {
            logger_->AddEntry(this);
        }
    }
    return tc;
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
