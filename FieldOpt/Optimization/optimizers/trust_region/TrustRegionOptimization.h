/***********************************************************
Created by thiagols on 27.11.18
Copyright (C) 2018
Thiago Lima Silva <thiagolims@gmail.com>

Modified 2018-2019 Mathias Bellout
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

#ifndef FIELDOPT_TRUSTREGION_H
#define FIELDOPT_TRUSTREGION_H

#include "Optimization/optimizer.h"
#include "Optimization/optimizers/trust_region/TrustRegionModel.h"
#include <Eigen/Core>

using namespace Eigen;

namespace Optimization {
namespace Optimizers {

// =========================================================
/*!
 * @brief Implementation of Trust Region (TR) algorithm
 * for derivative-free optimization
 *
 * \todo Derivative-free Optimization algorithm
 * \todo TR model
 * \todo TR sub-problem optimization
 * \todo Implement tests
 */
class TrustRegionOptimization : public Optimizer {

 public:

  TerminationCondition IsFinished() override;

  TrustRegionOptimization(
    Settings::Optimizer *settings,
    Case *base_case,
    Model::Properties::VarPropContainer *variables,
    Reservoir::Grid::Grid *grid,
    Logger *logger,
    CaseHandler *case_handler = 0,
    Constraints::ConstraintHandler *constraint_handler=0
  );

  TrustRegionModel *getTrustRegionModel() { return tr_model_; };

  int GetNumInitPoints() { return n_initial_points_; };

  int getNumIterations() { return iteration_; };

  bool ensureImprovementPostProcessing();

 protected:
  void handleEvaluatedCase(Case *c) override;
  void iterate() override;

 private:

  // -------------------------------------------------------
  // TR constructs
  // Matrix<double, Dynamic, Dynamic> initial_points_;
  // RowVectorXd initial_fvalues_;
  VectorXd lb_, ub_; //!< Upper and lower bounds
  TrustRegionModel *tr_model_;
  int n_initial_points_{};

  double fval_current_{}; //!<Function value of current point>
  VectorXd x_current_; //!<Current point>
  double rho_; //!<Agreement factor>
  double ared_{};//!<Actual reduction>
  double sum_rho_; //!<Cumulative rho>
  double sum_rho_sqr_; //!<Cumulative squared rho>
  double delay_reduction_; //!<Delay reduction>
  int mchange_flag_{};
  double gamma_dec_;

  VectorXd trial_point_;
  VectorXd trial_step_;
  double fval_trial_{};
  double predicted_red_{};
  double criticality_init_radius_;

  bool iteration_model_fl_ = false;
  bool criticality_step_performed_ = false;
  bool improve_model_ = false;
  bool criticality_step_execution_ongoing_ = false;

  QString tr_log_path_;
  Settings::VerbParams vp_;

  // -------------------------------------------------------
  // FO constructs
  Settings::Optimizer *settings_;
  Model::Properties::VarPropContainer *variables_;
  Case *base_case_;

  void computeInitialPoints();
  void updateRadius();

  void setLowerUpperBounds();
  void projectToBounds(VectorXd *point);
  void printIteration(double fval_current);
  // \todo make new version of printIteration function
  //  such that function rints out to file conditionally
  void printIteration(double fval_current, string ifile);
  void createLogFile();

  // -------------------------------------------------------
  class ConfigurationSummary : public Loggable {
   public:
    ConfigurationSummary(TrustRegionOptimization *opt) { opt_ = opt; }
    LogTarget GetLogTarget() override;
    QUuid GetId() override;
   private:
    TrustRegionOptimization *opt_;
  };
};

}  // Optimizers
}  // Optimization
#endif //FIELDOPT_TRUSTREGIONOPTIMIZATION_H
