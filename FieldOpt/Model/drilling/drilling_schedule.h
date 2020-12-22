/********************************************************************
 Copyright (C) 2020-
 Mathias Bellout <chakibbb-pcg@gmail.com>
 Thiago L. Silva <thiagolims@gmail.com>

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
*********************************************************************/

#ifndef FIELDOPT_DRILLING_SCHEDULE_H
#define FIELDOPT_DRILLING_SCHEDULE_H

#include "Model/properties/variable_property_container.h"
#include "drilling_schedule.h"
#include "Settings/model.h"
#include "Settings/optimizer.h"
#include <iostream>

namespace Model {
namespace Drilling {

class DrillingSchedule {
 private:
  void assignDrillingPoints(QMap<int, QList<Settings::Model::Drilling::DrillingPoint>> drilling_points_settings);
  void assignOptimizationTriggers(QMap<int, Settings::Model::Drilling::OptimizationTrigger> opt_triggers);

 public:
  struct DrillingPoint {
    DrillingPoint() {}
    double x, y, z;
    bool is_variable = false;
  };

  struct OptimizationTrigger {
    OptimizationTrigger() {}
    double min_model_deviation;
    double max_model_deviation;
    double max_objective_improvement;
  };

  enum DrillingOperation : int {StartDrilling=1, Drilling=2, PullingOutOfHole=3};
  enum ModelType: int {TrueModel=1, Surrogate=2};
  enum Execution: int {Serial=1, Parallel=2};

  DrillingSchedule(Settings::Model *settings, Properties::VariablePropertyContainer *variables);

  QList<int> getSteps() { return drilling_steps_; }
  QMap<int, double> getTimeSteps() { return time_steps_; }
  QMap<int, QList<DrillingPoint>> getDrillingPoints() { return drilling_points_; }
  QMap<int, DrillingOperation> getDrillingOperations() { return drilling_operations_; }
  QMap<int, ModelType> getModelTypes() { return model_types_; }
  QMap<int, Execution> getExecutionModes() { return execution_modes_; }
  QMap<int, Settings::Optimizer*> getOptimizerSettings() { return optimizer_settings_;}
  QMap<int, OptimizationTrigger> getOptimizationTriggers() { return optimization_triggers_; }

  QMap<int, bool> isVariableDrillingPoints() { return is_variable_drilling_points_; }
  QMap<int, bool> isVariableCompletions() { return is_variable_completions_; }
  QMap<int, bool> isModelUpdates() { return is_model_updates_;}

 private:
  QList<int> drilling_steps_;
  QMap<int, double> time_steps_; //!< Indexed by the drilling steps
  QMap<int, QList<DrillingPoint>> drilling_points_; //!< Indexed by the drilling steps
  QMap<int, OptimizationTrigger> optimization_triggers_;  //!< Indexed by the drilling steps

  QMap<int, DrillingOperation> drilling_operations_; //!< Indexed by the drilling steps
  QMap<int, ModelType> model_types_;                 //!< Indexed by the drilling steps
  QMap<int, Execution> execution_modes_;             //!< Indexed by the drilling steps
  QMap<int, Settings::Optimizer*> optimizer_settings_;//!< Indexed by the drilling steps
  QMap<int, bool> is_variable_drilling_points_;      //!< Indexed by the drilling steps
  QMap<int, bool> is_variable_completions_;          //!< Indexed by the drilling steps
  QMap<int, bool> is_model_updates_;                 //!< Indexed by the drilling steps

  void printDrillingPoints();
};

}
}
#endif //FIELDOPT_DRILLING_SCHEDULE_H
