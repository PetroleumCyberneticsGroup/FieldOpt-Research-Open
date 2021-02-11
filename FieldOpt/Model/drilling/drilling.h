/********************************************************************
 Copyright (C) 2020-
 Thiago L. Silva <thiagolims@gmail.com>
 Mathias Bellout <chakibbb-pcg@gmail.com>


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

#ifndef FIELDOPT_DRILLING_H
#define FIELDOPT_DRILLING_H

#include <QString>
#include <QList>

#include "Model/properties/variable_property_container.h"
#include "Model/model_synchronization_object.h"
#include "drilling_schedule.h"
#include "Settings/model.h"
#include "Settings/settings.h"


#include "main_runner.h"

namespace Model {
namespace Drilling {

class Drilling {

 public:
  Drilling(Settings::Model *settings, Properties::VariablePropertyContainer *drilling_variables);

  int getCurrentStep() { return current_step_; }

  QString getWellName(){ return well_name_; }

  Settings::Optimizer* getLocalOptimizerSettings() { return local_optimizer_settings_;}
  Settings::Optimizer* getGlobalOptimizerSettings() { return global_optimizer_settings_;}

  Properties::VariablePropertyContainer* getVariables() { return drilling_variables_;}

  QMap<int, std::map<string, QHash<QUuid, double>>> getOptimalVariables() { optimal_variables_; }
  QMap<int, std::map<string, std::vector<double>>> getOptimalValues() { optimal_values_; }
  QMap<int, Runner::RuntimeSettings*> getRuntimeSettings() { return runtime_settings_; }

  DrillingSchedule* getDrillingSchedule() { return drilling_schedule_; }

  QString GetStatusString(int drilling_step) const;
  QString GetStatusString() const;
  QString GetStatusStringHeader() const;

  void setOptRuntimeSettings(int drilling_step, int argc, const char** argv);
  void setOptRuntimeSettings(int drilling_step, Runner::RuntimeSettings* rts);
  void maintainRuntimeSettings(int drilling_step);

  void modelUpdate(int drilling_step);
  void runOptimization(int drilling_step);
  void createLogFile(int drilling_step);

  double getBestObjective() { return best_objective_;}



 private:
  int current_step_;
  int current_model_;

  bool skip_optimization_;
  double best_objective_;

  string original_output_dir_;
  QString dw_log_path_;

  Properties::VariablePropertyContainer *drilling_variables_;
  Settings::Model  *settings_;

  DrillingSchedule *drilling_schedule_;
  QString well_name_;

  Settings::Optimizer* local_optimizer_settings_;
  Settings::Optimizer* global_optimizer_settings_;

  QMap<int, QHash<QUuid, double>> optimal_variables_;  //!< Optimal variables per drilling step
  QMap<int, QHash<QUuid, int>> optimal_int_variables_;  //!< Optimal variables per drilling step
  QMap<int, QHash<QUuid, bool>> optimal_bin_variables_;  //!< Optimal variables per drilling step

  QMap<int, std::map<string, std::vector<double>>>  optimal_values_;     //!< Optimal values per drilling step
  QMap<int, Runner::RuntimeSettings*> runtime_settings_;                    //!< Optimization runtime settings per drilling step

  Optimization::Case* best_case_;
  Optimization::Case* base_case_;
  ModelSynchronizationObject* mso_;

  void setWellOptimalVariables(const QHash<QUuid, bool>& opt_bin_var,  const QHash<QUuid, int>& opt_int_var, const QHash<QUuid, double>& opt_var, int drilling_step);
  void setWellOptimizationValues(const std::map<string, std::vector<double>>& opt_val, int drilling_step);


  void printIteration(int drilling_step, double model_deviation, double obj_improvement, bool skip_opt);

};

}
}
#endif //FIELDOPT_DRILLING_H
