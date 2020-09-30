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
***************************************************************** ****/


#include "drilling.h"
#include "paths.h"
#include "Optimization/case.h"

namespace Model {
namespace Drilling {

Drilling::Drilling(Settings::Model *settings, Properties::VariablePropertyContainer *drilling_variables) {
  drilling_variables_  = drilling_variables;
  settings_   = settings;


  well_name_ = settings->drilling().well_name;
  drilling_schedule_ = new DrillingSchedule(settings_, drilling_variables_);

  current_step_ = 0;
  current_model_ = 0;

  best_case_ = 0;
  mso_ = 0;
};

void Drilling::setWellOptimalVariables(const QHash<QUuid, bool>& opt_bin_var,  const QHash<QUuid, int>& opt_int_var, const QHash<QUuid, double>& opt_var, int drilling_step) {
  optimal_int_variables_.insert(drilling_step, opt_int_var);
  optimal_bin_variables_.insert(drilling_step, opt_bin_var);
  optimal_variables_.insert(drilling_step, opt_var);
}

void Drilling::setWellOptimizationValues(const std::map<string, std::vector<double>>& opt_val, int drilling_step) {
  optimal_values_.insert(drilling_step, opt_val);
}

QString Drilling::GetStatusStringHeader() const
{
  return QString("%1,%2,%3,%4,%5,%6\n")
      .arg("DrillingStep")
      .arg("NumberIterations")
      .arg("TerminationCriteria")
      .arg("TentativeBestCaseOFValue");
}


QString Drilling::GetStatusString(int drilling_step) const
{
  return QString("%1,%2,%3,%4 \n")
      .arg(drilling_step)
      .arg(optimal_values_.value(drilling_step).at("IterNr")[0])
      .arg(optimal_values_.value(drilling_step).at("TermCond")[0])
      .arg(optimal_values_.value(drilling_step).at("CBOFnV")[0]);
}


QString Drilling::GetStatusString() const
{
  return GetStatusString(current_step_);
}

void Drilling::maintainRuntimeSettings(int drilling_step) {
  Runner::RuntimeSettings* rts = runtime_settings_.value(drilling_step);
  runtime_settings_.value(drilling_step)->paths().SetPath(Paths::OUTPUT_DIR, original_output_dir_);
  runtime_settings_.insert(drilling_step+1, rts);
}

void Drilling::setOptRuntimeSettings(int drilling_step, int argc, const char** argv) {
  runtime_settings_.insert(drilling_step, new Runner::RuntimeSettings(argc, argv));
  original_output_dir_ = runtime_settings_.value(drilling_step)->paths().GetPath(Paths::OUTPUT_DIR);
}

void Drilling::setOptRuntimeSettings(int drilling_step, Runner::RuntimeSettings* rts) {
  if (settings_->drilling().drilling_schedule.execution_modes.value(drilling_step) == settings_->drilling().drilling_schedule.Serial) {
    rts->setRunnerType(rts->SERIAL);
  } else if ((settings_->drilling().drilling_schedule.execution_modes.value(drilling_step) == settings_->drilling().drilling_schedule.Parallel)) {
    rts->setRunnerType(rts->MPISYNC);
  }

  runtime_settings_.insert(drilling_step, rts);
  original_output_dir_ = runtime_settings_.value(drilling_step)->paths().GetPath(Paths::OUTPUT_DIR);
}

void Drilling::modelUpdate(int drilling_step) {
   /* TODO: implement interface to SLB code for model update
    * Remember to check if the model update should be performed in a the full-scale model or
    * in a surrogate model instead. Update the information in the drilling object afterwards.
    *
    * The updated model should be put in the folder:  QString output_dir = QString::fromStdString(rts->paths().GetPath(Paths::OUTPUT_DIR)) + QString("/drilling_step_%1").arg(drilling_step)
    *
    */

  std::string sim_driver_file = runtime_settings_.value(drilling_step)->paths().GetPath(Paths::SIM_DRIVER_FILE);
  std::string grid_file = runtime_settings_.value(drilling_step)->paths().GetPath(Paths::GRID_FILE);
  std::string driver_file = runtime_settings_.value(drilling_step)->paths().GetPath(Paths::DRIVER_FILE);

  std::string realization = "r" + std::to_string(current_model_);

  sim_driver_file.replace(sim_driver_file.find(realization), realization.size(), "r" + std::to_string(current_model_+1));
  grid_file.replace(grid_file.find(realization), realization.size(), "r" + std::to_string(current_model_+1));
  driver_file.replace(driver_file.find(realization), realization.size(), "r" + std::to_string(current_model_+1));

  Runner::RuntimeSettings* rts = runtime_settings_.value(drilling_step);
  runtime_settings_.insert(drilling_step+1, rts);
  runtime_settings_.value(drilling_step+1)->paths().SetPath(Paths::SIM_DRIVER_FILE, sim_driver_file);
  runtime_settings_.value(drilling_step+1)->paths().SetPath(Paths::GRID_FILE, grid_file);
  runtime_settings_.value(drilling_step+1)->paths().SetPath(Paths::DRIVER_FILE, driver_file);

  current_model_++;
}

void Drilling::runOptimization(int drilling_step) {
  current_step_ = drilling_step;

  //TODO: allow change of deck and EGRID files in each drilling step (folder structure with models)
  QString output_dir = QString::fromStdString(original_output_dir_) + QString("/drilling_step_%1").arg(drilling_step);
  Utilities::FileHandling::CreateDirectory(output_dir);

  runtime_settings_.value(drilling_step)->paths().SetPath(Paths::OUTPUT_DIR,output_dir.toStdString());

  // warm-starting optimization
  Runner::MainRunner* runner = new Runner::MainRunner(runtime_settings_.value(drilling_step), best_case_, mso_);

  //   Configure the optimization with the drilling workflow settings.
  for (auto well: runner->getSettings()->model()->wells()) {
    for (int i=0; i < well.completions.size(); i++) {
      // current drilling step
      if (drilling_schedule_->isVariableDrillingPoints().value(i))
        well.completions.at(i).variable_placement = true;

      if (drilling_schedule_->isVariableCompletions().value(i))
        well.completions.at(i).variable_strength = true;
    }
  }
  runner->Execute();

  Optimization::Optimizer* opt = runner->getOptimizer();

  setWellOptimalVariables(opt->GetOptimalBinaryVariables(), opt->GetOptimalIntegerVariables(), opt->GetOptimalVariables(), drilling_step);
  setWellOptimizationValues(opt->GetOptimalValues(), drilling_step);

  best_case_ = opt->GetTentativeBestCase();
  mso_ = new ModelSynchronizationObject(runner->getModel());

}
}
}