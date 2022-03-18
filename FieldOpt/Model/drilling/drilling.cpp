/***********************************************************
Copyright (C) 2020-
Thiago L. Silva <thiagolims@gmail.com>
Mathias Bellout <chakibbb-pcg@gmail.com>

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

#include "drilling.h"
#include "paths.h"
#include "Optimization/case.h"

namespace Model {
namespace Drilling {

Drilling::Drilling(Settings::Drilling *settings,
                   Properties::VarPropContainer *drilling_variables) {
  drilling_variables_ = drilling_variables;
  setd_ = settings;

  well_name_ = setd_->well_name;

  drilling_schedule_ = new DrillingSchedule(setd_, drilling_variables_);

  local_optimizer_settings_ = setd_->local_optimizer_settings;
  global_optimizer_settings_ = setd_->global_optimizer_settings;

  current_step_ = 0;
  current_model_ = 0;

  skip_optimization_ = false;

  best_case_ = nullptr;
  base_case_ = nullptr;
  mso_ = nullptr;
  best_objective_ = 0.0;
};

void Drilling::setWellOptimalVariables(const QList<QPair<QUuid, bool>>& opt_bin_var,
                                       const QList<QPair<QUuid, int>>& opt_int_var,
                                       const QList<QPair<QUuid, double>>& opt_var,
                                       int drilling_step) {
  optimal_int_variables_.insert(drilling_step, opt_int_var);
  optimal_bin_variables_.insert(drilling_step, opt_bin_var);
  optimal_variables_.insert(drilling_step, opt_var);
}

void Drilling::setWellOptimizationValues(
  const std::map<string, std::vector<double>>& opt_val,
  int drilling_step) {
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

void Drilling::setOptRuntimeSettings(int drilling_step,
                                     int argc, const char** argv) {
  runtime_settings_.insert(drilling_step, new Runner::RuntimeSettings(argc, argv));
  original_output_dir_ = runtime_settings_.value(drilling_step)->paths().GetPath(Paths::OUTPUT_DIR);
}

void Drilling::setOptRuntimeSettings(int drilling_step,
                                     Runner::RuntimeSettings* rts) {
  if (setd_->drilling_sched().execution_modes.value(drilling_step) == setd_->drilling_sched().Serial) {
    rts->setRunnerType(rts->SERIAL);
  } else if ((setd_->drilling_sched().execution_modes.value(drilling_step) == setd_->drilling_sched().Parallel)) {
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
  Utilities::FileHandling::CreateDir(output_dir, setd_->verbParams());

  runtime_settings_.value(drilling_step)->paths().SetPath(Paths::OUTPUT_DIR,output_dir.toStdString());

  Runner::MainRunner* runner;
  if (drilling_schedule_->isWarmStart().value(current_step_)) {
    // warm-starting optimization
    runner = new Runner::MainRunner(runtime_settings_.value(drilling_step), best_case_, mso_);
  } else {
    runner = new Runner::MainRunner(runtime_settings_.value(drilling_step));
  }

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

  if (drilling_schedule_->getOptimizerSettings().value(drilling_step) != nullptr)
    runner->ReplaceOptimizer(drilling_schedule_->getOptimizerSettings().value(drilling_step));

  base_case_ = runner->getOptimizer()->GetTentBestCase();

  //!<Trigger #1: model deviation>
  double model_dev = abs((base_case_->objf_value() - best_objective_) / best_objective_);
  if (drilling_schedule_->getOptimizationTriggers().contains(drilling_step)) {
    double min_model_dev = drilling_schedule_->getOptimizationTriggers().value(drilling_step).min_model_deviation;
    double max_model_dev = drilling_schedule_->getOptimizationTriggers().value(drilling_step).max_model_deviation;
    if (min_model_dev >=0) {
      runner->getOptimizer()->setMinModelDeviation(min_model_dev);

      if (model_dev <= min_model_dev) {
        skip_optimization_ = true;
      } else if ((model_dev > min_model_dev)  && (model_dev <= max_model_dev ))   {
        runner->ReplaceOptimizer(local_optimizer_settings_);
      } else if ((model_dev > max_model_dev) && (max_model_dev >=0)) {
        runner->ReplaceOptimizer(global_optimizer_settings_);
      }
    }
  }

  //!<Trigger #2: sufficient improvement>
  if (drilling_schedule_->getOptimizationTriggers().contains(drilling_step)) {
    double max_obj_improvement = drilling_schedule_->getOptimizationTriggers().value(drilling_step).max_objective_improvement;
    if (max_obj_improvement >= 0)
      runner->getOptimizer()->setSufficientImprovementTolerance(max_obj_improvement);
  }

  if (!skip_optimization_) {
    runner->Execute();
  }

  Optimization::Optimizer *opt = runner->getOptimizer();

  setWellOptimalVariables(opt->GetOptimalBinaryVariables(),
                          opt->GetOptimalIntegerVariables(),
                          opt->GetOptimalVariables(),
                          drilling_step);
  setWellOptimizationValues(opt->GetOptimalValues(), drilling_step);

  best_case_ = opt->GetTentBestCase();
  best_objective_ = best_case_->objf_value();

  double obj_improvement = abs((best_objective_-base_case_->objf_value())/base_case_->objf_value());

  printIteration(drilling_step, model_dev, obj_improvement, skip_optimization_);


  mso_ = new ModelSynchronizationObject(runner->getModel());

  skip_optimization_ = false;
}

void Drilling::createLogFile(int drilling_step) {
  QString output_dir = QString::fromUtf8(runtime_settings_.value(drilling_step)->paths().GetPath(Paths::OUTPUT_DIR).c_str());

  if (output_dir.length() > 0) {
    Utilities::FileHandling::CreateDir(output_dir, setd_->verbParams());
  }

  dw_log_path_ = output_dir + "/log_drilling_workflow.csv";

  // Delete existing logs if --force flag is on
  if (Utilities::FileHandling::FileExists(dw_log_path_, setd_->verbParams())) {
    Utilities::FileHandling::DeleteFile(dw_log_path_, setd_->verbParams());
  }

  const QString tr_log_header = "ds dev obj skip";
  Utilities::FileHandling::WriteLineToFile(tr_log_header, dw_log_path_);

}

void Drilling::printIteration(int drilling_step,
                              double model_deviation,
                              double obj_improvement,
                              bool skip_opt) {
    stringstream ss;

    ss << setw(4) << right << drilling_step << setprecision(4)
       << setw(6) << right << model_deviation << setprecision(8)
       << setw(6) << right << obj_improvement << setprecision(8)
       << setw(6)  << right << skip_opt;

    string str = ss.str();
    Utilities::FileHandling::WriteLineToFile(QString::fromStdString(str), dw_log_path_);
}

}
}