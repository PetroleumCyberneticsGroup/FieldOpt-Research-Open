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
#include <algorithm>

namespace Model {
namespace Drilling {

Drilling::Drilling(Settings::Model *settings, Properties::VariablePropertyContainer *drilling_variables) {
  drilling_variables_  = drilling_variables;
  settings_   = settings;


  well_name_ = settings->drilling().well_name;
  drilling_schedule_ = new DrillingSchedule(settings_, drilling_variables_);

  current_step_ = 0;
};

void Drilling::setWellOptimalVariables(const std::map<string, QHash<QUuid, double>>& opt_var, int drilling_step) {
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

void Drilling::modelUpdate(int drilling_step) {
   /* TODO: implement interface to SLB model update code.
    * Remember to check if the model update should be performed in a the full-scale model or
    * in a surrogate model instead. Update the information in the drilling object afterwards.
    *
    * The updated model should be put in the folder:  QString output_dir = QString::fromStdString(rts->paths().GetPath(Paths::OUTPUT_DIR)) + QString("/drilling_step_%1").arg(drilling_step)
    *
    */
  //
}

void Drilling::runOptimization(int drilling_step) {
  // TODO: prepare the JSON/variables vector to launch the optimization
  // TODO: change the deck and EGRID files in each drilling step

  current_step_ = drilling_step;

  int argc = 16;
  const char *argv[16] = {"FieldOpt",
                          TestResources::ExampleFilePaths::driver_5pot_icds.c_str(),
                          TestResources::ExampleFilePaths::directory_output_.c_str(),
                          "-g", TestResources::ExampleFilePaths::grid_5spot_icds.c_str(),
                          "-s", TestResources::ExampleFilePaths::deck_5spot_icds.c_str(),
                          "-b", "./",
                          "-r", "serial",
                          "-f",
                          "-v", "0",
                          "-t", "1000"};

  Runner::RuntimeSettings *rts = new Runner::RuntimeSettings(argc, argv);

  QString output_dir = QString::fromStdString(rts->paths().GetPath(Paths::OUTPUT_DIR)) + QString("/drilling_step_%1").arg(drilling_step);
  Utilities::FileHandling::CreateDirectory(output_dir);
  rts->paths().SetPath(Paths::OUTPUT_DIR,output_dir.toStdString());

  Runner::SerialRunner serial_runner = Runner::SerialRunner(rts);
  serial_runner.Execute();

  //TODO: save path to the case (?)
  Optimization::Optimizer* opt = serial_runner.getOptimizer();

  map<string, QHash<QUuid, double>> opt_variables = opt->GetOptimalVariables();
  map<string, vector<double>> opt_values = opt->GetOptimalValues();

  setWellOptimalVariables(opt_variables, drilling_step);
  setWellOptimizationValues(opt_values, drilling_step);
}

}
}