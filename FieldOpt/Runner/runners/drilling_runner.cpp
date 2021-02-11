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

#include <Utilities/time.hpp>
#include "drilling_runner.h"


namespace Runner {

DrillingRunner::DrillingRunner(Runner::RuntimeSettings *runtime_settings)
    : AbstractRunner(runtime_settings)
{
  InitializeLogger();
  InitializeSettings();
  InitializeModel();
  InitializeSimulator();

  EvaluateBaseModel();
  InitializeObjectiveFunction();

  InitializeBaseCase();
  InitializeOptimizer();
  InitializeBookkeeper();

  InitializeDrillingWorkflow();

  FinalizeInitialization(true);
}


void DrillingRunner::InitializeDrillingWorkflow() {
  Settings::Model::Drilling drilling_settings = settings_->model()->drilling();
  drilling_ = new Model::Drilling::Drilling(settings_->model(), nullptr);
}

void DrillingRunner::Execute() {
  QList<int> drilling_steps = drilling_->getDrillingSchedule()->getSteps();

  drilling_->setOptRuntimeSettings(0, runtime_settings_);

  drilling_->createLogFile(0);

  for (int i: drilling_steps) {
    double ts = drilling_->getDrillingSchedule()->getTimeSteps().value(i);

    // Optimization
    if ((drilling_->getDrillingSchedule()->isVariableDrillingPoints().value(i)) || (drilling_->getDrillingSchedule()->isVariableCompletions().value(i))) {
      drilling_->runOptimization(i);
    }

    // Model update
    if ((i < drilling_steps.size()-1) && (drilling_->getDrillingSchedule()->isModelUpdates().value(i))) {
      drilling_->modelUpdate(i);
    } else {
      drilling_->maintainRuntimeSettings(i);
    }

  }

}


}