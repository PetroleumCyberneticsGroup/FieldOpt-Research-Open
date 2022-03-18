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

#include <gtest/gtest.h>
#include <drilling/drilling.h>
#include "Model/tests/test_resource_model.h"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Optimization/tests/test_resource_test_functions.h"

#include "drilling_runner.h"

namespace {

class DrillingTest :
  public ::testing::Test,
  public TestResources::TestResourceOptimizer,
  public TestResources::TestResourceGrids {
 protected:
  DrillingTest() = default;

  void runOptimization(int drilling_step);

};

TEST_F(DrillingTest, Constructor) {
  EXPECT_TRUE(true);
}

TEST_F(DrillingTest, ChildObjects) {
  EXPECT_NO_THROW(settings_drilling_);
  EXPECT_NO_THROW(settings_drilling_->well_name);
  EXPECT_NO_THROW(settings_drilling_->drilling_sched());
}

TEST_F(DrillingTest, DrillingClasses) {
  EXPECT_NO_THROW(new Model::Drilling::Drilling(settings_drilling_, nullptr));
  EXPECT_NO_THROW(new Model::Drilling::DrillingSchedule(settings_drilling_, nullptr));
}

TEST_F(DrillingTest, DrillingScheduleDataStructure) {
  settings_drilling_->readDrilling(json_settings_drilling_);
  Settings::Drilling* drilling_settings = settings_drilling_;

  auto *drilling = new Model::Drilling::Drilling(settings_drilling_, nullptr);
  auto *schedule = new Model::Drilling::DrillingSchedule(settings_drilling_, nullptr);

  bool is_drilling_object_correct = true;
  bool dbg = false;

  QList<int> steps = schedule->getSteps();
  if (steps.size() != drilling_settings->drilling_sched().drilling_steps.size())
    is_drilling_object_correct = false;

  for (int i = 0; i < steps.size(); i++) {
    // check drilling steps and time steps
    if (steps.value(i) != drilling_settings->drilling_sched().drilling_steps.value(i))
      is_drilling_object_correct = false;

    if (schedule->getTimeSteps().value(i) != drilling_settings->drilling_sched().time_steps.value(i))
      is_drilling_object_correct = false;

    // check model types
    if ((int)schedule->getModelTypes().value(i) != (int)drilling_settings->drilling_sched().model_types.value(i))
      is_drilling_object_correct = false;

    // check drilling completions
    if (schedule->isVariableCompletions().value(i) != drilling_settings->drilling_sched().is_variable_completions.value(i))
      is_drilling_object_correct = false;

    // check model updates
    if (schedule->isModelUpdates().value(i) != drilling_settings->drilling_sched().is_model_updates.value(i))
      is_drilling_object_correct = false;

    // check the drilling points
    QList<Model::Drilling::DrillingSchedule::DrillingPoint> points = schedule->getDrillingPoints().value(i);
    QList<Settings::Drilling::DrillingPoint> points_settings = drilling_settings->drilling_sched().drilling_points.value(i);
    if (points.size() == points_settings.size()) {
      for (int j = 0; j < points.size(); j++) {

        if ((points.value(j).x != points_settings.value(j).x) || (points.value(j).y != points_settings.value(j).y) || (points.value(j).z != points_settings.value(j).z)) {
          is_drilling_object_correct = false;
          break;
        } else {
          if (dbg) {
            cout << "point " << j << ": ";
            cout << "[" << points.value(j).x << "," << points.value(j).y << "," << points.value(j).z << "]" << endl;

            cout << "point_settings " << j << ": ";
            cout << "[" << points_settings.value(j).x << "," << points_settings.value(j).y << ","
                 << points_settings.value(j).z << "]" << endl;
          }
        }
      }
    } else {
      is_drilling_object_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_drilling_object_correct);
}

TEST_F(DrillingTest, ParseJson) {
  EXPECT_NO_THROW(settings_drilling_->readDrilling(json_settings_drilling_));
  bool dbg = false;
  if (dbg) {
    Settings::Drilling* drilling = settings_drilling_;
    cout << "wellName:" << drilling->well_name.toStdString() << endl;

    Settings::Drilling::DrillingSchedule schedule = drilling->drilling_sched();
    QList<int> steps = schedule.drilling_steps;
    for (int i = 0; i < steps.size(); i++) {
      cout << "step: " << i << endl;
      cout << "TimeStep: " << schedule.time_steps.value(i) << endl;

      cout << "drilling points:" << endl;
      QList<Settings::Drilling::DrillingPoint> points = schedule.drilling_points.value(i);
      cout << "points.size()=" << points.size() << endl;
      for (int j = 0; j < points.size(); j++) {
        cout << "point " << j << ": ";
        cout << "[" << points.value(j).x << "," << points.value(j).y << "," << points.value(j).z << "]" << endl;
      }

      cout << "Operation:";
      switch (schedule.drilling_operations.value(i)) {
        case Settings::Drilling::DrillingSchedule::DrillingOperation::StartDrilling:
          cout << "StartDrilling" << endl;
          break;
        case Settings::Drilling::DrillingSchedule::DrillingOperation::Drilling:
          cout << "Drilling" << endl;
          break;
        case Settings::Drilling::DrillingSchedule::DrillingOperation::PullingOutOfHole:
          cout << "PullingOutOfHole" << endl;
          break;

        default:cout << "undefined" << endl;
      }

      cout << "ModelType:";
      switch (schedule.model_types.value(i)) {
        case Settings::Drilling::DrillingSchedule::ModelType::TrueModel:cout << "TrueModel" << endl;
          break;
        case Settings::Drilling::DrillingSchedule::ModelType::Surrogate:cout << "Surrogate" << endl;
          break;
        default:cout << "undefined" << endl;
      }

      cout << "OptimizeDrillingPoints:" << schedule.is_variable_drilling_points.value(i) << endl;
      cout << "OptimizeCompletion:" << schedule.is_variable_completions.value(i) << endl;
      cout << "ModelUpdate:" << schedule.is_model_updates.value(i) << endl;
    }
  }
}

TEST_F(DrillingTest, DrillingRunner) {
  settings_drilling_->readDrilling(json_settings_drilling_);
  Settings::Drilling* drilling_settings = settings_drilling_;
  bool dbg = false;

  Model::Drilling::Drilling *drilling = new Model::Drilling::Drilling(settings_drilling_, nullptr);
  Model::Drilling::DrillingSchedule *schedule = drilling->getDrillingSchedule();

  int argc = 16;
  const char *argv[16] = {"FieldOpt",
                          TestResources::ExampleFilePaths::driver_5pot_icds_.c_str(),
                          TestResources::ExampleFilePaths::directory_output_icds_.c_str(),
                          "-g", TestResources::ExampleFilePaths::grid_5spot_icds_.c_str(),
                          "-s", TestResources::ExampleFilePaths::deck_5spot_icds_.c_str(),
                          "-b", "./",
                          "-r", "serial",
                          "-f",
                          "-v", "0",
                          "-t", "1000"};

  bool is_drilling_workflow_completed = false;

  QList<int> drilling_steps = schedule->getSteps();

  if (dbg) {
    cout << drilling->GetStatusStringHeader().toStdString() << endl;
  }
  drilling->setOptRuntimeSettings(0, argc, argv);
  for (int i: drilling_steps) {
    double ts = schedule->getTimeSteps().value(i);

    if (dbg) {
      cout << "drilling_step:" << i << endl;
    }

    // Optimization
    if ((schedule->isVariableDrillingPoints().value(i)) || (schedule->isVariableCompletions().value(i))) {
      drilling->runOptimization(i);
    }

    // Model update
    if ((i < drilling_steps.size()-1) && (schedule->isModelUpdates().value(i))) {
      drilling->modelUpdate(i);
    } else {
      drilling->maintainRuntimeSettings(i);
    }

    if (dbg) {
      cout << drilling->GetStatusString().toStdString() << endl;
    }
  }

  is_drilling_workflow_completed = true;
  EXPECT_TRUE(is_drilling_workflow_completed);
}

TEST_F(DrillingTest, DrillingRunnerWorkflow) {
  int argc = 16;
  const char *argv[16] = {"FieldOpt",
                          TestResources::ExampleFilePaths::driver_5pot_womud_drilling_triggers.c_str(),
                          TestResources::ExampleFilePaths::directory_output_womud_.c_str(),
                          "-g", TestResources::ExampleFilePaths::grid_5spot_icds_.c_str(),
                          "-s", TestResources::ExampleFilePaths::deck_5spot_icds_.c_str(),
                          "-b", "./",
                          "-r", "drillingWorkflow",
                          "-f",
                          "-v", "0",
                          "-t", "1000"};


  auto* rts = new Runner::RuntimeSettings(argc, argv);
  auto* drilling_runner = new Runner::DrillingRunner(rts);

  drilling_runner->Execute();
}


}






