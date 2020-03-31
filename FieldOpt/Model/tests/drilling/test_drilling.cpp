/******************************************************************************
   Copyright (C) 2020- Thiago Lima Silva <thiagolims@gmail.com>

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

#include <gtest/gtest.h>
#include <drilling/drilling.h>
#include "Model/tests/test_resource_model.h"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Optimization/tests/test_resource_test_functions.h"
#include "serial_runner.h"

namespace {

class DrillingTest :
    public ::testing::Test,
    public TestResources::TestResourceOptimizer,
    public TestResources::TestResourceGrids {
 protected:
  DrillingTest() {}

  void runOptimization(int drilling_step);

};

TEST_F(DrillingTest, Constructor) {
  EXPECT_TRUE(true);
}

TEST_F(DrillingTest, ChildObjects) {
  EXPECT_NO_THROW(settings_model_->drilling());
  EXPECT_NO_THROW(settings_model_->drilling().well_name);
  EXPECT_NO_THROW(settings_model_->drilling().drilling_schedule);
}

TEST_F(DrillingTest, DrillingClasses) {
  EXPECT_NO_THROW(new Model::Drilling::Drilling(settings_model_, nullptr));
  EXPECT_NO_THROW(new Model::Drilling::DrillingSchedule(settings_model_, nullptr));
}

TEST_F(DrillingTest, DrillingScheduleDataStructure) {
  settings_model_->readDrilling(json_settings_drilling_);
  Settings::Model::Drilling drilling_settings = settings_model_->drilling();

  Model::Drilling::Drilling *drilling = new Model::Drilling::Drilling(settings_model_, nullptr);
  Model::Drilling::DrillingSchedule *schedule = new Model::Drilling::DrillingSchedule(settings_model_, nullptr);

  bool is_drilling_object_correct = true;
  bool dbg = false;

  QList<int> steps = schedule->getSteps();
  if (steps.size() != drilling_settings.drilling_schedule.drilling_steps.size())
    is_drilling_object_correct = false;

  for (int i = 0; i < steps.size(); i++) {
    // check drilling steps and time steps
    if (steps.value(i) != drilling_settings.drilling_schedule.drilling_steps.value(i))
      is_drilling_object_correct = false;

    if (schedule->getTimeSteps().value(i) != drilling_settings.drilling_schedule.time_steps.value(i))
      is_drilling_object_correct = false;

    // check model types
    if (schedule->getModelTypes().value(i) != drilling_settings.drilling_schedule.model_types.value(i))
      is_drilling_object_correct = false;

    // check drilling completions
    if (schedule->isVariableCompletions().value(i) != drilling_settings.drilling_schedule.is_variable_completions.value(i))
      is_drilling_object_correct = false;

    // check model updates
    if (schedule->isModelUpdates().value(i) != drilling_settings.drilling_schedule.is_model_updates.value(i))
      is_drilling_object_correct = false;

    // check the drilling points
    QList<Model::Drilling::DrillingSchedule::DrillingPoint> points = schedule->getDrillingPoints().value(i);
    QList<Settings::Model::Drilling::DrillingPoint> points_settings = drilling_settings.drilling_schedule.drilling_points.value(i);
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
  EXPECT_NO_THROW(settings_model_->readDrilling(json_settings_drilling_));
  bool dbg = false;
  if (dbg) {
    Settings::Model::Drilling drilling = settings_model_->drilling();
    cout << "wellName:" << drilling.well_name.toStdString() << endl;

    Settings::Model::Drilling::DrillingSchedule schedule = drilling.drilling_schedule;
    QList<int> steps = schedule.drilling_steps;
    for (int i = 0; i < steps.size(); i++) {
      cout << "step: " << i << endl;
      cout << "TimeStep: " << schedule.time_steps.value(i) << endl;

      cout << "drilling points:" << endl;
      QList<Settings::Model::Drilling::DrillingPoint> points = schedule.drilling_points.value(i);
      cout << "points.size()=" << points.size() << endl;
      for (int j = 0; j < points.size(); j++) {
        cout << "point " << j << ": ";
        cout << "[" << points.value(j).x << "," << points.value(j).y << "," << points.value(j).z << "]" << endl;
      }

      cout << "Operation:";
      switch (schedule.drilling_operations.value(i)) {
        case drilling.drilling_schedule.StartDrilling:
          cout << "StartDrilling" << endl;
          break;
        case drilling.drilling_schedule.Drilling:
          cout << "Drilling" << endl;
          break;
        case drilling.drilling_schedule.PullingOutOfHole:
          cout << "PullingOutOfHole" << endl;
          break;

        default:cout << "undefined" << endl;
      }

      cout << "ModelType:";
      switch (schedule.model_types.value(i)) {
        case drilling.drilling_schedule.TrueModel:cout << "TrueModel" << endl;
          break;
        case drilling.drilling_schedule.Surrogate:cout << "Surrogate" << endl;
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
  settings_model_->readDrilling(json_settings_drilling_);
  Settings::Model::Drilling drilling_settings = settings_model_->drilling();

  Model::Drilling::Drilling *drilling = new Model::Drilling::Drilling(settings_model_, nullptr);
  Model::Drilling::DrillingSchedule *schedule = new Model::Drilling::DrillingSchedule(settings_model_, nullptr);
  bool is_drilling_workflow_completed = false;

  QList<int> drilling_steps = schedule->getSteps();
  QList<int> time_steps = schedule->getSteps();
  for (int i: drilling_steps) {
    cout << "drilling_step:" << i << endl;
    int ts = time_steps.value(i);

    // Model update
    if (schedule->isModelUpdates().value(i)) {
      drilling->modelUpdate();
    }

    // Optimization
    if ((schedule->isVariableDrillingPoints().value(i)) || (schedule->isVariableCompletions().value(i))) {
      runOptimization(i);
    }
  }

  is_drilling_workflow_completed = true;
  EXPECT_TRUE(is_drilling_workflow_completed);
}

void DrillingTest::runOptimization(int drilling_step) {
  // TODO: prepare the JSON/variables vector to launch the optimization
  // TODO: change the deck and EGRID files in each drilling step

  int argc = 16;
  const char *argv[16] = {"FieldOpt",
                          TestResources::ExampleFilePaths::driver_5pot_icds.c_str(),
                          TestResources::ExampleFilePaths::directory_output_.c_str(),
                          "-g", TestResources::ExampleFilePaths::grid_5spot_icds.c_str(),
                          "-s", TestResources::ExampleFilePaths::deck_5spot_icds.c_str(),
                          "-b", "./",
                          "-r", "serial",
                          "-f",
                          "-v", "2",
                          "-t", "1000"};

  Runner::RuntimeSettings *rts = new Runner::RuntimeSettings(argc, argv);

  QString output_dir = QString::fromStdString(rts->paths().GetPath(Paths::OUTPUT_DIR)) + QString("/drilling_step_%1").arg(drilling_step);
  Utilities::FileHandling::CreateDirectory(output_dir);
  rts->paths().SetPath(Paths::OUTPUT_DIR,output_dir.toStdString());

  Runner::SerialRunner serial_runner = Runner::SerialRunner(rts);
  serial_runner.Execute();



  //TODO: save optimal variables, optimal case, optimal results in the drilling object
  Optimization::Optimizer* opt = serial_runner.getOptimizer();
  //auto results = opt.GetValues(); // contains the results
}


}






