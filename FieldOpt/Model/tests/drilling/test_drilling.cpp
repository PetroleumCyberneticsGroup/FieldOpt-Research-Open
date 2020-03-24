#include <gtest/gtest.h>
#include "Model/tests/test_resource_model.h"





namespace {

class DrillingTest : public ::testing::Test,  public TestResources::TestResourceSettings {
 protected:
  DrillingTest() { }

};

TEST_F(DrillingTest, Constructor) {
    EXPECT_TRUE(true);
}

TEST_F(DrillingTest, ChildObjects) {
  EXPECT_NO_THROW(settings_model_->drilling());
  EXPECT_NO_THROW(settings_model_->drilling().well_name);
  EXPECT_NO_THROW(settings_model_->drilling().drilling_schedule);
}


TEST_F(DrillingTest, ParseJson) {
  EXPECT_NO_THROW(settings_model_->readDrilling(json_settings_drilling_));
  bool dbg = true;
  if (dbg) {
    Settings::Model::Drilling drilling = settings_model_->drilling();
    cout << "wellName:" << drilling.well_name.toStdString() << endl;

    Settings::Model::Drilling::DrillingSchedule schedule = drilling.drilling_schedule;
    QList<int> steps = schedule.drilling_steps_;
    for (int i=0; i < steps.size(); i++) {
      cout << "step: " << i << endl;
      cout << "TimeStep: " << schedule.time_steps_.value(i) << endl;

      cout << "drilling points:" << endl;
      QList<Settings::Model::Drilling::DrillingPoint> points = schedule.drilling_points_.value(i);
      cout << "points.size()=" << points.size() << endl;
      for (int j=0; j < points.size(); j++) {
        cout << "point " << j << ": ";
        cout << "[" << points.value(j).x << "," << points.value(j).y << "," << points.value(j).z << "]" << endl;
      }

      cout << "Operation:";
      switch(schedule.drilling_operations.value(i)) {
        case drilling.drilling_schedule.StartDrilling:
          cout << "StartDrilling" << endl;
          break;
        case drilling.drilling_schedule.Drilling:
          cout << "Drilling" << endl;
          break;

        case drilling.drilling_schedule.PullingOutOfHole:
          cout << "PullingOutOfHole" << endl;
        default:
          cout << "undefined" << endl;
      }

      cout << "ModelType:";
      switch (schedule.model_types.value(i)) {
        case drilling.drilling_schedule.TrueModel:
          cout << "TrueModel" << endl;
          break;
        case drilling.drilling_schedule.Surrogate:
          cout << "Surrogate" << endl;
          break;
        default:
          cout << "undefined" << endl;
      }

      cout << "OptimizeDrillingPoints:" << schedule.is_variable_drilling_points.value(i) << endl;
      cout << "OptimizeCompletion:" << schedule.is_variable_completions.value(i) << endl;
      cout << "ModelUpdate:" << schedule.is_model_update.value(i) << endl;
    }
  }
}

}
