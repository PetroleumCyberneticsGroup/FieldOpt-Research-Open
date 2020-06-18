/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb-pcg@gmail.com>

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
#include <QString>

#include "Settings/tests/test_resource_settings.hpp"
#include "Settings/tests/test_resource_multispline_wells_settings.hpp"

namespace {

//==========================================================
class ModelSettingsTest : public ::testing::Test,
                          public TestResources::TestResourceSettings {
 protected:
  ModelSettingsTest() {}

  virtual ~ModelSettingsTest() {}
  
  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(ModelSettingsTest, ControlTimes) {
  EXPECT_EQ(4, settings_model_->control_times().size());
  EXPECT_EQ(365, settings_model_->control_times().last());
}

TEST_F(ModelSettingsTest, ProducerWell) {
  Settings::Model::Well producer = settings_model_->wells().first();
  EXPECT_STREQ("PROD", producer.name.toLatin1().constData());
  EXPECT_EQ(Settings::Model::WellType::Producer, producer.type);
  EXPECT_EQ(Settings::Model::WellDefinitionType::WellBlocks, producer.definition_type);
  EXPECT_EQ(Settings::Model::PreferredPhase::Oil, producer.preferred_phase);
  EXPECT_FLOAT_EQ(0.75, producer.wellbore_radius);
  EXPECT_EQ(Settings::Model::Direction::X, producer.direction);
}

TEST_F(ModelSettingsTest, ProducerCompletions) {
  Settings::Model::Well producer = settings_model_->wells().first();
  EXPECT_EQ(4, producer.well_blocks.size());
  
  for (int i = 0; i < producer.well_blocks.size(); ++i) {
    EXPECT_EQ(i + 1, producer.well_blocks[i].i);
    EXPECT_EQ(4, producer.well_blocks[i].j);
    EXPECT_EQ(1, producer.well_blocks[i].k);
    QString expected_name = QString("WellBlock#PROD#%1").arg(i);
    EXPECT_STREQ(expected_name.toLatin1().constData(),
        producer.well_blocks[i].name.toLatin1().constData());
    
    if (i == 1) {
      EXPECT_FALSE(producer.well_blocks[i].is_variable);
      EXPECT_TRUE(producer.well_blocks[i].has_completion);
    } else {
      EXPECT_TRUE(producer.well_blocks[i].is_variable);
      EXPECT_TRUE(producer.well_blocks[i].has_completion);
      EXPECT_FLOAT_EQ(1.0, producer.well_blocks[i].completion.transmissibility_factor);
      
      EXPECT_EQ(Settings::Model::WellCompletionType::Perforation,
          producer.well_blocks[i].completion.type);
      QString expected_name = QString("Transmissibility#PROD#%1").arg(i);
      EXPECT_STREQ(expected_name.toLatin1().constData(),
                   producer.well_blocks[i].completion.name.toLatin1().constData());
    }
  }
}

TEST_F(ModelSettingsTest, ProducerControls) {
  Settings::Model::Well producer = settings_model_->wells().first();
  EXPECT_EQ(3, producer.controls.length());
  for (auto control : producer.controls) {
    QString expected_name = QString("BHP#PROD#%1").arg(control.time_step);
    
    EXPECT_STREQ(expected_name.toLatin1().constData(),
        control.name.toLatin1().constData());
    EXPECT_EQ(Settings::Model::ControlMode::BHPControl, control.control_mode);
    EXPECT_FLOAT_EQ(100.0, control.bhp);
    EXPECT_LE(0, control.time_step);
  }
}

TEST_F(ModelSettingsTest, InjectorWell) {
  Settings::Model::Well injector = settings_model_->wells()[1];
  EXPECT_STREQ("INJ", injector.name.toLatin1().constData());
  EXPECT_EQ(Settings::Model::WellType::Injector, injector.type);
  EXPECT_EQ(Settings::Model::PreferredPhase::Water, injector.preferred_phase);
  EXPECT_FLOAT_EQ(0.75, injector.wellbore_radius);
  EXPECT_EQ(Settings::Model::WellDefinitionType::WellSpline,
      injector.definition_type);
}


TEST_F(ModelSettingsTest, InjectorSpline) {
  Settings::Model::Well injector = settings_model_->wells()[1];
  
  EXPECT_FLOAT_EQ(12.0, injector.spline_heel.x);
  EXPECT_FLOAT_EQ(12.0, injector.spline_heel.y);
  EXPECT_FLOAT_EQ(1712.0, injector.spline_heel.z);
  EXPECT_TRUE(injector.spline_heel.is_variable);
  EXPECT_STREQ("SplinePoint#INJ#heel",
      injector.spline_heel.name.toLatin1().constData());
  
  EXPECT_FLOAT_EQ(36.0, injector.spline_toe.x);
  EXPECT_FLOAT_EQ(12.0, injector.spline_toe.y);
  EXPECT_FLOAT_EQ(1712.0, injector.spline_toe.z);
  EXPECT_TRUE(injector.spline_toe.is_variable);
  EXPECT_STREQ("SplinePoint#INJ#toe",
      injector.spline_toe.name.toLatin1().constData());
}

TEST_F(ModelSettingsTest, InjectorCompletions) {
  Settings::Model::Well injector = settings_model_->wells()[1];
  // The injector should not have any completions or well blocks
  EXPECT_EQ(0, injector.well_blocks.length());
}

TEST_F(ModelSettingsTest, InjectorControls) {
  Settings::Model::Well injector = settings_model_->wells()[1];
  EXPECT_EQ(1, injector.controls.length());
  EXPECT_EQ(0, injector.controls[0].time_step);
  EXPECT_EQ(Settings::Model::InjectionType::WaterInjection, injector.controls[0].injection_type);
  EXPECT_EQ(Settings::Model::WellState::WellOpen, injector.controls[0].state);
  EXPECT_EQ(Settings::Model::ControlMode::LRATControl, injector.controls[0].control_mode);
  EXPECT_FLOAT_EQ(1200, injector.controls[0].liq_rate);
  EXPECT_FALSE(injector.controls[0].is_variable);
  EXPECT_STREQ("Rate#INJ#0", injector.controls[0].name.toLatin1().constData());
}

TEST_F(ModelSettingsTest, MultisplineWell) {
  
  Paths paths_;
  QJsonObject partial_deck = model_adtl_pts;
  paths_.SetPath(Paths::SIM_DRIVER_FILE, TestResources::ExampleFilePaths::norne_deck_ );
  paths_.SetPath(Paths::SIM_SCH_FILE, TestResources::ExampleFilePaths::norne_sch_);
  // QJsonObject sim_json_;
  QJsonObject mod_json_;
  // sim_json_ = partial_deck["Simulator"].toObject();
  mod_json_ = partial_deck["Model"].toObject();
  // auto sim_settings = dc_6Settings::Simulator(sim_json_, paths_);
  auto mod_settings = Settings::Model(mod_json_, paths_);
}

}
