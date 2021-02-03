/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb.pcg@gmail.com>

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
#include "Model/tests/test_resource_model.h"
#include "Simulation/simulator_interfaces/driver_file_writers/ix_driver_file_writer.h"
#include <Settings/tests/test_resource_example_file_paths.hpp>
#include <simulator_interfaces/simulator.h>
#include "Utilities/stringhelpers.hpp"

using namespace Simulation;
using namespace TestResources;
namespace {

class IXSimulatorTest : public testing::Test, TestResourceModel {
 protected:
  IXSimulatorTest() {
    Paths paths;
    paths.SetPath(Paths::DRIVER_FILE, ExampleFilePaths::norne_driver_example_);
    paths.SetPath(Paths::OUTPUT_DIR, ExampleFilePaths::norne_test_output_);
    paths.SetPath(Paths::GRID_FILE, ExampleFilePaths::norne_grid_);
    paths.SetPath(Paths::BUILD_DIR, ExampleFilePaths::bin_dir_);
    paths.SetPath(Paths::SIM_DRIVER_FILE, ExampleFilePaths::norne_deck_);

    settings_norne_full_ = new Settings::Settings(paths);
    settings_norne_optimizer_ = settings_norne_full_->optimizer();
    settings_norne_simulator_ = settings_norne_full_->simulator();
    settings_norne_model_ = settings_norne_full_->model();
    norne_model_ = new Model::Model(*settings_norne_full_, logger_norne_);
//        simulator_ = new IXSimulator(settings_norne_full_, norne_model_);

  }
  Settings::Settings *settings_norne_full_;
  Settings::Optimizer *settings_norne_optimizer_;
  Settings::Simulator *settings_norne_simulator_;
  Settings::Model *settings_norne_model_;
  Model::Model *norne_model_;
};

TEST_F(IXSimulatorTest, DriverFileWriter) {
  Simulation::IXDriverFileWriter dfw = Simulation::IXDriverFileWriter(settings_norne_full_, norne_model_);
  dfw.WriteDriverFile("");
}

TEST_F(IXSimulatorTest, Evaluate) {
//    simulator_->Evaluate();
}
//
//TEST_F(IXSimulatorTest, CleanUp) {
//    //simulator_->CleanUp();
//}

}

