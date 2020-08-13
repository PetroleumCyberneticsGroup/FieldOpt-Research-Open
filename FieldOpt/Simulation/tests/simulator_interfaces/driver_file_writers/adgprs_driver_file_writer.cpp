/***********************************************************
Copyright (C) 2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020 Mathias Bellout
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
#include "Model/tests/test_resource_model.h"
#include "Simulation/simulator_interfaces/adgprssimulator.h"

namespace {

class AdgprsDriverFileWriterTest : public ::testing::Test, TestResources::TestResourceModel {
 protected:
  AdgprsDriverFileWriterTest() {
    settings_full_->paths().SetPath(Paths::BUILD_DIR, TestResources::ExampleFilePaths::bin_dir_);
    settings_full_->paths().SetPath(Paths::SIM_DRIVER_FILE, TestResources::ExampleFilePaths::gprs_drv_5spot_);
    settings_full_->paths().SetPath(Paths::SIM_DRIVER_DIR, GetParentDirPath(settings_full_->paths().GetPath(Paths::SIM_DRIVER_FILE)));
    simulator_ = new Simulation::AdgprsSimulator(settings_full_, model_);
  }
  virtual ~AdgprsDriverFileWriterTest() {}
  virtual void SetUp() {}
  Simulation::AdgprsSimulator *simulator_;
};

TEST_F(AdgprsDriverFileWriterTest, Initialization) {
}

}


