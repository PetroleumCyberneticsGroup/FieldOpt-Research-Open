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
#include "Simulation/simulator_interfaces/flowsimulator.h"

namespace {

class FlowDriverFileWriterTest
  : public ::testing::Test, TestResources::TestResourceModel {
 protected:
  FlowDriverFileWriterTest() {}
  virtual ~FlowDriverFileWriterTest() {}
  virtual void SetUp() {}
  Simulation::FlowSimulator *simulator_;
};

TEST_F(FlowDriverFileWriterTest, Initialization) {
}

TEST_F(FlowDriverFileWriterTest, WriteDriverFile) {
//        simulator_->Evaluate();
}

}



