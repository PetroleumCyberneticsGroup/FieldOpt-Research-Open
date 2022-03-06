/***********************************************************
Copyright (C) 2016
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

#include <Model/tests/test_resource_model.h>
#include <gtest/gtest.h>
#include "Simulation/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/welspecs.h"

using namespace ::Simulation::ECLDriverParts;

namespace {

class DriverPartWelspecsTest : public ::testing::Test,
                               public TestResources::TestResourceModel {
 protected:
  DriverPartWelspecsTest(){
    welspecs_ = new Welspecs(model_->wells());
  }
  virtual ~DriverPartWelspecsTest(){}

  Welspecs *welspecs_;
};

TEST_F(DriverPartWelspecsTest, Constructor) {
  //std::cout << welspecs_->GetPartString().toStdString() << std::endl;
}

}
