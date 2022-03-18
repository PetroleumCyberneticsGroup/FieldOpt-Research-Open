/***********************************************************
Copyright (C) 2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<mathias.bellout@gmail.com>

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
#include "Runner/runtime_settings.h"
#include "Settings/tests/test_resource_example_file_paths.hpp"

namespace {

class RuntimeSettingsTest : public testing::Test {
  int argc = 16;
  const char *argv[16] = {"FieldOpt",
                          TestResources::ExampleFilePaths::driver_5spot_.c_str(),
                          TestResources::ExampleFilePaths::directory_output_.c_str(),
                          "-g", TestResources::ExampleFilePaths::grid_flow_5spot_.c_str(),
                          "-s", TestResources::ExampleFilePaths::deck_flow_5spot_.c_str(),
                          "-b", ".",
                          "-r", "mpisync",
                          "-f",
                          "-v", "2",
                          "-t", "1000"
  };

 protected:

  RuntimeSettingsTest() {
    rts = new Runner::RuntimeSettings(argc, argv);
  }
  ~RuntimeSettingsTest(){}
  Runner::RuntimeSettings *rts;
};

//    TEST_F(RuntimeSettingsTest, Verbosity) {
//        EXPECT_EQ(rts->verbosity_level(), 2);
//    }

}
