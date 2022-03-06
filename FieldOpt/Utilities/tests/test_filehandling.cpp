/***********************************************************
Created: 28.10.2015 2015 by einar
Copyright (C) 2015
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

#include "filehandling.hpp"
#include "Settings/tests/test_resource_example_file_paths.hpp"
#include "Utilities/verbosity.h"

namespace {

using ::Utilities::FileHandling::FileExists;
using ::Utilities::FileHandling::DirExists;
using ::Utilities::FileHandling::ReadFileToStringList;
using TestResources::ExampleFilePaths::driver_example_;

class FileHandlingTest : public ::testing::Test {
 protected:
  FileHandlingTest() {}
  virtual ~FileHandlingTest() {}

  virtual void SetUp() {}
  virtual void TearDown() {}

  Settings::VerbParams vp_;
  string md_ = "Utilities::tests";
  string cl_ = "FileHandlingTest";
};

TEST_F(FileHandlingTest, Existance) {
  EXPECT_TRUE(FileExists(driver_example_, vp_));
  EXPECT_FALSE(FileExists(driver_example_ + "wrong", vp_));

  EXPECT_FALSE(DirExists(driver_example_, vp_, md_, cl_));
}

TEST_F(FileHandlingTest, FileReading) {
  EXPECT_LE(155, ReadFileToStringList(QString::fromStdString(driver_example_))->size());
}


}
