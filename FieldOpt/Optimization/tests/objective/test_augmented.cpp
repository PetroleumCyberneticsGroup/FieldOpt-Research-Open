/***********************************************************
Created by bellout on 8/28/20.

Copyright (C) 2017-2020 Mathias Bellout
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

#include "Simulation/results/eclresults.h"
#include "Optimization/objective/augmented.h"

#include "Settings/tests/test_resource_example_file_paths.hpp"
#include "Settings/tests/test_resource_settings.hpp"
#include "Simulation/tests/test_resource_results.h"

using namespace Optimization::Objective;
using namespace Simulation::Results;
using namespace TestResources::ExampleFilePaths;
using namespace Printer;

namespace {

using ECLProp = Simulation::Results::ECLResults::Property;

class AugmentedObjTest : public TestResources::TestResourceResults,
                         public TestResources::RunnerResources {

 protected:
  AugmentedObjTest() {
    obj_olympr37_ = new Augmented(settings_olympr37_opt_,
                                  results_olympr37_,
                                  model_olympr37_);

    model_olympr37_ = new Model::Model(*settings_olympr37_,
                                       logger_olympr37_T01_);

    obj_olympr37_ext_ = new Augmented(settings_olympr37_opt_,
                                  results_olympr37_ext_,
                                  model_olympr37_ext_);

    model_olympr37_ext_ = new Model::Model(*settings_olympr37_,
                                       logger_olympr37_T02_);

  }

  ~AugmentedObjTest() override = default;

  Augmented *obj_olympr37_;
  Model::Model *model_olympr37_;

  Augmented *obj_olympr37_ext_;
  Model::Model *model_olympr37_ext_;

};

TEST_F(AugmentedObjTest, Value) {
  obj_olympr37_->setDbgFileName("dbg_aug_obj.py");
  obj_olympr37_->value();
}

TEST_F(AugmentedObjTest, ValueExt) {
  obj_olympr37_ext_->setDbgFileName("dbg_aug_obj_ext.py");
  obj_olympr37_ext_->value();
}

}