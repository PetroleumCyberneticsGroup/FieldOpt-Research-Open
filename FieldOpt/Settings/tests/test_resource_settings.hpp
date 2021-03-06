/***********************************************************
Copyright (C) 2015-2018
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

#ifndef FIELDOPT_TEST_RESOURCE_SETTINGS_H
#define FIELDOPT_TEST_RESOURCE_SETTINGS_H

#include "Settings/settings.h"
#include "Settings/tests/test_resource_example_file_paths.hpp"
#include "Runner/tests/test_resource_runner.hpp"

namespace TestResources {

class TestResourceSettings {
 protected:
  TestResourceSettings() {
    Settings::VerbParams vp_;

    // -----------------------------------------------------
    paths_.SetPath(Paths::DRIVER_FILE,
                   ExampleFilePaths::driver_example_);

    paths_.SetPath(Paths::OUTPUT_DIR,
                   ExampleFilePaths::directory_output_);

    paths_.SetPath(Paths::SIM_DRIVER_FILE,
                   ExampleFilePaths::deck_horzwel_);

    paths_.SetPath(Paths::GRID_FILE,
                   ExampleFilePaths::grid_5spot_);

    if(vp_.vSET >= 5) {
      string im = "test: settings_full_(";
      im += paths_.GetPath(Paths::DRIVER_FILE) + ")";
      ext_warn(im, "", "TestResourceSettings", vp_.lnw, 1);
    }

    settings_full_ = new Settings::Settings(paths_);
    settings_optimizer_ = settings_full_->optimizer();
    settings_simulator_ = settings_full_->simulator();
    settings_model_ = settings_full_->model();
    settings_well_ = settings_full_->model()->wells().at(0);

    // -----------------------------------------------------
    paths_hybridopt_.SetPath(Paths::DRIVER_FILE,
                             ExampleFilePaths::hybridopt_driver_example_);

    paths_hybridopt_.SetPath(Paths::OUTPUT_DIR,
                             ExampleFilePaths::directory_output_);

    paths_hybridopt_.SetPath(Paths::SIM_DRIVER_FILE,
                             ExampleFilePaths::deck_horzwel_);

    paths_hybridopt_.SetPath(Paths::GRID_FILE,
                             ExampleFilePaths::grid_5spot_);

    if(vp_.vSET >= 5) {
      string im = "test: settings_hybridopt_full_(";
      im += paths_hybridopt_.GetPath(Paths::DRIVER_FILE) + ")";
      ext_warn(im, "", "TestResourceSettings", vp_.lnw, 1);
    }

    settings_hybridopt_full_ = new Settings::Settings(paths_hybridopt_);
    settings_hybridopt_optimizer_ = settings_hybridopt_full_->optimizer();
    settings_hybridopt_simulator_ = settings_hybridopt_full_->simulator();
    settings_hybridopt_model_ = settings_hybridopt_full_->model();

    // -----------------------------------------------------
    runner_resources_ = new TestResources::RunnerResources();

    auto en_paths = runner_resources_->rts_en_5spot_->paths();
    if(vp_.vSET >= 5) {
      string im = "test: settings_en_5spot_full_(";
      im += en_paths.GetPath(Paths::DRIVER_FILE) + ")";
      ext_warn(im, "", "TestResourceSettings", vp_.lnw, 1);
    }
    settings_en_5spot_full_ = new Settings::Settings(en_paths);

    // -----------------------------------------------------
    paths_olympr37_ = runner_resources_->rts_olympr37_T01_->paths();
    if(vp_.vSET >= 5) {
      string im = "test: paths_olympr37_full_(";
      im += paths_olympr37_.GetPath(Paths::DRIVER_FILE) + ")";
      ext_warn(im, "", "TestResourceSettings", vp_.lnw, 1);
    }
    settings_olympr37_ = new Settings::Settings(paths_olympr37_);
    settings_olympr37_opt_ = settings_olympr37_->optimizer();
    settings_olympr37_sim_ = settings_olympr37_->simulator();
    settings_olympr37_mod_ = settings_olympr37_->model();
  }

  TestResources::RunnerResources *runner_resources_;
  Settings::Settings *settings_en_5spot_full_;

  Settings::Settings *settings_full_;
  Settings::Optimizer *settings_optimizer_;
  Settings::Simulator *settings_simulator_;
  Settings::Model *settings_model_;
  Settings::Model::Well settings_well_;
  Paths paths_;

  Settings::Settings *settings_hybridopt_full_;
  Settings::Optimizer *settings_hybridopt_optimizer_;
  Settings::Simulator *settings_hybridopt_simulator_;
  Settings::Model *settings_hybridopt_model_;
  Paths paths_hybridopt_;

  Settings::Settings *settings_olympr37_;
  Settings::Optimizer *settings_olympr37_opt_;
  Settings::Simulator *settings_olympr37_sim_;
  Settings::Model *settings_olympr37_mod_;
  Paths paths_olympr37_;

};

}

#endif //FIELDOPT_TEST_RESOURCE_SETTINGS_H
