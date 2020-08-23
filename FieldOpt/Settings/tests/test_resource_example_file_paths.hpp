/***********************************************************
Created by einar on 6.2.16
Copyright (C) 2016
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020- Mathias Bellout
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

#ifndef FIELDOPT_TEST_RESOURCE_EXAMPLE_FILE_PATHS_H
#define FIELDOPT_TEST_RESOURCE_EXAMPLE_FILE_PATHS_H

#include "Utilities/system.hpp"
#include <string>

namespace TestResources {
namespace ExampleFilePaths {

using std::string;

inline std::string base_path() {
  char const* fieldopt_path = getenv("FIELDOPT_BUILD_ROOT");
  if (!is_env_var_set("FIELDOPT_BUILD_ROOT")) {
    return "..";
  }
  else {
    return get_env_var_value("FIELDOPT_BUILD_ROOT");
  }
}

static std::string bin_dir_                  = base_path() + "/bin/";
static std::string driver_example_           = base_path() + "/examples/driver.json";
static std::string hybridopt_driver_example_ = base_path() + "/examples/driver_hybridopt.json";
static std::string directory_output_         = base_path() + "/fieldopt-output";
static std::string norne_test_output_        = base_path() + "/fieldopt-output/norne_test";
static std::string norne_driver_example_     = base_path() + "/examples/ECLIPSE/norne-simplified/fo_driver.1_rate.cs.json";
static std::string norne_grid_               = base_path() + "/examples/ECLIPSE/norne-simplified/NORNE_SIMPLIFIED.EGRID";
static std::string norne_atw_grid_           = base_path() + "/examples/Flow/norne/NORNE_ATW2013.EGRID";
static std::string norne_deck_               = base_path() + "/examples/ECLIPSE/norne-simplified/NORNE_SIMPLIFIED.DATA";
static std::string norne_sch_                = base_path() + "/examples/ECLIPSE/norne-simplified/INCLUDE/BC0407_HIST01122006.SCH";
static std::string deck_horzwel_             = base_path() + "/examples/ECLIPSE/HORZWELL/HORZWELL.DATA";
static std::string grid_horzwel_             = base_path() + "/examples/ECLIPSE/HORZWELL/HORZWELL.EGRID";
static std::string ecl_base_horzwell         = base_path() + "/examples/ECLIPSE/HORZWELL/HORZWELL";

static std::string grid_5spot_               = base_path() + "/examples/ADGPRS/5spot/ECL_5SPOT.EGRID";
static std::string driver_5spot_             = base_path() + "/examples/ADGPRS/5spot/fo_driver_5vert_wells.json";

static std::string grid_flow_5spot_          = base_path() + "/examples/Flow/5spot/5SPOT.EGRID";
static std::string deck_flow_5spot_          = base_path() + "/examples/Flow/5spot/5SPOT.DATA";

static std::string gprs_drv_5spot_           = base_path() + "/examples/ADGPRS/5spot/5SPOT.gprs";
static std::string gprs_smry_json_5spot_     = base_path() + "/examples/ADGPRS/5spot/5SPOT.json";
static std::string gprs_smry_hdf5_5spot_     = base_path() + "/examples/ADGPRS/5spot/5SPOT.vars.h5";
static std::string gprs_base_5spot_          = base_path() + "/examples/ADGPRS/5spot/5SPOT";
static std::string gprs_smry_gpt_hdf5_5spot_ = base_path() + "/examples/ADGPRS/5spot/5SPOT-gpt.vars.h5";

static std::string trajectories_             = base_path() + "/examples/ECLIPSE/norne-simplified/trajectories";
static std::string cube_grid_                = base_path() + "/examples/ECLIPSE/cube_9x9/CUBE.EGRID";
static std::string schedule_inset_           = base_path() + "/examples/ECLIPSE/schedule_inset.txt";

static std::string olympr37_dir_    = base_path() + "/examples/ECLIPSE/olympr37";
static std::string olympr37_egrid_  = olympr37_dir_ + "/OLPZ_BCXX_R37_F37_W01.EGRID";
static std::string olympr37_driver_ = olympr37_dir_ + "/fo-drv.r001.c02-olymr37.icd-opt.apps.aug.json";

}
}

#endif //FIELDOPT_TEST_RESOURCE_EXAMPLE_FILE_PATHS_H
