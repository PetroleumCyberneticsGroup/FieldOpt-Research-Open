//
// Created by einar on 6/2/16.
//

#ifndef FIELDOPT_TEST_RESOURCE_EXAMPLE_FILE_PATHS_H
#define FIELDOPT_TEST_RESOURCE_EXAMPLE_FILE_PATHS_H

#include "Utilities/system.hpp"
#include <string>

namespace TestResources {
namespace ExampleFilePaths {

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
static std::string directory_output_icds_    = base_path() + "/fieldopt-output/5spot-ICD";
static std::string dir_5spot_SLB_output_     = base_path() + "/fieldopt-output/5spot-SLB";
static std::string norne_test_output_        = base_path() + "/fieldopt-output/norne_test";
static std::string norne_driver_example_     = base_path() + "/examples/ECLIPSE/norne-simplified/fo_driver.1_rate.cs.json";
static std::string norne_grid_               = base_path() + "/examples/ECLIPSE/norne-simplified/NORNE_SIMPLIFIED.EGRID";
static std::string norne_atw_grid_           = base_path() + "/examples/Flow/norne/NORNE_ATW2013.EGRID";
static std::string norne_deck_               = base_path() + "/examples/ECLIPSE/norne-simplified/NORNE_SIMPLIFIED.DATA";
static std::string norne_sch_                = base_path() + "/examples/ECLIPSE/norne-simplified/INCLUDE/BC0407_HIST01122006.SCH";
static std::string deck_horzwel_             = base_path() + "/examples/ECLIPSE/HORZWELL/HORZWELL.DATA";
static std::string grid_horzwel_             = base_path() + "/examples/ECLIPSE/HORZWELL/HORZWELL.EGRID";
static std::string ecl_base_horzwell         = base_path() + "/examples/ECLIPSE/HORZWELL/HORZWELL";

static std::string deck_ECL_5spot_           = base_path() + "/examples/ECLIPSE/5spot/ECL_5SPOT.DATA";
static std::string grid_ECL_5spot_           = base_path() + "/examples/ECLIPSE/5spot/ECL_5SPOT.EGRID";
static std::string driver_ECL_5pot_          = base_path() + "/examples/ECLIPSE/5spot/fo_driver_2_horz_icd_apps.json";

static std::string deck_5spot_icds_          = base_path() + "/examples/ECLIPSE/5spot-ICD/r0/ECL_5SPOT_C05_OPT.DATA";
static std::string grid_5spot_icds_          = base_path() + "/examples/ECLIPSE/5spot-ICD/r0/ECL_5SPOT_C05_OPT.DATA__R001.EGRID";
static std::string driver_5pot_icds_         = base_path() + "/examples/ECLIPSE/5spot-ICD/r0/fo_driver_2_horz_icd_apps.json";

static std::string deck_5spot_SLB            = base_path() + "/examples/ECLIPSE/5spot-SLB/mods/5spot_c04_flow_en10_perm/r001/ECL_5SPOT_C04_FLOW__R001.DATA ";
static std::string grid_5spot_SLB            = base_path() + "/examples/ECLIPSE/5spot-SLB/mods/5spot_c04_flow_en10_perm/r001/ECL_5SPOT_C04_FLOW__R001.EGRID";
static std::string driver_5spot_SLB          = base_path() + "/examples/ECLIPSE/5spot-SLB/mods/5spot_c04_flow_en10_perm/r001/fo-drv.r001.c04-5spot-flow.wbhp-opt.apps-t01.npv.json";
static std::string aux_dir_5spot_SLB         = base_path() + "/examples/ECLIPSE/5spot-SLB/mods/5spot_c04_flow";

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

}
}

#endif //FIELDOPT_TEST_RESOURCE_EXAMPLE_FILE_PATHS_H
