/***********************************************************
Created by einar on 4/4/17.
Copyright (C) 2017
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

#ifndef FIELDOPT_TEST_RESOURCE_RUNNER_HPP_H
#define FIELDOPT_TEST_RESOURCE_RUNNER_HPP_H

#include "../runtime_settings.h"
#include "../logger.h"
#include "../../Settings/tests/test_resource_example_file_paths.hpp"

namespace TestResources {

using namespace TestResources::ExampleFilePaths;

class RunnerResources {

 private:

  // -----------------------------------------------------
  // string olympr37_RLZ_driver_nm = "fo-drv.r001.c02-olympr37.icd-opt.apps.aug";
  // string olympr37_RLZ_dir = "../examples/ECLIPSE/olympr37-dwf-rlz-flat/BC01/W01/";
  // string olympr37_RLZ_driver = olympr37_RLZ_dir + "../drivers/" + olympr37_RLZ_driver_nm + ".json";
  //
  // string olympr37_RLZ_simexe = olympr37_RLZ_dir + "../drivers/bash_e300_slbs.sh";
  // string olympr37_RLZ_xdir = "../fieldopt-output/olympr37/" + olympr37_RLZ_driver_nm;
  //
  // string olympr37_RLZ_egrid = olympr37_RLZ_dir + "/OLPZ_BCXX_R37_F37_W01.EGRID";
  // string olympr37_RLZ_data = olympr37_RLZ_dir + "/OLPZ_BCXX_R37_F37_W01.DATA";
  //
  // string olympr37_RLZ_ensf = olympr37_RLZ_dir + "../drivers/olympr37.ens";
  //
  // int argc_olympr37 = 20;
  // const char *argv_olympr37[20] = {
  //   "FieldOpt", olympr37_RLZ_driver.c_str(),
  //   olympr37_RLZ_xdir.c_str(),
  //   "-r", "serial",
  //   "-g", olympr37_RLZ_egrid.c_str(),
  //   "-s", olympr37_RLZ_data.c_str(),
  //   "-f",
  //   "-e", olympr37_RLZ_simexe.c_str(),
  //   "--ensemble-path", olympr37_RLZ_ensf.c_str(),
  //   "-v", "0",
  //   "-n", "1", // threads-per-sim
  //   "-t", "1000", // sim-timeout
  // };

  // -----------------------------------------------------
  string olympr37_T01_driver_nm = "fo-drv.r001.c02-olympr37.icd-opt.apps.aug";
  string olympr37_T01_dir = "../examples/ECLIPSE/olympr37-dwf-rlz-flat/BC01/T01/";
  string olympr37_T01_driver = olympr37_T01_dir + "../../drivers/" + olympr37_T01_driver_nm + ".json";

  string olympr37_T01_simexe = olympr37_T01_dir + "../../drivers/bash_e300_slbs.sh";
  string olympr37_T01_xdir = "../fieldopt-output/olympr37/" + olympr37_T01_driver_nm;

  string olympr37_T01_egrid = olympr37_T01_dir + "/OLPZ_BCXX_R37_F37_W01.EGRID";
  string olympr37_T01_data = olympr37_T01_dir + "/OLPZ_BCXX_R37_F37_W01.DATA";

  int argc_olympr37_T01 = 18;
  const char *argv_olympr37_T01[18] = {
    "FieldOpt", olympr37_T01_driver.c_str(),
    olympr37_T01_xdir.c_str(),
    "-r", "serial",
    "-g", olympr37_T01_egrid.c_str(),
    "-s", olympr37_T01_data.c_str(),
    "-f",
    "-e", olympr37_T01_simexe.c_str(),
    "-v", "0",
    "-n", "1", // threads-per-sim
    "-t", "1000", // sim-timeout
  };

  // -----------------------------------------------------
  // string driver_nm = "fo-drv.r001.c04-5spot-flow.wbhp-opt.apps-t01.npv";
  string driver_nm = "fo-drv.r001.c04-5spot-flow.wbhp-opt.trdfo-t01.npv";
  string en_5spot_dir = "../examples/ECLIPSE/5spot_c04_flow_en10_perm";
  string en_5spot_driver = en_5spot_dir + "/drivers/" + driver_nm + ".json";

  string en_5spot_simexe = en_5spot_dir + "/drivers/bash_flw_en_5spot.sh";
  // string en_5spot_simexe = en_5spot_dir + "/drivers/bash_ecl_en_5spot.sh";
  // string en_5spot_simexe = en_5spot_dir + "/drivers/bash_ecl_en_5spot_xe.sh";

  string en_5spot_xdir = "../fieldopt-output/en_5spot/" + driver_nm;
  string en_5spot_aux = en_5spot_dir + "/include";

  string en_5spot_egrid = en_5spot_dir +
    "/r001/ECL_5SPOT_C04_FLOW__R001.EGRID";
  string en_5spot_data = en_5spot_dir +
    "/r001/ECL_5SPOT_C04_FLOW__R001.DATA";

  // string en_5spot_ensf = en_5spot_dir + "/en_5spot_dbg.ens";
  string en_5spot_ensf = en_5spot_dir + "/en_5spot.ens";

  int argc_en_5spot = 22;
  const char *argv_en_5spot[22] = {
    "FieldOpt", en_5spot_driver.c_str(),
    en_5spot_xdir.c_str(),
    "-r", "serial",
    "-g", en_5spot_egrid.c_str(),
    "-s", en_5spot_data.c_str(),
    "-f",
    "-e", en_5spot_simexe.c_str(),
    "--sim-aux", en_5spot_aux.c_str(),
    "--ensemble-path", en_5spot_ensf.c_str(),
    "-v", "0",
    "-n", "1", // threads-per-sim
    "-t", "1000", // sim-timeout
  };

  // -----------------------------------------------------
  const int argc = 16;
  const char *argv[16] = {
    "FieldOpt",
    "../examples/ADGPRS/5spot/fo_driver_5vert_wells.json",
    "../fieldopt-output",
    "-g", "../examples/Flow/5spot/5SPOT.EGRID",
    "-s", "../examples/Flow/5spot/5SPOT.DATA",
    "-b", ".",
    "-r", "mpisync",
    "-f",
    "-v", "0",
    "-t", "1000"
  };

 public:

  Runner::RuntimeSettings *rts_ =
    new Runner::RuntimeSettings(argc, argv);

  Runner::RuntimeSettings *rts_en_5spot_ =
    new Runner::RuntimeSettings(argc_en_5spot, argv_en_5spot);

  Runner::RuntimeSettings *rts_olympr37_T01_ =
    new Runner::RuntimeSettings(argc_olympr37_T01, argv_olympr37_T01);

 protected:
  Logger *logger_ = new Logger(rts_,
                               QString::fromStdString(directory_output_),
                               false);

  Logger *logger_norne_ = new Logger(rts_,
                                     QString::fromStdString(norne_test_output_),
                                     false);

  Logger *logger_olympr37_T01_ = new Logger(rts_olympr37_T01_,
                                     QString::fromStdString(olympr37_T01_test_out_),
                                     false);

};

}

#endif //FIELDOPT_TEST_RESOURCE_RUNNER_HPP_H
