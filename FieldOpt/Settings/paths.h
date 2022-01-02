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

#ifndef FIELDOPT_PATHS_H
#define FIELDOPT_PATHS_H

#include "Utilities/filehandling.hpp"
#include "Utilities/verbosity.h"

using namespace Utilities::FileHandling;
using std::string;
using std::pair;
using std::map;

class Paths {
 public:
  Paths();

  /*!
   * Path enum. Directory paths should have values < 0; files >= 0.
   */
  enum Path : int {
    DRIVER_FILE = 0, SIM_DRIVER_FILE = 1, GRID_FILE = 2,
    SIM_EXEC_SCRIPT_FILE = 3, SIM_SCH_FILE = 4,
    SIM_OUT_DRIVER_FILE = 5, SIM_OUT_SCH_FILE = 6,
    SIM_HDF5_FILE = 7, ENSEMBLE_FILE = 8,
    SIM_SCH_INSET_FILE = 9, RESTART_FILE = 10,
    BUILD_DIR = -1, OUTPUT_DIR = -2, SIM_DRIVER_DIR = -3,
    SIM_WORK_DIR = -4, SIM_AUX_DIR = -5, TRAJ_DIR = -6,
    CASE_ROOT_DIR = -7, CASE_DRVR_DIR = -8,
    SIM_EXEC_DIR = -9, OPTMZD_DIR = -10
  };

  const string &GetPathDesc(Path path) const;

  void SetPath(Path path, const string& path_string,
               bool skip_check = false,
               Settings::VerbParams vp={});

  void CopyPath(Path path0, Path path1);

  bool IsSet(Path path);

  string GetPath(Path path);
  QString GetPathQstr(Path path);

  void ShowPaths();

 private:
  map<Path, string> path_descriptions = {
    pair<Path, string> {DRIVER_FILE, "Driver file"},
    pair<Path, string> {SIM_DRIVER_FILE, "Sim driver file"},
    pair<Path, string> {RESTART_FILE, "restart file"},
    pair<Path, string> {GRID_FILE, "Grid file"},
    pair<Path, string> {SIM_EXEC_SCRIPT_FILE, "Sim exe script"},
    pair<Path, string> {SIM_SCH_FILE, "Sim schedule file"},
    pair<Path, string> {SIM_OUT_DRIVER_FILE, "Sim output driver file"},
    pair<Path, string> {SIM_OUT_SCH_FILE, "Sim output schedule file"},
    pair<Path, string> {SIM_SCH_INSET_FILE, "Sim schedule inset file"},
    pair<Path, string> {SIM_HDF5_FILE, "HDF5 Summary file"},
    pair<Path, string> {ENSEMBLE_FILE, "Ensemble description file"},
    pair<Path, string> {BUILD_DIR, "Build directory"},
    pair<Path, string> {OUTPUT_DIR, "Output directory"},
    pair<Path, string> {OPTMZD_DIR, "Optimized result dir"},
    pair<Path, string> {SIM_EXEC_DIR, "Sim exec dir"},
    pair<Path, string> {SIM_DRIVER_DIR, "Sim driver parent dir"},
    pair<Path, string> {SIM_WORK_DIR, "Sim work directory"},
    pair<Path, string> {SIM_AUX_DIR, "Auxiliary files for simulation dir"},
    pair<Path, string> {TRAJ_DIR, "Dir w/ trajectory files for import"},
    pair<Path, string> {CASE_ROOT_DIR, "Case root dir"},
    pair<Path, string> {CASE_DRVR_DIR, "Case driver dir (relative to ROOT)"}
  };

  // Examples:
  // Optimized result dir          : MAIN/xdir/20210513-5spot-icd-opt-sqp.bc01.xrun-sqp.t200/fo.5spot_c09_bc01.bc01.xrun-sqp.NPV-SQP-X00-A01_20210513-T121955/optcs
  // Sim exec dir                  : MAIN/xdir/20210513-5spot-icd-opt-sqp.bc01.xrun-sqp.t200/drivers
  // Case driver dir (relative to ROOT):
  // Sim driver parent dir         : MAIN/runs/20210210_e300.5spot.icd-opt-cong.test-case-dev-04/x.e300.c09_bc01-5spot_sc
  // Output directory              : MAIN/xdir/20210513-5spot-icd-opt-sqp.bc01.xrun-sqp.t200/fo.5spot_c09_bc01.bc01.xrun-sqp.NPV-SQP-X00-A01_20210513-T121955
  // Build directory               : MAIN/xdir/20210513-5spot-icd-opt-sqp.bc01.xrun-sqp.t200/drivers/
  // Driver file                   : MAIN/xdir/20210513-5spot-icd-opt-sqp.bc01.xrun-sqp.t200/drivers/fo-drv.bc01.xrun-sqp_NPV-SQP-X00-A01.json
  // Sim driver file               : MAIN/runs/20210210_e300.5spot.icd-opt-cong.test-case-dev-04/x.e300.c09_bc01-5spot_sc/ECL_5SPOT_C09_BC01_OPT.DATA
  // Grid file                     : MAIN/runs/20210210_e300.5spot.icd-opt-cong.test-case-dev-04/x.e300.c09_bc01-5spot_sc/ECL_5SPOT_C09_BC01_OPT.EGRID
  // Sim exe script                : MAIN/xdir/20210513-5spot-icd-opt-sqp.bc01.xrun-sqp.t200/drivers/bash_e300_slbs.sh
  // Sim schedule file             : MAIN/runs/20210210_e300.5spot.icd-opt-cong.test-case-dev-04/x.e300.c09_bc01-5spot_sc/fo_edits.INC

  map<Path, string> paths_;

  string md_ = "Settings";
  string cl_ = "Paths";
  Settings::VerbParams vp_;
};

#endif //FIELDOPT_PATHS_H
