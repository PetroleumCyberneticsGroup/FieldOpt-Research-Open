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
    SIM_HDF5_FILE=7, ENSEMBLE_FILE=8, SIM_SCH_INSET_FILE=9,
    BUILD_DIR = -1, OUTPUT_DIR = -2, SIM_DRIVER_DIR = -3,
    SIM_WORK_DIR = -4, SIM_AUX_DIR = -5, TRAJ_DIR = -6,
    CASE_ROOT_DIR = -7, CASE_DRVR_DIR = -8,
    SIM_EXEC_DIR = -9
  };

  const string &GetPathDesc(Path path) const;

  void SetPath(Path path, const string& path_string,
               bool skip_check = false,
               Settings::VerbParams vp={});

  bool IsSet(Path path);

  string GetPath(Path path);
  QString GetPathQstr(Path path);

  void ShowPaths();

 private:
  map<Path, string> path_descriptions = {
    pair<Path, string> {DRIVER_FILE, "Driver file"},
    pair<Path, string> {SIM_DRIVER_FILE, "Sim driver file"},
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
    pair<Path, string> {SIM_DRIVER_DIR, "Sim driver parent directory"},
    pair<Path, string> {SIM_WORK_DIR, "Sim work directory"},
    pair<Path, string> {SIM_AUX_DIR, "Auxiliary files for simulation directory"},
    pair<Path, string> {TRAJ_DIR, "Dir w/ trajectory files for import"},
    pair<Path, string> {CASE_ROOT_DIR, "Case root dir"},
    pair<Path, string> {CASE_DRVR_DIR, "Case driver dir (relative to ROOT)"}
  };

  map<Path, string> paths_;

  string md_ = "Settings";
  string cl_ = "Paths";
  Settings::VerbParams vp_;
};

#endif //FIELDOPT_PATHS_H
