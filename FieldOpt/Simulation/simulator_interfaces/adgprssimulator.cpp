/***********************************************************
Copyright (C) 2015-2017
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

#include <iostream>
// #include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>
#include "adgprssimulator.h"
#include "Utilities/execution.hpp"
#include "simulator_exceptions.h"
#include "Simulation/results/adgprsresults.h"

namespace Simulation {

AdgprsSimulator::AdgprsSimulator(Settings::Settings *settings,
                                 Model::Model *model)
  : Simulator(settings) {
  model_ = model;
  settings_ = settings;
  vp_ = settings_->global()->verbParams();

  driver_file_writer_ = new AdgprsDriverFileWriter(settings_,
                                                   model_);
  verifyOriginalDriverFileDirectory();
  results_ = new Simulation::Results::AdgprsResults(settings);
}

void AdgprsSimulator::Evaluate() {
  if (results_->isAvailable()) results()->DumpResults();
  copyDriverFiles();
  driver_file_writer_->WriteDriverFile(paths_.GetPathQstr(Paths::SIM_WORK_DIR));
  ::Utilities::Unix::ExecShellScript(paths_.GetPathQstr(Paths::SIM_EXEC_SCRIPT_FILE), script_args_, vp_);
  paths_.SetPath(Paths::SIM_HDF5_FILE,
                 paths_.GetPath(Paths::SIM_WORK_DIR) + "/"
                   + driver_file_name_.split(".").first().toStdString() + ".vars.h5"
  );
  results_->ReadResults(paths_.GetPathQstr(Paths::SIM_HDF5_FILE));
  updateResultsInModel();
}

void AdgprsSimulator::CleanUp() {
  DeleteFile(QString::fromStdString(paths_.GetPath(Paths::SIM_HDF5_FILE)), vp_);
}

void AdgprsSimulator::copyDriverFiles() {
  auto workdir = paths_.GetPath(Paths::OUTPUT_DIR) + driver_parent_dir_name_.toStdString();
  if (!DirExists(workdir, vp_)) {
    if (vp_.vSIM >= 1) {
      string im = "Output deck directory not found. Copying input deck:";
      im += paths_.GetPath(Paths::SIM_DRIVER_DIR) + " -> " + workdir;
      ext_info(im, "Simulation", "ADGPRSSimulator" );
    }
    CreateDir(workdir, vp_);
    CopyDir(paths_.GetPath(Paths::SIM_DRIVER_DIR), workdir, true);
  }
  paths_.SetPath(Paths::SIM_WORK_DIR, workdir);
}

void AdgprsSimulator::verifyOriginalDriverFileDirectory() {
  QStringList critical_files = {paths_.GetPathQstr(Paths::SIM_DRIVER_DIR) + "/include/compdat.in",
                                paths_.GetPathQstr(Paths::SIM_DRIVER_DIR) + "/include/controls.in",
                                paths_.GetPathQstr(Paths::SIM_DRIVER_DIR) + "/include/wells.in",
                                paths_.GetPathQstr(Paths::SIM_DRIVER_DIR) + "/include/welspecs.in"};
  for (auto file : critical_files) {
    if (!Utilities::FileHandling::FileExists(file, vp_))
      throw DriverFileDoesNotExistException(file);
  }
}

void AdgprsSimulator::UpdateFilePaths() {
  script_args_ = (QStringList() << QString::fromStdString(paths_.GetPath(Paths::SIM_WORK_DIR ))
                                << QString::fromStdString(paths_.GetPath(Paths::SIM_WORK_DIR ))
                                  + "/" + driver_file_name_);
}

bool AdgprsSimulator::Evaluate(int timeout, int threads) {
  script_args_[2] = QString::number(threads);
  int t = timeout;
  if (timeout < 10) t = 10; // Always let simulations run for at least 10 seconds
  if (results_->isAvailable()) results()->DumpResults();
  copyDriverFiles();
  driver_file_writer_->WriteDriverFile(QString::fromStdString(paths_.GetPath(Paths::SIM_WORK_DIR )));
  std::cout << "Starting monitored simulation with timeout " << timeout << std::endl;
  bool success = ::Utilities::Unix::ExecShellScriptTimeout(
    QString::fromStdString(paths_.GetPath(Paths::SIM_EXEC_SCRIPT_FILE)),
    script_args_, t, vp_);
  if (success) {
    paths_.SetPath(Paths::SIM_HDF5_FILE,
                   paths_.GetPath(Paths::SIM_WORK_DIR) + "/"
                     + driver_file_name_.split(".").first().toStdString() + ".vars.h5"
    );
    results_->ReadResults(output_h5_summary_file_path_);
  }
  updateResultsInModel();
  return success;
}

void AdgprsSimulator::WriteDriverFilesOnly() {
  copyDriverFiles();
  driver_file_writer_->WriteDriverFile(QString::fromStdString(paths_.GetPath(Paths::SIM_WORK_DIR )));
}

bool AdgprsSimulator::Evaluate(const Settings::Ensemble::Realization &realization, int timeout, int threads) {
  throw std::runtime_error("Ensemble optimization not yet implemented for the AD-GPRS reservoir simulator.");
}

}
