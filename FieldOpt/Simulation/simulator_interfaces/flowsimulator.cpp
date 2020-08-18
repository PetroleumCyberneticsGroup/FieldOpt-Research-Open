/***********************************************************
Copyright (C) 2015-2017
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

#include <results/eclresults.h>
#include <Utilities/execution.hpp>
#include <iostream>
#include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>
#include "flowsimulator.h"
#include "simulator_exceptions.h"

namespace Simulation {

FlowSimulator::FlowSimulator(Settings::Settings *settings, Model::Model *model)
  : Simulator(settings) {
  model_ = model;
  driver_file_writer_ = new FlowDriverFileWriter(settings_, model_);
  vp_ = settings_->global()->verbParams();

  verifyOriginalDriverFileDirectory();

  results_ = new Results::ECLResults(settings->simulator());
  try {
    results()->ReadResults(driver_file_writer_->output_driver_file_name_);
  } catch (...) {} // At this stage we don't really care if the results can be read, we just want to set the path.
}

void FlowSimulator::Evaluate() {
  if (results_->isAvailable()) results_->DumpResults();
  copyDriverFiles();
  driver_file_writer_->WriteDriverFile(QString::fromStdString(paths_.GetPath(Paths::SIM_WORK_DIR )));
  ::Utilities::Unix::ExecShellScript(QString::fromStdString(paths_.GetPath(Paths::SIM_EXEC_SCRIPT_FILE)), script_args_, vp_);
  results_->ReadResults(driver_file_writer_->output_driver_file_name_);
}

void FlowSimulator::CleanUp() {
  QStringList file_endings_to_delete{"EGRID", "INIT", "PRT", "RFT",
                                     "RSSPEC", "UNRST"};
  QString base_file_path = driver_file_writer_->output_driver_file_name_.split(".DATA").first();
  for (QString ending : file_endings_to_delete) {
    Utilities::FileHandling::DeleteFile(base_file_path + "." + ending, vp_);
  }
}

void FlowSimulator::verifyOriginalDriverFileDirectory() {
  QStringList critical_files = {QString::fromStdString(paths_.GetPath(Paths::SIM_DRIVER_DIR)) + "/include/compdat.in",
                                QString::fromStdString(paths_.GetPath(Paths::SIM_DRIVER_DIR)) + "/include/controls.in",
                                QString::fromStdString(paths_.GetPath(Paths::SIM_DRIVER_DIR)) + "/include/wells.in",
                                QString::fromStdString(paths_.GetPath(Paths::SIM_DRIVER_DIR)) + "/include/welspecs.in"};
  for (auto file : critical_files) {
    if (!Utilities::FileHandling::FileExists(file, vp_))
      throw DriverFileDoesNotExistException(file);
  }
}

void FlowSimulator::copyDriverFiles() {
  auto workdir = paths_.GetPath(Paths::OUTPUT_DIR) + driver_parent_dir_name_.toStdString();
  if (!DirExists(workdir, vp_, md_, cl_)) {
    if (VERB_SIM >= 1) {
      Printer::ext_info("Output deck directory not found. Copying input deck:"
                          + paths_.GetPath(Paths::SIM_DRIVER_DIR) + " -> " + workdir, "Simulation", "FlowSimulator" );
    }
    CreateDir(workdir, vp_);
    CopyDir(paths_.GetPath(Paths::SIM_DRIVER_DIR), workdir, true);
  }
  paths_.SetPath(Paths::SIM_WORK_DIR, workdir);
}

void FlowSimulator::UpdateFilePaths() {
  script_args_ = (QStringList() << QString::fromStdString(paths_.GetPath(Paths::SIM_WORK_DIR ))
                                << QString::fromStdString(paths_.GetPath(Paths::SIM_WORK_DIR ))
                                  + "/" + driver_file_name_);
}

bool FlowSimulator::Evaluate(int timeout, int threads) {
  int t = timeout;
  script_args_[2] = QString::number(threads);
  if (timeout < 10) t = 10; // Always let simulations run for at least 10 seconds
  if (results_->isAvailable()) results()->DumpResults();
  copyDriverFiles();
  driver_file_writer_->WriteDriverFile(QString::fromStdString(paths_.GetPath(Paths::SIM_WORK_DIR)));
  std::cout << "Starting monitored simulation with timeout " << timeout << std::endl;
  bool success = ::Utilities::Unix::ExecShellScriptTimeout(
    QString::fromStdString(paths_.GetPath(Paths::SIM_EXEC_SCRIPT_FILE)),
    script_args_, t, vp_);
  if (success) {
    results_->ReadResults(driver_file_writer_->output_driver_file_name_);
  }
  return success;
}
void FlowSimulator::WriteDriverFilesOnly() {
  copyDriverFiles();
  driver_file_writer_->WriteDriverFile(paths_.GetPathQstr(Paths::SIM_WORK_DIR ));
}
bool FlowSimulator::Evaluate(const Settings::Ensemble::Realization &realization, int timeout, int threads) {
  throw std::runtime_error("Ensemble optimization not yet implemented for the FLOW reservoir simulator.");
}
}
