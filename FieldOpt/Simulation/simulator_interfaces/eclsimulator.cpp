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

#include <iostream>
#include <boost/algorithm/string.hpp>
#include <Utilities/printer.hpp>
// #include <Utilities/verbosity.h>

#include "eclsimulator.h"
#include "Utilities/execution.hpp"
#include "simulator_exceptions.h"
#include "Simulation/results/eclresults.h"

namespace Simulation {

using namespace Utilities::FileHandling;
using Utilities::Unix::ExecShellScript;
using Printer::ext_info;
using Printer::info;

ECLSimulator::ECLSimulator(Settings::Settings *settings,
                           Model::Model *model)
  : Simulator(settings) {
  model_ = model;
  settings_ = settings;
  vp_ = settings_->global()->verbParams();

  if (!paths_.IsSet(Paths::ENSEMBLE_FILE)) {
    deck_name_ = driver_file_name_.split(".DATA").first();
  } else {
    deck_name_ = "";
  }

  results_ = new Results::ECLResults(settings);
}

void ECLSimulator::EvalSecndOptmzr() {
  if (vp_.vSIM >= 2) {
    im_ = "Starting EvalSecndOptmzr [DFTR].";
    ext_info(im_,md_, cl_, vp_.lnw);
  }
}

void ECLSimulator::EvalXWrapper() {
  if (vp_.vSIM >= 2) {
    im_ = "Starting EvalXWrapper [sqp_snopt].";
    ext_info(im_,md_, cl_, vp_.lnw);
  }
  copyDriverFiles();
  UpdateFilePaths();
  script_args_ = (QStringList() << paths_.GetPathQstr(Paths::SIM_WORK_DIR) << deck_name_);

  auto driver_file_writer = EclDriverFileWriter(settings_, model_);

}

void ECLSimulator::Evaluate() {

  if (vp_.vSIM >= 2) {
    im_ = "Starting unmonitored evaluation. Copying driver files.";
    ext_info(im_,md_, cl_, vp_.lnw);
  }
  copyDriverFiles();

  if (vp_.vSIM >= 2) { info("Updating file paths.", vp_.lnw); }
  UpdateFilePaths();
  script_args_ = (QStringList() << paths_.GetPathQstr(Paths::SIM_WORK_DIR) << deck_name_);

  auto driver_file_writer = EclDriverFileWriter(settings_, model_);

  if (vp_.vSIM >= 2) { info("Writing schedule.", vp_.lnw); }
  driver_file_writer.WriteDriverFile(
    paths_.GetPathQstr(Paths::SIM_OUT_SCH_FILE));

  PreSimWork();
  if (vp_.vSIM >= 2) { info("Starting unmonitored sim.", vp_.lnw); }
  ExecShellScript(paths_.GetPathQstr(Paths::SIM_EXEC_SCRIPT_FILE),
                  script_args_,vp_);
  PostSimWork();

  if (vp_.vSIM >= 2) { info("Unmonitored sim done. Reading results.", vp_.lnw); }
  results_->ReadResults(paths_.GetPathQstr(Paths::SIM_OUT_DRIVER_FILE));
  updateResultsInModel();
}

bool ECLSimulator::Evaluate(int timeout, int threads=1) {
  copyDriverFiles();
  UpdateFilePaths();

  auto sim_work = paths_.GetPathQstr(Paths::SIM_WORK_DIR);
  script_args_ = (QStringList() << sim_work << deck_name_ << QString::number(threads));

  auto driver_file_writer = EclDriverFileWriter(settings_, model_);
  auto sch_file_path = paths_.GetPathQstr(Paths::SIM_OUT_SCH_FILE);
  driver_file_writer.WriteDriverFile(sch_file_path);

  int t = timeout; // Let sims run for at least 10 secs
  if ( timeout < 10 ) { t = 10; }
  PreSimWork();

  if (vp_.vSIM >= 2) { info("Starting monitored sim w/ timeout.", vp_.lnw); }
  auto script_path = paths_.GetPathQstr(Paths::SIM_EXEC_SCRIPT_FILE);
  bool success = ::Utilities::Unix::ExecShellScriptTimeout(script_path, script_args_, t, vp_);

  if (vp_.vSIM >= 2) { info("Monitored sim done.", vp_.lnw); }
  if ( success ) {
    if (vp_.vSIM >= 2) { info("Sim successful. Reading results.", vp_.lnw); }
    results_->DumpResults();
    PostSimWork();

    auto file_path = paths_.GetPathQstr(Paths::SIM_OUT_DRIVER_FILE);
    results_->ReadResults(file_path);
  }
  updateResultsInModel();
  return success;
}

bool ECLSimulator::
Evaluate(const Settings::Ensemble::Realization &realization,
         int timeout, int threads=1) {

  driver_file_name_ =
    QString::fromStdString(FileName(realization.data()));
  driver_parent_dir_name_ =
    QString::fromStdString(ParentDirName(realization.data()));
  deck_name_ = driver_file_name_.split(".").first();

  paths_.SetPath(Paths::SIM_DRIVER_FILE, realization.data());
  paths_.SetPath(Paths::SIM_DRIVER_DIR, GetParentDirPath(realization.data()));
  paths_.SetPath(Paths::SIM_SCH_FILE, realization.schedule());

  // return Evaluate(timeout, threads); // Hack
  Evaluate();
  return true;
}

void ECLSimulator::CleanUp() {
  UpdateFilePaths();
  QStringList file_endings_to_delete{
    "DBG", "ECLEND", "ECLRUN", "EGRID", "GRID",
    "h5", "INIT", "INSPEC", "MSG", "PRT",
    "RSSPEC", "UNRST"};
  QString base_file_path = QString::fromStdString(
    paths_.GetPath(Paths::SIM_OUT_DRIVER_FILE)).split(".DATA").first();

  for ( QString ending : file_endings_to_delete ) {
    DeleteFile(base_file_path + "." + ending, vp_);
  }
}

void ECLSimulator::UpdateFilePaths() {
  std::string workdir = paths_.GetPath(Paths::OUTPUT_DIR);
  std::string casedir = workdir;
  workdir += "/" + driver_parent_dir_name_.toStdString();
  casedir += "/" + case_parent_dir_name_.toStdString();

  if (settings_->simulator()->file_structure().type_
    == Settings::Simulator::FileStructType::Flat) {

    paths_.SetPath(Paths::SIM_WORK_DIR, workdir);
    paths_.SetPath(Paths::SIM_OUT_DRIVER_FILE,
                   paths_.GetPath(Paths::SIM_WORK_DIR) + "/" + driver_file_name_.toStdString());
  }

  if (settings_->simulator()->file_structure().type_
    == Settings::Simulator::FileStructType::Branched) {

    paths_.SetPath(Paths::SIM_WORK_DIR,
                   casedir + paths_.GetPath(Paths::CASE_DRVR_DIR));
    paths_.SetPath(Paths::SIM_OUT_DRIVER_FILE,
                   paths_.GetPath(Paths::SIM_WORK_DIR) + "/" + driver_file_name_.toStdString());
  }

  std::string tmp = paths_.GetPath(Paths::SIM_SCH_FILE);
  boost::algorithm::replace_first(tmp,
                                  paths_.GetPath(Paths::SIM_DRIVER_DIR),
                                  paths_.GetPath(Paths::SIM_WORK_DIR));
  paths_.SetPath(Paths::SIM_OUT_SCH_FILE, tmp);
}

void ECLSimulator::WriteDriverFilesOnly() {
  UpdateFilePaths();
  auto driver_file_writer = EclDriverFileWriter(settings_, model_);
  driver_file_writer.WriteDriverFile(paths_.GetPathQstr(Paths::SIM_OUT_SCH_FILE));
}

using FT = Settings::Simulator::FileStructType;

void ECLSimulator::copyDriverFiles() {
  std::string workdir = paths_.GetPath(Paths::OUTPUT_DIR);
  std::string casedir = workdir;
  workdir += "/" + driver_parent_dir_name_.toStdString();
  casedir += "/" + case_parent_dir_name_.toStdString();

  // Flat file structure
  if (!DirExists(workdir, vp_, md_, cl_)
    && settings_->simulator()->file_structure().type_ == FT::Flat) {
    if (vp_.vSIM >= 1) {
      im_ = "Output deck directory not found. Copying input deck:";
      im_ += paths_.GetPath(Paths::SIM_DRIVER_DIR) + " -> " + workdir;
      ext_info(im_,md_,cl_);
    }
    CreateDir(workdir, vp_);
    CopyDir(paths_.GetPath(Paths::SIM_DRIVER_DIR), workdir, false, false, vp_);
    paths_.SetPath(Paths::SIM_WORK_DIR, workdir);
  }

  // Branched file structure
  if (!DirExists(casedir, vp_, md_, cl_)
    && settings_->simulator()->file_structure().type_ == FT::Branched) {
    if (vp_.vSIM >= 1) {
      im_ = "Output dir not found. Copying input deck:";
      im_ += paths_.GetPath(Paths::CASE_ROOT_DIR) + " -> " + casedir;
      ext_info(im_, md_, cl_);
    }
    CreateDir(casedir, vp_);
    CopyDir(paths_.GetPath(Paths::CASE_ROOT_DIR), casedir, false, true, vp_);
    paths_.SetPath(Paths::SIM_WORK_DIR, casedir + paths_.GetPath(Paths::CASE_DRVR_DIR));
  }

  // Copy sim aux dir
  if ( paths_.IsSet(Paths::SIM_AUX_DIR)) {
    std::string auxdir = paths_.GetPath(Paths::OUTPUT_DIR)
      + "/" + FileName(paths_.GetPath(Paths::SIM_AUX_DIR));

    if ( !DirExists(auxdir, vp_, md_, cl_)) {
      if (vp_.vSIM >= 1) {
        im_ = "Copying simulation aux. directory:";
        im_ += paths_.GetPath(Paths::SIM_AUX_DIR) + " -> " + auxdir;
        ext_info(im_, md_, cl_ );
      }
      CreateDir(auxdir, vp_);
      CopyDir(paths_.GetPath(Paths::SIM_AUX_DIR), auxdir, false);
    }
  }

  // Copy pre-/postsim scripts to workdir
  if (settings_->simulator()->add_sim_scripts()) {
    string post_script_src = paths_.GetPath(Paths::SIM_EXEC_DIR) + "/FO_POSTSIM.sh";
    string pre_script_src = paths_.GetPath(Paths::SIM_EXEC_DIR) + "/FO_PRESIM.sh";

    string post_script_trg = paths_.GetPath(Paths::SIM_WORK_DIR) + "/FO_POSTSIM.sh";
    string pre_script_trg = paths_.GetPath(Paths::SIM_WORK_DIR) + "/FO_PRESIM.sh";

    CopyFile(post_script_src, post_script_trg, true, vp_);
    CopyFile(pre_script_src, pre_script_trg, true, vp_);
  }

  if ( vp_.vSIM >= 2 ) {
    im_ = "Done copying directories. Set working directory to: " + workdir;
    ext_info(im_, md_, cl_);
  }
}

}
