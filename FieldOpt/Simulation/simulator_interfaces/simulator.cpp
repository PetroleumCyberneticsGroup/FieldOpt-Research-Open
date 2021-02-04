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

#include "simulator.h"
#include "simulator_exceptions.h"
#include "Utilities/execution.hpp"
#include "Simulation/results/json_results.h"

#include "QByteArray"
#include "QJsonDocument"

namespace Simulation {

using std::string;
using namespace Utilities::FileHandling;
using Utilities::FileHandling::FileExists;
//using Utilities::FileHandling::FileNameRoot;
using Utilities::Unix::ExecShellScript;

using Printer::info;
using Printer::ext_info;
using Printer::ext_warn;

Simulator::Simulator(Settings::Settings *settings) {
  settings_ = settings;
  vp_ = settings_->global()->verbParams();
  paths_ = settings_->paths();

  if (!paths_.IsSet(Paths::ENSEMBLE_FILE)) { // single realization
    driver_file_name_ = FileNameQstr(paths_.GetPath(Paths::SIM_DRIVER_FILE));
    driver_parent_dir_name_ = ParentDirectoryNameQstr(paths_.GetPath(Paths::SIM_DRIVER_FILE));

    if (paths_.IsSet(Paths::CASE_ROOT_DIR)) {
      case_parent_dir_name_ = ParentDirectoryNameQstr(paths_.GetPath(Paths::CASE_ROOT_DIR));
    }

  } else { // multiple realizations
    driver_file_name_ = "";
    driver_parent_dir_name_ = "";
  }

  // Use custom execution script if provided in runtime
  // settings, else use the one from json driver file
  if (!paths_.IsSet(Paths::SIM_EXEC_SCRIPT_FILE)) {
    string exec_script_path = paths_.GetPath(Paths::BUILD_DIR);
    auto sn = settings->simulator()->script_name();
    if (sn.isEmpty()) {
      E("No execution script given in JSON.");
    } else {
      exec_script_path += ExecutionScripts::GetScriptPath(sn).toStdString();
      paths_.SetPath(Paths::SIM_EXEC_SCRIPT_FILE, exec_script_path);
      paths_.SetPath(Paths::SIM_EXEC_DIR, GetParentDirPath(exec_script_path));
    }
  }

  script_args_ = (QStringList() << paths_.GetPathQstr(Paths::OUTPUT_DIR)
                                << paths_.GetPathQstr(Paths::OUTPUT_DIR) + "/" + driver_file_name_
                                << QString::number(1));

  control_times_ = settings->model()->control_times();
}

::Simulation::Results::Results *Simulator::results() {
  return results_;
}

void Simulator::PreSimWork() {
  sim_wrk_dir_ = GetAbsoluteFilePath(paths_.GetPathQstr(Paths::SIM_WORK_DIR));
  QString presim_script = paths_.GetPathQstr(Paths::SIM_WORK_DIR) + "/FO_PRESIM.sh";
  QStringList* presim_args = settings_->simulator()->pre_sim_args();
  presim_args->prepend(sim_wrk_dir_);

  if (settings_->simulator()->use_pre_sim_script() && FileExists(presim_script, vp_, md_, cl_)) {
    ExecShellScript(presim_script, *presim_args, vp_);
  } else { ext_warn("Presim script not found."); }
}

void Simulator::PostSimWork() {
  ::Simulation::Results::Results::EclAjdGData adgG;

  if (settings_->simulator()->use_post_sim_script()) {

    ext_info("Executing postsim script", md_, cl_, vp_.lnw);
    QString expected_script_path = paths_.GetPathQstr(Paths::SIM_WORK_DIR) + "/FO_POSTSIM.sh";
    QStringList post_sim_args = *settings_->simulator()->post_sim_args();
    post_sim_args.prepend(sim_wrk_dir_);

    if (vp_.vSIM >= 3) {
      stringstream ss; ss << "post_sim_args: ";
      for (int ii=0; ii < post_sim_args.size(); ++ii) {
        ss << post_sim_args.at(ii).toStdString() + " ";
      }
    }

    if (FileExists(expected_script_path, vp_, md_, cl_)) { //  && post_sim_args->size() > 2
      ExecShellScript(expected_script_path, post_sim_args, vp_);
    } else {
      ext_warn("Postsim script not found.", md_, cl_, vp_.lnw);
    }
  }

  if (settings_->simulator()->read_external_json_results()) {
    std::string ext_file = paths_.GetPath(Paths::SIM_WORK_DIR) + "/FO_EXT_RESULTS.json";

    if (vp_.vSIM >= 2) { ext_info("Reading external JSON results: " + ext_file, md_, cl_); }
    auto json_results = Simulation::Results::JsonResults(ext_file, *settings_->simulator());
    results_->SetJsonResults(json_results);
  }

  if(settings_->simulator()->read_adj_grad_data()) {
    string adjg_json = sim_wrk_dir_.toStdString() + "/" + FileNameRoot(driver_file_name_.toStdString()) + ".JSON";

    if (vp_.vSIM >= 2) { ext_info("Reading adj grad data: " + adjg_json, md_, cl_); }
    auto grad_data = Simulation::Results::JsonResults(adjg_json, *settings_->simulator());
    results_->SetJsonResults(grad_data);
  }
}

// Check if obsolete: Reading e300 gradients done by PostSimWork:
//::Simulation::Results::Results::EclAjdGData Simulator::ReadEclAdjG() {
//  QString data_root = FileNameQstr(paths_.GetPath(Paths::SIM_OUT_DRIVER_FILE));
//  string adjg_file = data_root.split(".DATA").first().toStdString() + ".GRD";
//
//  string tm = "Reading gradient from e300 file: " + adjg_file;
//  ext_info(tm, md_, cl_, vp_.lnw);
//}

void Simulator::updateResultsInModel() {
  model_->SetResult("Time", results_->GetValueVector(Results::Results::Property::Time));
  model_->SetResult("FGPT", results_->GetValueVector(Results::Results::Property::FieldGasProdTotal));
  model_->SetResult("FOPT", results_->GetValueVector(Results::Results::Property::FieldOilProdTotal));
  model_->SetResult("FWPT", results_->GetValueVector(Results::Results::Property::FieldWatProdTotal));
}

}
