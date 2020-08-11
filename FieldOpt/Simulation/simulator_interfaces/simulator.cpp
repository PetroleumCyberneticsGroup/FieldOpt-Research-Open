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

#include "simulator.h"
#include "simulator_exceptions.h"
#include "Utilities/execution.hpp"
#include "Simulation/results/json_results.h"

using std::string;

namespace Simulation {

using namespace Utilities::FileHandling;

Simulator::Simulator(Settings::Settings *settings) {
  settings_ = settings;
  paths_ = settings_->paths();

  if (!paths_.IsSet(Paths::ENSEMBLE_FILE)) { // single realization
    driver_file_name_ = FileNameQstr(paths_.GetPath(Paths::SIM_DRIVER_FILE));
    driver_parent_dir_name_ = ParentDirectoryNameQstr(paths_.GetPath(Paths::SIM_DRIVER_FILE));
    case_parent_dir_name_ = ParentDirectoryNameQstr(paths_.GetPath(Paths::CASE_ROOT_DIR));

  } else { // multiple realizations
    driver_file_name_ = "";
    driver_parent_dir_name_ = "";
  }

  // Use custom execution script if provided in runtime
  // settings, else use the one from json driver file
  if (!paths_.IsSet(Paths::SIM_EXEC_SCRIPT_FILE)) {
    std::string exec_script_path = paths_.GetPath(Paths::BUILD_DIR)
      + ExecutionScripts::GetScriptPath(settings->simulator()->script_name()).toStdString();
    paths_.SetPath(Paths::SIM_EXEC_SCRIPT_FILE, exec_script_path);
    paths_.SetPath(Paths::SIM_EXEC_DIR, GetParentDirPath(exec_script_path));
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
  string presim_script = paths_.GetPath(Paths::SIM_WORK_DIR) + "/FO_PRESIM.sh";
  QStringList* presim_args = settings_->simulator()->pre_sim_args();

//  QString data_root = FileNameQstr(paths_.GetPath(Paths::SIM_OUT_DRIVER_FILE));
//  QString opt_inc_file = data_root.split(".DATA").first() + ".INC";
//  presim_args->append(data_root);

  QString sim_wrk_dir = FileNameQstr(paths_.GetPath(Paths::SIM_WORK_DIR));
  presim_args->prepend(sim_wrk_dir);

  if (settings_->simulator()->use_pre_sim_script() &&
    Utilities::FileHandling::FileExists(presim_script)) {

    Utilities::Unix::ExecShellScript(QString::fromStdString(presim_script),
                                     *presim_args);

  } else {
    Printer::ext_warn("Presim script not found.");
  }
}

::Simulation::Results::Results::EclAjdGData Simulator::ReadEclAdjG() {
  QString data_root = FileNameQstr(paths_.GetPath(Paths::SIM_OUT_DRIVER_FILE));
  string adjg_file = data_root.split(".DATA").first().toStdString() + ".GRD";

  cout << "adjg_file: " << adjg_file << endl;

}

void Simulator::PostSimWork() {

  ::Simulation::Results::Results::EclAjdGData adgG;

  if (settings_->simulator()->use_post_sim_script()) {
    string expected_scr_path = paths_.GetPath(Paths::SIM_WORK_DIR) + "/FO_POSTSIM.sh";
    QStringList* post_sim_args = settings_->simulator()->post_sim_args();

    if (VERB_SIM >= 2) {
      Printer::ext_info("Looking for PostSimWork script at " + expected_scr_path,
                        "Simulation", "Simulator");
    }

    if (Utilities::FileHandling::FileExists(expected_scr_path)) {
      if (VERB_SIM >= 2) {
        Printer::ext_info("PostSimWork script found. Executing... ",
                          "Simulation", "Simulator");
      }

      Utilities::Unix::ExecShellScript(QString::fromStdString(expected_scr_path),
                                       *post_sim_args);

      if(post_sim_args[0].contains("AdjG")) {
        adgG = ReadEclAdjG();
      }


      if (VERB_SIM >= 2) {
        Printer::ext_info("Done executing PostSimWork script.",
                          "Simulation", "Simulator");
      }

    } else {
      Printer::ext_warn("PostSimWork script not found.");
    }
  }

  if (settings_->simulator()->read_external_json_results()) {
    std::string expected_res_path = paths_.GetPath(Paths::SIM_WORK_DIR) + "/FO_EXT_RESULTS.json";

    if (VERB_SIM >= 2) {
      Printer::ext_info("Reading external JSON results at " + expected_res_path,
                        "Simulation", "Simulator");
    }
    auto json_results = Simulation::Results::JsonResults(expected_res_path);
    results_->SetJsonResults(json_results);
  }
}

void Simulator::updateResultsInModel() {
  model_->SetResult("Time", results_->GetValueVector(Results::Results::Property::Time));
  model_->SetResult("FGPT", results_->GetValueVector(Results::Results::Property::CumulativeGasProduction));
  model_->SetResult("FOPT", results_->GetValueVector(Results::Results::Property::CumulativeOilProduction));
  model_->SetResult("FWPT", results_->GetValueVector(Results::Results::Property::CumulativeWaterProduction));
}

}
