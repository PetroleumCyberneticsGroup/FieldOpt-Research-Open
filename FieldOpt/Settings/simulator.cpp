/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#include <Utilities/printer.hpp>
#include "simulator.h"
#include "settings_exceptions.h"
#include "Utilities/filehandling.hpp"
#include "Settings/helpers.hpp"

namespace Settings {

using namespace Utilities::FileHandling;
using Printer::num2str;
using Printer::info;
using Printer::ext_warn;

Simulator::Simulator(const QJsonObject& json_simulator,
                     Paths &paths, VerbParams vp) {
  vp_ = vp;
  setStructure(json_simulator);
  setPaths(json_simulator, paths);
  setType(json_simulator);
  setParams(json_simulator);
  setCommands(json_simulator);
  setFluidModel(json_simulator);
}

void Simulator::setStructure(QJsonObject json_simulator) {
  if(vp_.vSET >= 5) { info("Parsing FileStructure.", vp_.lnw); }

  if (json_simulator.contains("FileStructure")) {
    QString type = json_simulator["FileStructure"].toString();
    if (QString::compare(type, "Branched") == 0) {
      file_structure_.type_ = FileStructType::Branched;
      set_opt_prop_int(file_structure_.levels_num_,
        json_simulator, "BranchLevels", vp_);
      for (int ii=0; ii < file_structure_.levels_num_; ii++) {
        file_structure_.levels_str_ = file_structure_.levels_str_ + "../";
      }
    }
  } else {
    file_structure_.type_ = FileStructType::Flat;
    file_structure_.levels_num_ = 0;
    file_structure_.levels_str_ = "";
  }
}

void Simulator::setPaths(QJsonObject json_simulator, Paths &paths) {
  if(vp_.vSET >= 5) { info("Setting paths.", vp_.lnw); }

  sched_file_name_ = json_simulator["ScheduleFile"].toString();

  if (!paths.IsSet(Paths::ENSEMBLE_FILE)) {
    is_ensemble_ = false;

    if(vp_.vSET >= 4) { info("Setting SIM_DRIVER_FILE, SIM_DRIVER_DIR", vp_.lnw); }
    if (!paths.IsSet(Paths::SIM_DRIVER_FILE)
        && json_simulator.contains("DriverPath")) {
      paths.SetPath(Paths::SIM_DRIVER_FILE,
                    json_simulator["DriverPath"].toString().toStdString());
    }

    QString sim_drvr_dir = GetParentDirPath(paths.GetPathQstr(Paths::SIM_DRIVER_FILE));
    paths.SetPath(Paths::SIM_DRIVER_DIR, sim_drvr_dir.toStdString());

    if(vp_.vSET >= 4) { info("Setting CASE_DRVR_DIR, CASE_ROOT_DIR, SIM_SCH_FILE", vp_.lnw); }
    QStringList parts = sim_drvr_dir.split("/");
    string root_drvr_path;
    for (int ii=parts.length() - file_structure_.levels_num_ ; ii < parts.length(); ii++) {
      root_drvr_path += "/" + parts.at(ii).toStdString();
    }
    paths.SetPath(Paths::CASE_DRVR_DIR, root_drvr_path, true);

    string schedule_path;
    if (!paths.IsSet(Paths::SIM_SCH_FILE)
        && json_simulator.contains("ScheduleFile")) {

      if (file_structure_.type_ == FileStructType::Flat) {
        schedule_path = paths.GetPath(Paths::SIM_DRIVER_DIR) + "/"
            + sched_file_name_.toStdString();

      } else if (file_structure_.type_ == FileStructType::Branched) {
        schedule_path = paths.GetPath(Paths::SIM_DRIVER_DIR) + "/"
            + file_structure_.levels_str_ + sched_file_name_.toStdString();
        paths.SetPath(Paths::CASE_ROOT_DIR,
                      GetParentDirPath(GetAbsFilePath(schedule_path)));
      }

      paths.SetPath(Paths::SIM_SCH_FILE, schedule_path, true);
    }

  } else { // Ensemble run
    is_ensemble_ = true;
    ensemble_ = Ensemble(paths.GetPath(Paths::ENSEMBLE_FILE), vp_);
    // Set the data file path to the first realization so that the deck parser can find it
    paths.SetPath(Paths::SIM_DRIVER_FILE,
                  ensemble_.GetRealization(ensemble_.GetAliases()[0]).data());
    if (json_simulator.contains("SelectRealizations")) {
      ensemble_.SetNSelect(json_simulator["SelectRealizations"].toInt());
    }
  }
}

void Simulator::setType(QJsonObject json_simulator) {
  if(vp_.vSET >= 5) { info("Setting Simulator type.", vp_.lnw); }

  QString type = json_simulator["Type"].toString();
  if (QString::compare(type, "ECLIPSE") == 0) {
    type_ = SimulatorType::ECLIPSE;
  } else if (QString::compare(type, "ADGPRS") == 0) {
    type_ = SimulatorType::ADGPRS;
  } else if (QString::compare(type, "Flow", Qt::CaseInsensitive) == 0) {
    type_ = SimulatorType::Flow;
  } else if (QString::compare(type, "IX", Qt::CaseInsensitive) == 0) {
    type_ = SimulatorType::INTERSECT;
  } else {
    string em = "Simulator type " + type.toStdString() + " not recognized.";
    throw SimulatorTypeNotRecognizedException(em);
  }
}

void Simulator::setParams(QJsonObject json_simulator) {
  if(vp_.vSET >= 5) { info("Setting Simulator parameters.", vp_.lnw); }
  set_opt_prop_int(max_minutes_, json_simulator, "MaxMinutes", vp_);
  set_opt_prop_bool(ecl_use_actionx_, json_simulator, "UseACTIONX", vp_);
  set_opt_prop_bool(use_post_sim_script_, json_simulator, "UsePostSimScript", vp_);
  set_opt_prop_bool(use_pre_sim_script_, json_simulator, "UsePreSimScript", vp_);
  set_opt_prop_bool(add_sim_scripts_, json_simulator, "AddSimScripts", vp_);
  set_opt_prop_bool(read_external_json_results_, json_simulator, "ReadExternalJsonResults", vp_);
  set_opt_prop_bool(read_adj_grad_data_, json_simulator, "ReadAdjGrads", vp_);
}

void Simulator::setCommands(QJsonObject json_simulator) {
  if(vp_.vSET >= 5) { info("Setting Simulator commands.", vp_.lnw); }
  script_name_ = "";
  QJsonArray sim_exec_cmds = json_simulator["SimExecCmds"].toArray();
  QJsonArray pre_sim_args = json_simulator["PreSimArgs"].toArray();
  QJsonArray post_sim_args = json_simulator["PostSimArgs"].toArray();
  string tm;

  // ExecScript
  if ((json_simulator.contains("ExecScript")
    && json_simulator["ExecScript"].toString().size() > 0)) {
    script_name_ = json_simulator["ExecScript"].toString();
    tm += "| ExecScript: " + script_name_.toStdString() + " ";
  }

  // ExecScript (backwards compatibility)
  if ((json_simulator.contains("ExecutionScript")
    && json_simulator["ExecutionScript"].toString().size() > 0)) {
    script_name_ = json_simulator["ExecutionScript"].toString();
    tm += "| ExecScript: " + script_name_.toStdString() + " ";
  }

  // SimExecCmds
  if (json_simulator.contains("SimExecCmds")
      && !sim_exec_cmds.isEmpty()) {
    sim_exec_cmds_ = new QStringList();
    for (auto && sim_exec_cmd : sim_exec_cmds) {
      sim_exec_cmds_->append(sim_exec_cmd.toString());
    }
    tm += "| SimExecCmds: " + sim_exec_cmds_->join(" ").toStdString() + " ";
  }

  // PreSimArgs
  if (json_simulator.contains("PreSimArgs")
      && !pre_sim_args.isEmpty()) {
    pre_sim_args_ = new QStringList();
    for (auto && pre_sim_arg : pre_sim_args) {
      pre_sim_args_->append("\"" + pre_sim_arg.toString() + "\"");
    }
    tm += "| PreSimArgs: " + pre_sim_args_->join(" ").toStdString() + " ";
  }

  // PostSimArgs
  if (json_simulator.contains("PostSimArgs")
      && !post_sim_args.isEmpty()) {
    post_sim_args_ = new QStringList();
    for (auto && post_sim_arg : post_sim_args) {

      if(post_sim_arg.toString().isEmpty()) {
        post_sim_args_->append(QString::fromStdString(num2str(post_sim_arg.toInt(),0)));
      } else {
        post_sim_args_->append( "\""+ post_sim_arg.toString() + "\"");
      }
    }
    tm += "| PostSimArgs: " + post_sim_args_->join(" ").toStdString() + " ";
  }

  if (script_name_.length() == 0 && sim_exec_cmds.isEmpty() && vp_.vSET >= 3) {
    string wm = "No execution commands or scripts given in driver file. ";
    wm += "Relying on script path being passed as runtime argument.";
    ext_warn(wm,md_, cl_, vp_.lnw);
  }

  if(vp_.vSET >= 3 && !tm.empty()) {
    ext_info(tm, md_, cl_, vp_.lnw);
  }
}

void Simulator::setFluidModel(QJsonObject json_simulator) {
  if (json_simulator.contains("FluidModel")) {
    QString fluid_model = json_simulator["FluidModel"].toString();
    if (QString::compare(fluid_model, "DeadOil") == 0) {
      fluid_model_ = SimulatorFluidModel::DeadOil;
    } else if (QString::compare(fluid_model, "BlackOil") == 0) {
      fluid_model_ = SimulatorFluidModel::BlackOil;
    }
  } else {
    fluid_model_ = SimulatorFluidModel::BlackOil;
  }
}

}
