/***********************************************************
Created: 28.09.2015 2015 by einar
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
#include <sstream>
#include <cmath>
#include <string>
#include <boost/lexical_cast.hpp>

#include "model.h"
#include "settings_exceptions.h"
#include "trajectory_importer.h"
#include "Settings/helpers.hpp"

#include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>
#include "Utilities/filehandling.hpp"

namespace Settings {

using Printer::ext_warn;
using Printer::ext_info;
using Printer::info;
using Printer::num2str;

Model::Model(QJsonObject json_model, Paths &paths, VerbParams vp) {
  vp_ = vp;

  // [RESERVOIR]
  if (!paths.IsSet(Paths::ENSEMBLE_FILE)) {
    try {
      QJsonObject json_reservoir = json_model["Reservoir"].toObject();
      readReservoir(json_reservoir, paths);
    }
    catch (std::exception const &ex) {
      throw UnableToParseReservoirModelSectionException(
        "Unable to parse reservoir model section: "
          + std::string(ex.what()));
    }
  }

  // [STARTDATE] -------------------------------------------
  start_date_ = QList<int>();
  if (json_model.contains("StartDate")) {
    QJsonArray json_start_date = json_model["StartDate"].toArray();
    for (int i = 0; i < json_start_date.size(); ++i) {
      start_date_.append(json_start_date.at(i).toInt());
    }
  } else {
    start_date_ = { 1, 1, 1980 };
    // throw UnableToParseModelSectionException("StartDate must be specified.");
  }

  // [TSTEP REFINEMENT] ------------------------------------
  if (json_model.contains("TStepRefinement")) {
    tstep_refinement_ = json_model["TStepRefinement"].toInt();
  } else {
    tstep_refinement_ = 1;
  }

  // [CONTROL TIMES] ---------------------------------------
  if (!json_model.contains("ControlTimes") || !json_model["ControlTimes"].isArray()) {
    string em = "ControlTimes array must be defined with at least one time step.";
    throw UnableToParseModelSectionException(em);
  }

  // [CONTROL TIMES] -> inserting NPV-defined control times
  control_times_ = QList<double>();
  if (json_model.contains("NPVInterval")) {
    if (json_model["NPVInterval"].toString().compare("Yearly") == 0) {
      if (json_model.contains("NPVYears")) {
        for (int i = 0; i < json_model["NPVYears"].toInt(); ++i) {
          control_times_.append(365 * i);
        }
      } else {
        throw UnableToParseModelSectionException(
          "Unable to parse NPVYears");
      }
    } else if (json_model["NPVInterval"].toString().compare("Monthly") == 0) {
      if (json_model.contains("NPVMonths")) {
        for (int i = 0; i < json_model["NPVMonths"].toInt(); ++i) {
          control_times_.append(30 * i);
        }
      } else {
        throw UnableToParseModelSectionException("Unable to parse NPVMonths");
      }
    }
  }

  // [CONTROL TIMES] -> adding remaining control times in ControlTimes array
  for (int i = 0; i < json_model["ControlTimes"].toArray().size(); ++i) {
    if (!control_times_.contains(json_model["ControlTimes"].toArray().at(i).toInt())) {
      control_times_.append(json_model["ControlTimes"].toArray().at(i).toInt());
    }
  }
  qSort(control_times_);

  // [WELLS]
  wells_ = QList<Well>();

  // IMPORT WELLS ------------------------------------------
  // [WELLS] -> import trajectories, segmentation
  if (json_model.contains("Import")) {

    auto json_import = json_model["Import"].toObject();
    if (!paths.IsSet(Paths::SIM_DRIVER_FILE)) {
      string em = "SchedulePath must be specified (relative";
      em += " to DriverPath) to use the Import feature.";
      throw std::runtime_error(em);
    }

    // Trajectory import
    if (json_import.contains("ImportTrajectories")) {
      auto json_import_well_names = json_import["ImportTrajectories"].toArray();
      std::vector<std::string> import_well_names;
      for (auto jwn : json_import_well_names) {
        import_well_names.push_back(jwn.toString().toStdString());
      }
      std::string trajectories_path;
      if (!paths.IsSet(Paths::TRAJ_DIR)) {
        trajectories_path = paths.GetPath(Paths::SIM_DRIVER_DIR) + "/trajectories";
        assert(DirExists(trajectories_path, vp_));
      } else {
        trajectories_path = paths.GetPath(Paths::TRAJ_DIR);
      }
      auto traj_importer = TrajectoryImporter(trajectories_path, import_well_names, vp_);

      // set list in well objects
      for (auto wname : import_well_names) {
        for (int i = 0; i < wells_.size(); ++i) {
          if (QString::compare(wells_[i].name, QString::fromStdString(wname)) == 0) {
            wells_[i].imported_wellblocks_ = traj_importer.GetImportedTrajectory(wname);
            wells_[i].definition_type = WellDefinitionType::WellSpline;
            wells_[i].convert_well_blocks_to_spline = false;
          }
        }
      }
    }

    // Segmentation import
    if (json_model.contains("Wells") && json_model["Wells"].isArray()) {
      for (auto jwell : json_model["Wells"].toArray()) { // Go through list of wells in json file
        if (jwell.toObject().contains("Segmentation")) { // Check if the well has the Segmentation keyword
          for (int j = 0; j < wells_.size(); ++j) { // Loop through parsed wells
            if (QString::compare(wells_[j].name, jwell.toObject()["Name"].toString()) == 0) { // Check if well names match
              wells_[j].use_segmented_model = true;
              parseSegmentation(jwell.toObject()["Segmentation"].toObject(), wells_[j]);
              break;
            }
          }
        }
      }
    }

  // READ WELLS --------------------------------------------
  // [WELLS] -> read wells using readSingleWell
  } else {
    try {
      QJsonArray json_wells = json_model["Wells"].toArray();
      for (int i = 0; i < json_wells.size(); ++i) {
        QJsonObject json_well = json_wells[i].toObject();
        wells_.append(readSingleWell(json_well));
      }
    }
    catch (std::exception const &ex) {
      throw UnableToParseWellsModelSectionException(
        "Unable to parse wells model section: " + std::string(ex.what()));
    }
  }
}

void Model::readReservoir(QJsonObject json_reservoir, Paths &paths) {
  // Reservoir grid path
  if (!paths.IsSet(Paths::GRID_FILE) && json_reservoir.contains("Path")) {
    paths.SetPath(Paths::GRID_FILE, json_reservoir["Path"].toString().toStdString());
  }
}

// WELL ----------------------------------------------------
#include "model/model__01_read_well.cpp"

bool Model::controlTimeIsDeclared(int time) const {
  return control_times_.contains(time);
}

std::string Model::Well::ControlEntry::toString() {
  std::stringstream ce;
  ce << "ControlEntry - Timestep:  " << time_step << "\n";
  ce << "               State:     " << (state == WellState::WellOpen ? "Open" : "Shut") << "\n";
  ce << "               Mode:      " << (control_mode == ControlMode::LRATControl ? "Rate" : "BHP") << "\n";
  ce << "               BHP:       " << boost::lexical_cast<std::string>(bhp) << "\n";
  ce << "               Liq. Rate: " << boost::lexical_cast<std::string>(liq_rate) << "\n";
  ce << "               Oil Rate:  " << boost::lexical_cast<std::string>(oil_rate) << "\n";
  ce << "               Gas Rate:  " << boost::lexical_cast<std::string>(gas_rate) << "\n";
  ce << "               Wat. Rate: " << boost::lexical_cast<std::string>(wat_rate) << "\n";
  ce << "               Res. Rate: " << boost::lexical_cast<std::string>(res_rate) << "\n";
  ce << "               Inj. type: " << (injection_type == InjectionType::WaterInjection ? "Water" : "Gas/UNKWN") << "\n";
  ce << "               Variable:  " << (is_variable ? "Yes" : "No") << "\n";
  return ce.str();
}

std::string Model::Well::toString(std::string sp) {
  std::stringstream ce;
  ce << "Well - Name:           " << name.toStdString() << sp;
  ce << "       Type:           " << (type == WellType::Injector ? "Injector" : "Producer") << sp;
  ce << "       Group:          " << group.toStdString() << sp;
  ce << "       Radius:         " << wellbore_radius << sp;
  ce << "       Direction:      " << direction << sp;
  ce << "       Pref. phase:    " << preferred_phase << sp;
  ce << "       Def. type:      " << definition_type << sp;
  ce << "       N. well blocks: " << well_blocks.size() << sp;
  ce << "       N. controls:    " << controls.size() << sp;
  return ce.str();
}

bool Model::Well::ControlEntry::isDifferent(ControlEntry other) {
  if (state != other.state)
    return true;
  if (liq_rate != other.liq_rate)
    return true;
  if (bhp != other.bhp)
    return true;
  if (control_mode != other.control_mode)
    return true;
  return false; // Assume they're equal if none of the above hits.
}

// COMPARTMENTS --------------------------------------------
#include "model/model__02_read_comp.cpp"

}

