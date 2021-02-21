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

#include "settings.h"
#include "settings_exceptions.h"
#include "Utilities/filehandling.hpp"

#include <QJsonDocument>
#include <iostream>

namespace Settings {

Settings::Settings(Paths &paths) {
  if(vp_.vSET >= 6) {
    im_ = "JSON -> " + paths.GetPath(Paths::DRIVER_FILE);
    ext_info(im_, md_, cl_, vp_.lnw);
  }
  paths_ = paths;
  readDriverFile();
}

QString Settings::GetLogCsvString() const {

  QStringList header = QStringList();
  QStringList content = QStringList();
  header << "name"
         << "maxevals"
         << "initstep"
         << "minstep";
  content << QString::fromStdString(global_->name_)
          << QString::number(optimizer_->parameters().max_evaluations)
          << QString::number(optimizer_->parameters().initial_step_length)
          << QString::number(optimizer_->parameters().minimum_step_length);

  return QString("%1\n%2").arg(header.join(",")).arg(content.join(","));
}

void Settings::readDriverFile() {

  if(vp_.vSET >= 5) {
    im_ = "Reading driver file.";
    ext_info(im_, md_, cl_, vp_.lnw);
  }
  QFile *file = new QFile(paths_.GetPathQstr(Paths::DRIVER_FILE));

  if (!file->open(QIODevice::ReadOnly)) {
    em_ = "Unable to open the driver file";
    throw DriverFileReadException(em_);
  }
  QByteArray data = file->readAll();

  QJsonDocument json = QJsonDocument::fromJson(data);
  if (json.isNull()) {
    em_ = "Unable to parse the input file to JSON.";
    throw DriverFileJsonParsingException(em_);
  }

  if (!json.isObject()) {
    em_= "Driver file format incorrect. Must be a JSON object.";
    throw DriverFileFormatException(em_);
  }
  json_driver_ = new QJsonObject(json.object());

  readGlobalSection();
  readSimulatorSection();
  readModelSection();
  readOptimizerSection();

  file->close();
  if(vp_.vSET >= 5) {
    im_ = "Finished reading JSON.";
    ext_info(im_, md_, cl_, vp_.lnw, 2, 2);
  }
}

void Settings::readGlobalSection() {
  if(vp_.vSET >= 5) {
    im_ = "Reading Global section.";
    ext_info(im_, md_, cl_, vp_.lnw);
  }
  try {
    QJsonObject json_global = json_driver_->value("Global").toObject();
    global_ = new Global(json_global);
  } catch (std::exception const &ex) {
    em_ = "Unable to parse Global section: " + string(ex.what());
    throw UnableToParseGlobalSectionException(em_);
  }
}

void Settings::readSimulatorSection() {
  if(vp_.vSET >= 5) {
    im_ = "Reading Simulator section.";
    ext_info(im_, md_, cl_, vp_.lnw);
  }
  try {
    QJsonObject json_simulator = json_driver_->value("Simulator").toObject();
    simulator_ = new Simulator(json_simulator, paths_,global()->verbParams());
  } catch (std::exception const &ex) {
    em_ = "Unable to parse Simulator section: " + string(ex.what());
    throw UnableToParseSimulatorSectionException(em_);
  }
}

void Settings::readOptimizerSection() {
  if(vp_.vSET >= 5) {
    im_ = "Reading Optimizer section.";
    ext_info(im_, md_, cl_, vp_.lnw);
  }
  try {
    QJsonObject optimizer = json_driver_->value("Optimizer").toObject();
    optimizer_ = new Optimizer(optimizer,global()->verbParams());
  } catch (std::exception const &ex) {
    em_ = "Unable to parse Optimizer section: " + string(ex.what());
    throw UnableToParseOptimizerSectionException(em_);
  }
}

void Settings::readModelSection() {
  if(vp_.vSET >= 5) {
    im_ = "Reading Model section.";
    ext_info(im_, md_, cl_, vp_.lnw);
  }
  try {
    QJsonObject model = json_driver_->value("Model").toObject();
    model_ = new Model(model, paths_,global()->verbParams());
  } catch (std::exception const &ex) {
    em_ = "Unable to parse Model section: " + string(ex.what());
    throw UnableToParseModelSectionException(em_);
  }
}

}

