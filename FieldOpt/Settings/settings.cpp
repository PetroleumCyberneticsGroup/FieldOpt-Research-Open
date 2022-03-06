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
#include <qdebug.h>

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

  if(vp_.vSET >= 5)
    ext_info("Reading driver file.", md_, cl_, vp_.lnw);
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

  if(paths_.IsSet(Paths::RESTART_FILE))
    readRestartFile();

  file->close();
  if(vp_.vSET >= 5)
    ext_info("Finished reading JSON.", md_, cl_, vp_.lnw, 2, 2);
}

void Settings::readRestartFile() {
  if(vp_.vSET >= 5)
    ext_info("Reading Restart file.", md_, cl_, vp_.lnw);

  QFile *file = new QFile(paths_.GetPathQstr(Paths::RESTART_FILE));
  if (!file->open(QIODevice::ReadOnly))
    ext_info("Restart file no found.", md_, cl_, vp_.lnw);

  QJsonDocument json = QJsonDocument::fromJson(file->readAll());
  optimizer_->restart_json_ = new QJsonObject(json.object());
  optimizer_->restart_ = true;

  // dbg
  // qDebug() << json.toJson();
  // QTextStream ts(stdout);
  // ts << "json.toJson(): " << endl;
  // ts << json.toJson();

  // keep for ref.:
  // cout << "rjson0:" << endl; // QJsonDocument::Indented
  // QString strJson0(QJsonDocument(rjson0).toJson(QJsonDocument::Compact));
  // cout << strJson0.toStdString() << endl;

  // QJsonArray array = seto_->restartJson()->value("BestPoint").toArray();
  // foreach (const QJsonValue & v, array)
  //     cout << v.toObject().keys()[0].toStdString() << endl;
}

void Settings::readGlobalSection() {
  if(vp_.vSET >= 5)
    ext_info("Reading Global section.", md_, cl_, vp_.lnw);

  try {
    QJsonObject json_global = json_driver_->value("Global").toObject();
    global_ = new Global(json_global);
  } catch (std::exception const &ex) {
    em_ = "Unable to parse Global section: " + string(ex.what());
    throw UnableToParseGlobalSectionException(em_);
  }
}

void Settings::readSimulatorSection() {
  if(vp_.vSET >= 5)
    ext_info("Reading Simulator section.", md_, cl_, vp_.lnw);

  try {
    QJsonObject json_simulator = json_driver_->value("Simulator").toObject();
    simulator_ = new Simulator(json_simulator, paths_,global()->verbParams());
  } catch (std::exception const &ex) {
    em_ = "Unable to parse Simulator section: " + string(ex.what());
    throw UnableToParseSimulatorSectionException(em_);
  }
}

void Settings::readOptimizerSection() {
  if(vp_.vSET >= 5)
    ext_info("Reading Optimizer section.", md_, cl_, vp_.lnw);

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

