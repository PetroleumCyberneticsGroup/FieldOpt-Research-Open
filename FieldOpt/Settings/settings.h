/***********************************************************
Created: 28.09.2015 2015 by einar
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2019 Mathias Bellout
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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <QString>
#include <QJsonObject>

#include "Utilities/verbosity.h"

#include "global.h"
#include "simulator.h"
#include "optimizer.h"
#include "model.h"
#include "Settings/paths.h"

namespace Settings {

class Global;
class Simulator;
class Model;
class Optimizer;
class Drilling;

/*!
 * \brief The Settings class contains both general settings
 * for a FieldOpt run and pointers to objects containing
 * specific settings for the Model, Simulator and Optimizer.
 * Settings takes as input the path to a "driver file" in
 * the JSON format.
 *
 * \todo Remove bool verbose_, and associated functions,
 * since this functionality is now handled by int
 * verbosity_level_ in runtime_settings
 */
class Settings
{
  friend class Optimizer;

 public:
  Settings() {}
  Settings(Paths &paths);

  Global *global() const { return global_; }

  //!< Object containing model specific settings.
  Model *model() const { return model_; }

  //!< Object containing optimizer specific settings.
  Optimizer *optimizer() const { return optimizer_; }

  //!< Object containing simulator specific settings.
  Simulator *simulator() const { return simulator_; }

  Drilling *drilling() const { return drilling_; }

  //!< Get a string containing the CSV header and contents
  //!< for the log.
  QString GetLogCsvString() const;

  Paths &paths() { return paths_; }

  void setOptimizer(Optimizer *opt) { optimizer_ = opt; }

  void setSimulator(Simulator *sim) { simulator_ = sim; }

  void setModel(Model *mod) { model_ = mod; }

 private:
  Paths paths_;
  QJsonObject *json_driver_;

  string im_, wm_, em_;
  string md_ = "Settings/settings";
  string cl_ = "Settings";
  VerbParams vp_;

  Global *global_;
  Model *model_;
  Optimizer *optimizer_;
  Simulator *simulator_;
  Drilling *drilling_;

  void readDriverFile();
  void readGlobalSection();
  void readSimulatorSection();
  void readModelSection();
  void readOptimizerSection();

  void readRestartFile();
};

}

#endif // SETTINGS_H
