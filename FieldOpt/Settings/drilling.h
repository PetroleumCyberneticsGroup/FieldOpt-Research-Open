/***********************************************************
Copyright (C) 2020-
Thiago L. Silva <thiagolims@gmail.com>
Mathias Bellout <chakibbb-pcg@gmail.com>

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

#ifndef FIELDOPT_SETTINGS_DRILLING_H_
#define FIELDOPT_SETTINGS_DRILLING_H_

#include "settings.h"

namespace Settings {

class Optimizer;

class Drilling
{
  friend class Settings;

 public:
  Drilling(const QJsonObject& json_drilling, VerbParams vp);

  //!< Get the drilling settings in the model
  void readDrilling(QJsonObject json_drilling_workflow);

  struct DrillingPoint {
    DrillingPoint() {}
    double x, y, z;
    bool is_variable = false;
  };

  struct OptimizationTrigger {
    OptimizationTrigger() = default;
    double min_model_deviation = -1;
    double max_model_deviation = 100;
    double max_objective_improvement = -1;
  };

  struct DrillingSchedule {
    DrillingSchedule() = default;
    QList<int> drilling_steps;    //!< Assume the indexes are from ordered from 0 to N
    QMap<int, double> time_steps; //!< Indexed by the drilling steps
    QMap<int, QList<DrillingPoint>> drilling_points; //!< Indexed by the drilling steps

    enum DrillingOperation : int {
      StartDrilling=1, Drilling=2, PullingOutOfHole=3 };
    enum ModelType: int { TrueModel=1, Surrogate=2 };
    enum Execution: int { Serial=1, Parallel=2 };

    QMap<int, DrillingOperation> drilling_operations;
    QMap<int, ModelType> model_types;
    QMap<int, Execution> execution_modes;
    QMap<int, Optimizer*> optimizer_settings;
    QMap<int, OptimizationTrigger> optimization_triggers;

    QMap<int, bool> is_variable_drilling_points;
    QMap<int, bool> is_variable_completions;
    QMap<int, bool> is_model_updates;
    QMap<int, bool> is_warm_start;
  };

  DrillingSchedule drilling_sched() { return drilling_sched_; };
  QString wellName() { return well_name; };
  Optimizer* localOptzrSet() { return local_optimizer_settings; };
  Optimizer* globalOptzrSet() { return global_optimizer_settings; };
  VerbParams verbParams() { return vp_; };

 private:
  QString well_name = "";
  Optimizer* local_optimizer_settings = nullptr;
  Optimizer* global_optimizer_settings = nullptr;

  DrillingSchedule drilling_sched_;

  string im_, wm_, em_;
  string md_ = "Settings";
  string cl_ = "Drilling";
  VerbParams vp_;
};

}
#endif //FIELDOPT_SETTINGS_DRILLING_H_
