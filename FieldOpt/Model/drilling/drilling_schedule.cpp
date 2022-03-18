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

#include "drilling_schedule.h"

namespace Model {
namespace Drilling {
DrillingSchedule::DrillingSchedule(Settings::Drilling* setd,
                                   Properties::VarPropContainer *variables) {

  d_sched_ = setd->drilling_sched();
  drilling_steps_ = d_sched_.drilling_steps;
  time_steps_ = d_sched_.time_steps;

  for (int i: drilling_steps_) {
    model_types_.insert(i, (ModelType) d_sched_.model_types.value(i));
    execution_modes_.insert(i, (Execution) d_sched_.execution_modes.value(i));

    drilling_operations_.insert(i, (DrillingOperation) d_sched_.drilling_operations.value(i));
    optimizer_settings_.insert(i, (Settings::Optimizer*) d_sched_.optimizer_settings.value(i));

    is_variable_completions_.insert(i, d_sched_.is_variable_completions.value(i));
    is_variable_drilling_points_.insert(i, d_sched_.is_variable_drilling_points.value(i));

    is_model_updates_.insert(i, d_sched_.is_model_updates.value(i));
    is_warm_start_.insert(i, d_sched_.is_warm_start.value(i));
  }

  assignDrillingPoints(d_sched_.drilling_points);
  assignOptimizationTriggers(d_sched_.optimization_triggers);
}

void DrillingSchedule::assignDrillingPoints(QMap<int,
                                            QList<Settings::Drilling::DrillingPoint>> drilling_points_settings) {
  for (int i: drilling_steps_) {
    QList<DrillingSchedule::DrillingPoint> points_list;
    for (int j = 0; j < drilling_points_settings.value(i).size(); j++) {
      DrillingSchedule::DrillingPoint point;

      point.x = drilling_points_settings.value(i).at(j).x;
      point.y = drilling_points_settings.value(i).at(j).y;
      point.z = drilling_points_settings.value(i).at(j).z;
      point.is_variable = drilling_points_settings.value(i).at(j).is_variable;

      points_list.append(point);
    }
    drilling_points_.insert(i, points_list);
  }
}

void DrillingSchedule::assignOptimizationTriggers(QMap<int, Settings::Drilling::OptimizationTrigger> opt_triggers) {
  for (int i: drilling_steps_) {
    if (opt_triggers.contains(i)) {
      DrillingSchedule::OptimizationTrigger opt_trigger;

      if (opt_triggers.value(i).max_objective_improvement >=0) {
        opt_trigger.max_objective_improvement = opt_triggers.value(i).max_objective_improvement;
      } else {
        opt_trigger.max_objective_improvement = -1;
      }

      if (opt_triggers.value(i).min_model_deviation >=0) {
        opt_trigger.min_model_deviation = opt_triggers.value(i).min_model_deviation;
      } else {
        opt_trigger.min_model_deviation = -1;
      }


      if (opt_triggers.value(i).max_model_deviation >=0) {
        opt_trigger.max_model_deviation = opt_triggers.value(i).max_model_deviation;
      } else {
        opt_trigger.max_model_deviation = 100;
      }

      optimization_triggers_.insert(i, opt_trigger);
    }
  }
}


}
}







