/********************************************************************
 Copyright (C) 2020-
 Mathias Bellout <chakibbb-pcg@gmail.com>
 Thiago L. Silva <thiagolims@gmail.com>

 This file is part of the FieldOpt project.

 FieldOpt is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FieldOpt is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FieldOpt.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#include "drilling_schedule.h"
using namespace std;

namespace Model {
namespace Drilling {
DrillingSchedule::DrillingSchedule(Settings::Model *settings, Properties::VariablePropertyContainer *variables) {

  drilling_steps_ = settings->drilling().drilling_schedule.drilling_steps;
  time_steps_ = settings->drilling().drilling_schedule.time_steps;


  for (int i: drilling_steps_) {
    model_types_.insert(i, (ModelType) settings->drilling().drilling_schedule.model_types.value(i));
    drilling_operations_.insert(i, (DrillingOperation) settings->drilling().drilling_schedule.drilling_operations.value(i));
    is_variable_completions_.insert(i, settings->drilling().drilling_schedule.is_variable_completions.value(i));
    is_variable_drilling_points_.insert(i, settings->drilling().drilling_schedule.is_variable_drilling_points.value(i));
    is_model_updates_.insert(i, settings->drilling().drilling_schedule.is_model_updates.value(i));
  }

  assignDrillingPoints(settings->drilling().drilling_schedule.drilling_points);
}

void DrillingSchedule::assignDrillingPoints(QMap<int, QList<Settings::Model::Drilling::DrillingPoint>> drilling_points_settings) {
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

}
}







