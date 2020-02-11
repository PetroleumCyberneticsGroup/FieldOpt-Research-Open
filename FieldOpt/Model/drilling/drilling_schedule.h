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

#ifndef FIELDOPT_DRILLING_SCHEDULE_H
#define FIELDOPT_DRILLING_SCHEDULE_H

#include "settings.h"
#include "model.h"

namespace Model {
namespace Drilling {

class DrillingSchedule {

 public:
  DrillingSchedule(Settings::Model *settings, Properties::VariablePropertyContainer *variables);

 private:
  vector<int> drilling_steps_;
  map<int, double> time_steps_; //!< Indexed by the drilling steps
  //vector<Settings::Model::Well::WellBlock block
  map<int, tuple<double,double,double>> drilling_points_; //!< Indexed by the drilling steps

  enum DrillingOperation : int {StartDrilling=1, Drilling=2, PullingOutOfHole=3};
  enum ModelType: int {TrueModel=1, Surrogate=2};

  struct DrillingSettings {
    DrillingSettings() {}
    map<int, DrillingOperation> drilling_operations;
    map<int, ModelType> model_types;
    map<int, bool> is_variable_drilling_points;
    map<int, bool> is_variable_completions;
  };

};

}
}
#endif //FIELDOPT_DRILLING_SCHEDULE_H
