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

#ifndef FIELDOPT_DRILLING_H
#define FIELDOPT_DRILLING_H

#include <QString>
#include <properties/variable_property_container.h>
#include "drilling_schedule.h"
#include "Settings/settings.h"
#include "Settings/model.h"

namespace Model {
namespace Drilling {

class Drilling {

 public:
  Drilling(Settings::Model *settings, Properties::VariablePropertyContainer *variables);

  QString getWellName(){ return well_name_; }

  Properties::VariablePropertyContainer* getVariables() { return variables_;}

 private:
  Properties::VariablePropertyContainer *variables_;
  Settings::Model  *settings_;

  DrillingSchedule *drilling_schedule_;
  QString well_name_;
};

}
}
#endif //FIELDOPT_DRILLING_H
