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

class DrillingSchedule {

 public:
  DrillingSchedule(::Settings::Model::Well well,
  ::Model::Properties::VariablePropertyContainer *variables);

 private:
  Properties::ContinousProperty *time_step_;
  ::Settings::Model::DrillingOperation operation_;


};

#endif //FIELDOPT_DRILLING_SCHEDULE_H
