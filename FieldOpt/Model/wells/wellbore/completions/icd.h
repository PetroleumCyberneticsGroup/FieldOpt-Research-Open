/******************************************************************************
   Copyright (C) 2015-2018 Einar J.M. Baumann <einar.baumann@gmail.com>

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
******************************************************************************/

#ifndef FIELDOPT_ICD_H
#define FIELDOPT_ICD_H

#include "segmented_completion.h"
namespace Model {
namespace Wells {
namespace Wellbore {
namespace Completions {

class ICD : public SegmentedCompletion {
 public:
  ICD(const Settings::Model::Well::Completion &completion_settings,
      Properties::VariablePropertyContainer *variable_container);

  /*!
   * @brief Get the setting for the valve.
   *
   * This is a multiplier [0.0, 1.0] that is applied to the base strength.
   * It should applied using the INVERSE_SQUARE method.
   * It is calculated from the current- and min- and max valve size.
   */
  double setting() const;
  double valveSize() const { return valve_size_->value(); }
  double flowCoefficient() const { return flow_coefficient_; }
  int timeStep() const { return time_step_; }
  int segmentIdx() const { return segment_idx_; }
  std::string deviceName() const { return device_name_; }

  void setValveSize(const double size) { valve_size_->setValue(size); }

 private:
  Properties::ContinousProperty *valve_size_;
  double min_valve_size_, max_valve_size_;
  double flow_coefficient_;
  int time_step_;
  int segment_idx_;
  std::string device_name_;

};

}
}
}
}

#endif //FIELDOPT_ICD_H
