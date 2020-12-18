/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb-pcg@gmail.com>

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

#include <Utilities/printer.hpp>
#include <Utilities/verbosity.h>
#include "icd.h"

namespace Model {
namespace Wells {
namespace Wellbore {
namespace Completions {

ICD::ICD(const Settings::Model::Well::Completion &completion_settings,
         Properties::VarPropContainer *variable_container)
         : SegmentedCompletion(completion_settings, variable_container) {

  flow_coefficient_ = completion_settings.valve_flow_coeff;
  valve_size_ = new Properties::ContinuousProperty(completion_settings.valve_size);
  valve_size_->setName(completion_settings.name);

  min_valve_size_ = completion_settings.min_valve_size;
  max_valve_size_ = completion_settings.max_valve_size;

  assert(min_valve_size_ < max_valve_size_);
  assert(min_valve_size_ >= 0.0);

  time_step_ = completion_settings.time_step;
  segment_idx_ = completion_settings.icd_segment;
  device_name_ = completion_settings.icd_name;

  if (completion_settings.variable_strength
  || completion_settings.is_variable) {
    variable_container->AddVariable(valve_size_);
  }
}

ICD::ICD(const Settings::Model::Well::ICVGroup &icv_group_settings,
         Properties::VarPropContainer *variable_container)
         : SegmentedCompletion(icv_group_settings, variable_container) {

  flow_coefficient_ = icv_group_settings.valve_flow_coeff;
  valve_size_ = new Properties::ContinuousProperty(icv_group_settings.valve_size);
  valve_size_->setName(icv_group_settings.name);
  device_name_ = icv_group_settings.icv_group_name;

  min_valve_size_ = icv_group_settings.min_valve_size;
  max_valve_size_ = icv_group_settings.max_valve_size;
  assert(min_valve_size_ < max_valve_size_);
  assert(min_valve_size_ >= 0.0);

  time_step_ = icv_group_settings.time_step;
  segment_idxs_ = icv_group_settings.icd_segments;
  device_names_ = icv_group_settings.icvs;

  if (icv_group_settings.is_variable) {
    variable_container->AddVariable(valve_size_);
  }
}

double ICD::setting() {
  if (valveSize() < min_valve_size_ || valveSize() > max_valve_size_) {
    wm_ = "Valve size " + num2str(valveSize()) + "outside bounds, throwing exception.";
    ext_warn(wm_, md_, cl_);
    throw std::runtime_error("Valve size outside bounds.");
  }

  double setting = (valveSize() - min_valve_size_) / (max_valve_size_ - min_valve_size_);

  if (vp_.vMOD >= 1) {
    std::stringstream ss;
    ss << "Computed setting " << setting << " from valve size " << valveSize()
       << ". [min, max] = [ " << min_valve_size_ << " , " << max_valve_size_ << " ]";
    ext_info(ss.str(), md_, cl_);
  }

  return setting;
}

}
}
}
}
