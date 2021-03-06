/***********************************************************
Created: 07.10.2015 2015 by einar

Copyright (C) 2015-2015
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

#include "weightedsum.h"
#include "Model/model.h"

namespace Optimization {
namespace Objective {

using Printer::num2str;

WeightedSum::WeightedSum(Settings::Optimizer *settings,
                         Simulation::Results::Results *results,
                         Model::Model *model) : Objective(settings) {
  results_ = results;
  settings_ = settings;
  components_ = new QList<WeightedSum::Component *>();

  for (int i = 0; i < settings->objective().weighted_sum.size(); ++i) {
    auto *comp = new WeightedSum::Component();
    comp->property = results_->GetPropertyKey(settings->objective().weighted_sum.at(i).property);
    comp->coefficient = settings->objective().weighted_sum.at(i).coefficient;
    comp->time_step = settings->objective().weighted_sum.at(i).time_step;

    if (vp_.vOPT >= 5) {
      string im = "comp.prop: "+ results_->GetPropertyKey(comp->property);
      im += ", comp.time_step: " + num2str(comp->time_step);
      ext_info(im, md_, cl_, vp_.lnw);
    }

    if (settings->objective().weighted_sum.at(i).is_well_prop) {
      comp->is_well_property = true;
      comp->well = settings->objective().weighted_sum.at(i).well;
    } else comp->is_well_property = false;
    components_->append(comp);
  }
  well_economy_ = model->wellCostConstructor();
}

double WeightedSum::value() const {
  double value = 0;
  for (int i = 0; i < components_->size(); ++i) {
    value += components_->at(i)->resolveValue(results_);
  }
  if (well_economy_->use_well_cost) {
    for (auto well: well_economy_->wells_pointer) {
      if (well_economy_->separate) {
        value -= well_economy_->costXY * well_economy_->well_xy[well->name().toStdString()];
        value -= well_economy_->costZ * well_economy_->well_z[well->name().toStdString()];
      } else {
        value -= well_economy_->cost * well_economy_->well_lengths[well->name().toStdString()];
      }
    }
  }
  return value;
}

double WeightedSum::Component::resolveValue(Simulation::Results::Results *results) const {
  if (is_well_property) {
    if (time_step < 0) { // Final time step well property
      return coefficient * results->GetValue(property, well);
    } else { // Non-final time step well property
      return coefficient * results->GetValue(property, well, time_step);
    }
  } else {
    if (time_step < 0) { // Final time step field/misc property
      return coefficient * results->GetValue(property);
    }
    else { // Non-final time step field/misc property
      return coefficient * results->GetValue(property, time_step);
    }
  }
}

}
}
