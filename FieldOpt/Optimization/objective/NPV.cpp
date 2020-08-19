/***********************************************************
Copyright (C) 2015-2020
Brage S. Kristoffersen <brage_sk@hotmail.com>

Modified 2020 Mathias Bellout
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

#include "NPV.h"

#include <iostream>
#include <iomanip>
#include "weightedsum.h"
#include <stdlib.h>
#include <cmath>
#include "Model/model.h"
#include "Model/wells/well.h"
#include <Utilities/printer.hpp>

namespace Optimization {
namespace Objective {

using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using std::scientific;

using Printer::ext_warn;
using Printer::ext_info;
using Printer::info;
using Printer::error;

using namespace Simulation::Results;

NPV::NPV(Settings::Optimizer *settings,
         Simulation::Results::Results *results,
         Model::Model *model) {
  settings_ = settings;
  vp_ = settings->verbParams();

  results_ = results;
  components_ = new QList<NPV::Component *>();

  for (int i = 0; i < settings->objective().NPV_sum.size(); ++i) {
    auto *comp = new NPV::Component();

    if (settings->objective().NPV_sum[i].property.compare(0, 4, "EXT-") == 0 ) {
      comp->is_json_component = true;
      ext_info("Adding external NPV component.", md_, cl_, vp_.lnw);

      comp->property_name = settings->objective().NPV_sum[i].property.substr(4, std::string::npos);
      comp->interval = settings->objective().NPV_sum[i].interval;

    } else {
      comp->is_json_component = false;
      comp->property_name = settings->objective().NPV_sum.at(i).property;
      comp->property = results_->GetPropertyKey(QString::fromStdString(comp->property_name));
    }

    comp->coefficient = settings->objective().NPV_sum.at(i).coefficient;

    if (settings->objective().NPV_sum.at(i).usediscountfactor) {
      comp->interval = settings->objective().NPV_sum.at(i).interval;
      comp->discount = settings->objective().NPV_sum.at(i).discount;
      comp->usediscountfactor = settings->objective().NPV_sum.at(i).usediscountfactor;

    } else {
      comp->interval = "None";
      comp->discount = 0;
      comp->usediscountfactor = false;
    }
    components_->append(comp);
  }
  well_economy_ = model->wellCostConstructor();
}

double NPV::value() const {
  try {
    double value = 0;
    auto report_times = results_->GetValueVector(results_->Time);
    auto report_steps = results_->GetValueVector(results_->Step);
    // Bug? NPV_times needs to be cleared after each component
    // auto NPV_times = new QList<double>;
    auto NPV_times = new QList<QList<double>>;
    auto discount_factor_list = new QList<double>;

    for (int k = 0; k < components_->size(); ++k) {
      auto curr_comp = components_->at(k);
      auto NPV_times_temp = QList<double>();
      
      if (curr_comp->is_json_component) { continue; }

      if (curr_comp->interval == "Yearly") {
        double discount_factor;
        int j = 1;
        for (int i = 0; i < report_times.size(); i++) {
          if (i < report_times.size() - 1 &&  (report_times[i+1] - report_times[i]) > 365) {
            std::stringstream ss;
            ss << "Skipping assumed pre-simulation time step " << report_times[i]
               << ". Next time step: " << report_times[i+1]
               << ". Ignore if this is time 0 in a restart case.";
            ext_warn(ss.str(), md_, cl_, vp_.lnw);
            continue;
          }
          if (std::fmod(report_times.at(i), 365) == 0) {
            discount_factor = 1 / (1 * (pow(1 + curr_comp->discount, j - 1)));
            discount_factor_list->append(discount_factor);
            j += 1;
            NPV_times_temp.append(i);
          }
        }

      } else if (curr_comp->interval == "Monthly") {
        double discount_factor;
        double discount_rate = curr_comp->discount;
        double monthly_discount = Optimization::Objective::NPV::Component::yearlyToMonthly(discount_rate);
        int j = 1;
        for (int i = 0; i < report_times.size(); i++) {
          if (std::fmod(report_times.at(i), 30) == 0) {
            discount_factor = 1 / (1 * (pow(1 + monthly_discount, j - 1)));
            discount_factor_list->append(discount_factor);
            j += 1;
            NPV_times_temp.append(i);
          }
        }

      } else if (curr_comp->interval == "None" || curr_comp->interval == "ReportSteps") {
        QVector<double> qVec = QVector<double>::fromStdVector(report_steps);
        NPV_times_temp = QList<double>::fromVector(qVec);
      }
      NPV_times->append(NPV_times_temp);
    }

    for (int i = 0; i < components_->size(); ++i) {
      auto curr_comp = components_->at(i);
      if (curr_comp->is_json_component) { continue; }

      if (curr_comp->usediscountfactor) {
        for (int j = 1; j < NPV_times->at(i).size(); ++j) {
          auto va = curr_comp->resolveValueDiscount(results_, NPV_times->at(i).at(j));
          auto vb = curr_comp->resolveValueDiscount(results_, NPV_times->at(i).at(j - 1));
          auto prod_difference = va - vb;
          value += prod_difference * curr_comp->coefficient * discount_factor_list->at(i);
        }
      } else {
        Results::Property rate_prop = Results::Property::AuxProp;;

        if (curr_comp->property == Results::Property::FieldOilProdTotal) {
          rate_prop = Results::Property::FieldOilProdRate;

        } else if (curr_comp->property == Results::Property::FieldGasProdTotal) {
          rate_prop = Results::Property::FieldGasProdRate;

        } else if (curr_comp->property == Results::Property::FieldWatProdTotal) {
          rate_prop = Results::Property::FieldWatProdRate;

        } else if (curr_comp->property == Results::Property::FieldWatInjTotal) {
          rate_prop = Results::Property::FieldWatInjRate;
        }

        if (rate_prop != Results::Property::AuxProp) {
          // NPV using rates over report time steps
          for (int jj = 0; jj < NPV_times->at(i).size(); ++jj) {
            value += curr_comp->coefficient * NPV_times->at(i).at(jj) * results_->GetValue(rate_prop, jj);
          }
        } else {
          // NPV using endpoint values of cumulative vectors
          value += curr_comp->resolveValue(results_);
          std::string prop_name = curr_comp->property_name;
        }
      }
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

    for (int j = 0; j < components_->size(); ++j) {
      if (components_->at(j)->is_json_component) {
        if (components_->at(j)->interval == "Single" || components_->at(j)->interval == "None") {
          double extval = components_->at(j)->coefficient
            * results_->GetJsonResults().GetSingleValue(components_->at(j)->property_name);
          value += extval;

        } else {
          ext_warn("Unable to parse external component.", md_, cl_);
        }
      }
    }
    return value;

  } catch (std::exception const &ex) {
    error("Failed to compute NPV: " + string(ex.what()) + "  Returning 0.0");
    return 0.0;
  }
}

double NPV::Component::resolveValue(Simulation::Results::Results *results) {
  return coefficient * results->GetValue(property);
}

double NPV::Component::resolveValueDiscount(Simulation::Results::Results *results,
                                            double time_step) {
  int time_step_int = (int) time_step;
  return results->GetValue(property, time_step_int);
}

double NPV::Component::yearlyToMonthly(double discount_factor) {
  return pow((1 + discount_factor), 0.083333) - 1;
}

}

}
