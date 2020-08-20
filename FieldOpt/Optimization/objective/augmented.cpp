/***********************************************************
Created by bellout on 8/19/20.

Copyright (C) 2020
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

#include "augmented.h"
#include <Utilities/printer.hpp>

namespace Optimization {
namespace Objective {

using Printer::ext_warn;
using Printer::ext_info;
using Printer::info;
using Printer::error;
using Printer::num2str;

using namespace Simulation::Results;

Augmented::Augmented(Settings::Optimizer *settings,
                     Simulation::Results::Results *results,
                     Model::Model *model) : Objective(settings) {
  settings_ = settings;
  results_ = results;
  setUpAugTerms();
}

double Augmented::value() const {
  auto timeXd_ = results_->GetValueVectorXd(results_->Time);
  auto stepXd_ = results_->GetValueVectorXd(results_->Step);
  Results::Property rate_prop;
  double value = 0.0;

  try {
    for (int ii = 0; ii < terms_.size(); ++ii) {
      Augmented::Term *term = terms_[ii];
      getRateProp(rate_prop, term);

      if (rate_prop != Results::Property::AuxProp) {
        auto rate_vec = results_->GetValueVectorXd(rate_prop);
        auto t_val = term->coeff * (stepXd_ * rate_vec).sum();
        if (vp_.vOPT >= 3) {
          info("Calculating "+ term->prop_name + ": " + num2str(t_val), vp_.lnw);
        }
        value += t_val;
      }

    }
  } catch (std::exception const &ex) {
    error("Failed to compute Augmented function "
            + string(ex.what()) + "  Returning 0.0");
  }
  return value;
}

void Augmented::getRateProp(Results::Property &prop,
                            Term *term) const {
  if (term->prop == Results::Property::FieldOilProdTotal) {
    prop = Results::Property::FieldOilProdRate;

  } else if (term->prop == Results::Property::FieldGasProdTotal) {
    prop = Results::Property::FieldGasProdRate;

  } else if (term->prop == Results::Property::FieldWatProdTotal) {
    prop = Results::Property::FieldWatProdRate;

  } else if (term->prop == Results::Property::FieldWatInjTotal) {
    prop = Results::Property::FieldWatInjRate;

  } else {
    prop = Results::Property::AuxProp;
  }
  if (vp_.vOPT >= 3) {
    info("Prop: " + Results::GetPropertyKey(prop), vp_.lnw);
  }
}

void Augmented::setUpAugTerms() {
  for (int ii = 0; ii < settings_->objective().terms.size(); ++ii) {
    auto *term = new Augmented::Term();
    term->prop_name = settings_->objective().terms.at(ii).prop_name;
    term->prop = results_->GetPropertyKey(term->prop_name);
    term->coeff = settings_->objective().terms.at(ii).coefficient;
    terms_.push_back(term);
    if (vp_.vOPT >= 3) {
      info("Term: " + term->prop_name + ", coeff: " + num2str(term->coeff, 3), vp_.lnw);
    }
  }
}

}
}