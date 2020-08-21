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

using RProp = Simulation::Results::Results::Property;

Augmented::Augmented(Settings::Optimizer *settings,
                     Simulation::Results::Results *results,
                     Model::Model *model) : Objective(settings) {
  settings_ = settings;
  results_ = results;
  setUpAugTerms();
}

double Augmented::value() const {
  // auto timeXd_ = results_->GetValueVectorXd(results_->Time);
  auto stepXd_ = results_->GetValueVectorXd(results_->Step);
  double value = 0.0;

  try {
    for (int ii = 0; ii < terms_.size(); ++ii) {
      Augmented::Term *term = terms_[ii];
      double t_val = 0.0;

      if (term->prop_type == RProp::FieldTotal) {
        auto rate_vec = results_->GetValueVectorXd(term->prop_spec);
        assert(stepXd_.size() == rate_vec.size());
        t_val = term->coeff * (stepXd_.cwiseProduct(rate_vec)).sum();

      } else if (term->prop_type == RProp::GnGFunction) {
        if (term->prop_spec == RProp::WellWatBrkTmTotal) {

        }
      }

      if (vp_.vOPT >= 3) {
        string im = "prop_name: "+ term->prop_name_str;
        im += ", prop_type: "+ Results::GetPropertyKey(term->prop_name);
        im += ", prop_spec: "+ Results::GetPropertyKey(term->prop_spec);
        im += ", term value: " + num2str(t_val, 5, 1);
        info(im, vp_.lnw);
      }
      value += t_val;
    }

  } catch (std::exception const &ex) {
    error("Failed to compute Augmented function "
            + string(ex.what()) + "  Returning 0.0");
  }
  return value;
}

void Augmented::setUpAugTerms() {
  for (int ii = 0; ii < settings_->objective().terms.size(); ++ii) {
    auto *term = new Augmented::Term();
    term->prop_name_str = settings_->objective().terms.at(ii).prop_name;
    term->prop_name = results_->GetPropertyKey(term->prop_name_str);
    term->coeff = settings_->objective().terms.at(ii).coefficient;
    terms_.push_back(term);

    if (term->prop_name == RProp::FieldOilProdTotal) {
      term->prop_type = RProp::FieldTotal;
      term->prop_spec = RProp::FieldOilProdRate;

    } else if (term->prop_name == RProp::FieldGasProdTotal) {
      term->prop_type = RProp::FieldTotal;
      term->prop_spec = RProp::FieldGasProdRate;

    } else if (term->prop_name == RProp::FieldWatProdTotal) {
      term->prop_type = RProp::FieldTotal;
      term->prop_spec = RProp::FieldWatProdRate;

    } else if (term->prop_name == RProp::FieldWatInjTotal) {
      term->prop_type = RProp::FieldTotal;
      term->prop_spec = RProp::FieldWatInjRate;

    } else if (term->prop_name == RProp::WellWatBrkTmTotal) {
      term->prop_type = RProp::GnGFunction;
      term->prop_spec = RProp::WellWatBrkTmTotal;

    } else {
      term->prop_type = RProp::AuxProp;
      term->prop_spec = RProp::AuxProp;
    }

    if (vp_.vOPT >= 4) {
      string im = "Term: " + term->prop_name_str;
      im += ", coeff: " + num2str(term->coeff, 5);
      im += ", prop_type: " + Results::GetPropertyKey(term->prop_type);
      im += ", prop_spec: " + Results::GetPropertyKey(term->prop_spec);
      info(im, vp_.lnw);
    }
  }
}

}
}