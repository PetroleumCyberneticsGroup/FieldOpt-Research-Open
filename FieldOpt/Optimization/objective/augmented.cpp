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

using Printer::DBG_prntDbl;
using Printer::DBG_prntVecXd;
using Printer::DBG_prntMatXd;

using ECLProp = Simulation::Results::Results::Property;

Augmented::Augmented(Settings::Optimizer *settings,
                     Simulation::Results::Results *results,
                     Model::Model *model) : Objective(settings) {
  settings_ = settings;
  results_ = results;
  model_ = model;
  setUpAugTerms();
}

// ECLProp::WellSegOilFlowRate
// ECLProp::WellSegWatFlowRate
// ECLProp::WellSegPress
// ECLProp::WellSegPressDrop
// ECLProp::WellSegWaterCut
// ECLProp::WellSegXSecArea

double Augmented::value() const {
  stringstream ss;
  // auto timeXd_ = results_->GetValueVectorXd(results_->Time);
  auto stepXd_ = results_->GetValueVectorXd(results_->Step);
  auto ntstep = stepXd_.size();
  double value = 0.0;

  try {
    for (auto term : terms_) {
      double t_val = 0.0;

      if (term->prop_type == ECLProp::FieldTotal) {
        auto rate_vec = results_->GetValueVectorXd(term->prop_spec);
        assert(stepXd_.size() == rate_vec.size());
        t_val = term->coeff * (stepXd_.cwiseProduct(rate_vec)).sum();

      } else if (term->prop_type == ECLProp::GnGFunction) {
        if (term->prop_spec == ECLProp::SegWCTstd) {
          for (const auto& wn : term->wells) {

            // computing alternative well-liq-flow-total (wlft_s) as sum of segment-liquid-flow-totals (slft)
            // instead of using original wlft since wlft is computed using time-steps while slft's are computed
            // using report-steps (not read-in from UNRST), therefore the sum of slft does not match wlft
            // VectorXd wlft = results_->GetValueVectorXd(ECLProp::WellLiqProdTotal, wn); // not used, keep-for-ref
            vector<VectorXd> slft = results_->GetValVectorSegXd(ECLProp::WellSegLiqFlowTotal, wn);

            // computing wlft_s:
            double wlft_s = 0.0;
            for (int jj=0; jj < term->segments[wn].size(); ++jj) {
              wlft_s += slft[jj](ntstep-1);
            }

            // computing weight segment flow total (w_slft)
            VectorXd w_slft(term->segments[wn].size(), 1);
            for (int jj=0; jj < term->segments[wn].size(); ++jj) {
              w_slft(jj) = slft[jj](ntstep-1) / wlft_s;
            }

            ss << "weight segment flow total (w_slft): " + DBG_prntVecXd(w_slft, "", "% 4.3f");
            ss << "; sum: " + num2str(w_slft.sum(),3,0);
            info(ss.str(), vp_.lnw); ss.str("");

            // segment-water-cut
            vector<VectorXd> swct = results_->GetValVectorSegXd(ECLProp::WellSegWaterCut, wn);

            //
            VectorXd wlfr = results_->GetValueVectorXd(ECLProp::WellLiqProdRate, wn);





            // VectorXd wlft = stepXd_.cwiseProduct(wlfr));
            //
            // vector<VectorXd> sofr = results_->GetValVectorSegXd(ECLProp::WellSegOilFlowRate, wn);
            // vector<VectorXd> swfr = results_->GetValVectorSegXd(ECLProp::WellSegWatFlowRate, wn);

            // vector<VectorXd> slfr, slft, soft, swft;


            //

            //   slfr.push_back(sofr[ii] + swfr[ii]);
            //
            //   soft.push_back(stepXd_.cwiseProduct(sofr[ii]));
            //   swft.push_back(stepXd_.cwiseProduct(swfr[ii]));
            //   slft.push_back(soft[ii] + swft[ii]);
            //
            //   sftw(ii) = lrate.cwiseQuotient(wlfr).norm();
            // }

          }
        }

        //   // compute segment rate ratio



      } else if (term->prop_spec == ECLProp::WellWBTTotal) {



      }

      if (vp_.vOPT >= 3) {
        string im = "prop_name: "+ term->prop_name_str;
        im += ", prop_type: "+ results_->GetPropertyKey(term->prop_type);
        im += ", prop_spec: "+ results_->GetPropertyKey(term->prop_spec);
        im += ", term value: " + num2str(t_val, 5, 1);
        info(im, vp_.lnw);
      }
      value += t_val;
    }

  } catch (std::exception const &ex) {
    error("Failed to compute Augmented function "
            + string(ex.what()) + " Returning 0.0");
  }
  return value;
}

void Augmented::setUpAugTerms() {
  for (int ii = 0; ii < settings_->objective().terms.size(); ++ii) {
    auto *term = new Augmented::Term();
    term->prop_name_str = settings_->objective().terms.at(ii).prop_name;
    term->prop_name = results_->GetPropertyKey(term->prop_name_str);
    term->coeff = settings_->objective().terms.at(ii).coefficient;
    term->wells = settings_->objective().terms.at(ii).wells;
    term->segments = settings_->objective().terms.at(ii).segments;
    terms_.push_back(term);

    if (term->prop_name == ECLProp::FieldOilProdTotal) {
      term->prop_type = ECLProp::FieldTotal;
      term->prop_spec = ECLProp::FieldOilProdRate;

    } else if (term->prop_name == ECLProp::FieldGasProdTotal) {
      term->prop_type = ECLProp::FieldTotal;
      term->prop_spec = ECLProp::FieldGasProdRate;

    } else if (term->prop_name == ECLProp::FieldWatProdTotal) {
      term->prop_type = ECLProp::FieldTotal;
      term->prop_spec = ECLProp::FieldWatProdRate;

    } else if (term->prop_name == ECLProp::FieldWatInjTotal) {
      term->prop_type = ECLProp::FieldTotal;
      term->prop_spec = ECLProp::FieldWatInjRate;

    } else if (term->prop_name == ECLProp::SegWCTstd) {
      term->prop_type = ECLProp::GnGFunction;
      term->prop_spec = ECLProp::SegWCTstd;

    } else if (term->prop_name == ECLProp::WellWBTTotal) {
      term->prop_type = ECLProp::GnGFunction;
      term->prop_spec = ECLProp::WellWBTTotal;

    } else {
      term->prop_type = ECLProp::AuxProp;
      term->prop_spec = ECLProp::AuxProp;
    }

    if (vp_.vOPT >= 4) {
      string im = "Term: " + term->prop_name_str;
      im += ", coeff: " + num2str(term->coeff, 5);
      im += ", prop_type: " + results_->GetPropertyKey(term->prop_type);
      im += ", prop_spec: " + results_->GetPropertyKey(term->prop_spec);
      info(im, vp_.lnw);
    }
  }
}

}
}