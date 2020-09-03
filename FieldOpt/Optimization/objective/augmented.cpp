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
using Printer::DBG_prntToFile;

using ECLProp = Simulation::Results::Results::Property;

Augmented::Augmented(Settings::Optimizer *settings,
                     Simulation::Results::Results *results,
                     Model::Model *model) : Objective(settings) {
  settings_ = settings;
  results_ = results;
  model_ = model;
  setUpAugTerms();
}

double Augmented::value() const {
  stringstream ss;
  auto timeXd = results_->GetValueVectorXd(results_->Time);
  auto stepXd = results_->GetValueVectorXd(results_->Step);
  auto ntstep = stepXd.size();

  auto vstep = VectorXd::LinSpaced(ntstep, 1, ntstep);
  double value = 0.0;

  try {
    for (auto term : terms_) {
      double t_val = 0.0;

      if (term->prop_type == ECLProp::FieldTotal) {
        auto rate_vec = results_->GetValueVectorXd(term->prop_spec);
        assert(stepXd.size() == rate_vec.size());
        t_val = term->coeff * (stepXd.cwiseProduct(rate_vec)).sum();

      } else if (term->prop_type == ECLProp::GnGFunction) {
        if (term->prop_spec == ECLProp::SegWCTStdEnd || term->prop_spec == ECLProp::SegWCTStdWbt) {
          for (const auto& wn : term->wells) {

            int nsegs = term->segments[wn].size();
            VectorXd vsegs = VectorXd::LinSpaced(nsegs, 1, nsegs);

            // computing alternative well-liq-flow-total (wlft_s) as sum of segment-liquid-flow-totals (slft)
            // instead of using original wlft since wlft is computed using time-steps while slft's are computed
            // using report-steps (not read-in from UNRST), therefore the sum of slft does not match wlft
            // VectorXd wlft = results_->GetValueVectorXd(ECLProp::WellLiqProdTotal, wn); // not used, keep-for-ref
            vector<VectorXd> slft = results_->GetValVectorSegXd(ECLProp::WellSegLiqFlowTotal, wn);

            // computing wlft_s:
            double wlft_s = 0.0;
            for (int jj=0; jj < nsegs; ++jj) {
              wlft_s += slft[jj](ntstep-1);
            }

            // computing weight segment flow total (w_slft)
            VectorXd w_slft(nsegs, 1);
            for (int jj=0; jj < nsegs; ++jj) {
              w_slft(jj) = slft[jj](ntstep-1) / wlft_s;
            }
            dbg(0, w_slft);

            // segment/well-water-cut
            vector<VectorXd> swct = results_->GetValVectorSegXd(ECLProp::WellSegWaterCut, wn);
            VectorXd wwct = results_->GetValueVectorXd(ECLProp::WellWaterCut, wn);

            // (relative) standard deviation at each report time
            VectorXd vec(nsegs), sd_vec(ntstep), rsd_vec(ntstep);
            vec.fill(0); sd_vec.fill(0); rsd_vec.fill(0);
            for (int ii=0; ii < ntstep; ++ii) {
              for (int jj=0; jj < nsegs; ++jj) {
                vec(jj) = (w_slft(jj) * swct[jj](ii));
              }
              // sample standard deviation, i.e., x.sum()/(N-1)
              sd_vec(ii) = std::sqrt((vec.array() - vec.mean()).square().sum()/(nsegs-1));
              rsd_vec(ii) = sd_vec(ii) * 100 / vec.mean();
              dbg(1, vec);
            }
            dbg(2, timeXd, sd_vec, rsd_vec);

            // compute splines that approximate sd over time
            // (typically, sd will increase with time)
            int const spl_dim = 1;
            const vector<int> spl_deg = {1, 2};
            Spline<double, spl_dim> spline_l, spline_q;

            spline_l = SplineFitting< Spline<double, 1> >::Interpolate(
              sd_vec.transpose(), spl_deg[0], timeXd);

            spline_q = SplineFitting< Spline<double, 1> >::Interpolate(
              sd_vec.transpose(), spl_deg[1], timeXd);

            VectorXd sd_vec_spline_l(ntstep), sd_vec_spline_q(ntstep);
            for (int jj = 0; jj < timeXd.size(); ++jj) {
              sd_vec_spline_l(jj) = spline_l(timeXd(jj)).value();
              sd_vec_spline_q(jj) = spline_q(timeXd(jj)).value();
            }
            dbg(3, sd_vec_spline_l, sd_vec_spline_q);

            // select time point (@end or @wbt)
            int wbt_idx = 0;
            wbt_idx = (int)wwct.rows() - 1;

            if (term->prop_spec == ECLProp::SegWCTStdWbt) {
              for (int ii=0; ii < wwct.size(); ++ii) {
                if (wwct(ii) > .45) {
                  wbt_idx = ii;
                  break;
                }
              }
            }

            ss << "[ECLProp::SegWCTStdWbt] wbt_idx = " << wbt_idx << "; sd_vec(wbt_idx) = " << sd_vec(wbt_idx) << endl;
            t_val = sd_vec(wbt_idx);
            info(ss.str(), vp_.lnw); ss.str("");

            // dbg(4, VectorXd(t_val));
          }
        }

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

    } else if (term->prop_name == ECLProp::SegWCTStdEnd) {
      term->prop_type = ECLProp::GnGFunction;
      term->prop_spec = ECLProp::SegWCTStdEnd;

    } else if (term->prop_name == ECLProp::SegWCTStdWbt) {
      term->prop_type = ECLProp::GnGFunction;
      term->prop_spec = ECLProp::SegWCTStdWbt;

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

void Augmented::dbg(int dm,
                    const VectorXd& v0,
                    const VectorXd& v1,
                    const VectorXd& v2,
                    const VectorXd& v3) const {
  stringstream ss;

  // dbg(0, w_slft);
  if (dm==0 && vp_.vOPT >= 4) {
    // weights liq seg flow total
    ss << "weights_slft    = " + DBG_prntVecXd(v0, ",", "% 5.3f");
    ss << " # sum = " + num2str(v0.sum(), 3, 0);
    DBG_prntToFile(fl_, ss.str() + "\n", "w");
    info(ss.str(), vp_.lnw); ss.str("");
  }

  // dbg(2, timeXd, sd_vec, rsd_vec);
  if (dm==1 && vp_.vOPT >= 5) {
    ss << "weighted_swct   = " + DBG_prntVecXd(v0, ",", "% 8.6f");
    DBG_prntToFile(fl_, ss.str() + "\n");
    info(ss.str(), vp_.lnw); ss.str("");
  }

  // dbg(2, timeXd, sd_vec, rsd_vec);
  if (dm==2 && vp_.vOPT >= 4) {
    if (vp_.vOPT >= 5) {
      ss << "time            = " + DBG_prntVecXd(v0, ",", "% 5.3f");
      DBG_prntToFile(fl_, ss.str()+ "\n");
      info(ss.str(), vp_.lnw); ss.str("");
    }
    ss << "sd_vec_data     = " + DBG_prntVecXd(v1, ",", "% 5.3f");
    ss << " # sum: " + num2str(v1.sum(),3,0);
    ss << " mean: " + num2str(v1.mean(),3,0);
    if (vp_.vOPT >= 5) DBG_prntToFile(fl_, ss.str() + "\n");
    info(ss.str(), vp_.lnw); ss.str("");

    ss << "rsd_vec_data    = " + DBG_prntVecXd(v2, ",", "% 5.3f");
    info(ss.str(), vp_.lnw); ss.str("");
  }

  // dbg(3, sd_vec_spline_l, sd_vec_spline_q);
  if (dm==3 && vp_.vOPT >= 4) {
    ss << "sd_vec_spline_l = " + DBG_prntVecXd(v0, ",", "% 5.3f");
    DBG_prntToFile(fl_, ss.str() + "\n"); info(ss.str(), vp_.lnw); ss.str("");
    ss << "sd_vec_spline_q = " + DBG_prntVecXd(v1, ",", "% 5.3f");
    DBG_prntToFile(fl_, ss.str() + "\n"); info(ss.str(), vp_.lnw); ss.str("");
  }

  if (dm==4 && vp_.vOPT >= 4) {

  }
}

}
}