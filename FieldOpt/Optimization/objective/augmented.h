/***********************************************************
Created by bellout on 8/19/20.

Copyright (C) 2020
Mathias Bellout <chakibbb.pcg@gmail.com>

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

#ifndef FIELDOPT_OPTIMIZATION_OBJECTIVE_AUGMENTED_H_
#define FIELDOPT_OPTIMIZATION_OBJECTIVE_AUGMENTED_H_

#include "objective.h"
#include "Settings/model.h"
#include "Simulation/results/results.h"
#include "Utilities/verbosity.h"

#include <unsupported/Eigen/Splines>

namespace Optimization {
namespace Objective {

using Printer::ext_warn;
using Printer::ext_info;
using Printer::info;
using Printer::ext_error;
using Printer::num2str;
using Printer::pad_text;

using Printer::DBG_prntDbl;
using Printer::DBG_prntVecXd;
using Printer::DBG_prntMatXd;
using Printer::DBG_prntToFile;

using std::string;
using std::vector;
using std::accumulate;
using namespace Simulation::Results;

class Augmented : public Objective {
 public:
/*!
 * \brief Augmented
 * \param
 */

  Augmented(Settings::Optimizer *settings,
            Simulation::Results::Results *results,
            Model::Model *model);

  double value() const override { return -1.0; };
  double value(bool base_case) override;

  void setDbgFileName(string fl) { fl_ = fl; }

  void dbg(int dm=0,
           const VectorXd& v0 = VectorXd::Zero(0),
           const VectorXd& v1 = VectorXd::Zero(0),
           const VectorXd& v2 = VectorXd::Zero(0),
           const VectorXd& v3 = VectorXd::Zero(0),
           const double& d0 = 0.0) const;

 private:
  void setUpAugTerms();
  void setUpWaterCutLimit();

  struct Term {
    string prop_name_str;
    Results::Property prop_name;
    Results::Property prop_type;
    Results::Property prop_spec;
    double coeff;
    bool active;
    string scaling;
    vector<string> wells;
    map<string, vector<int>> segments;
  };

  Simulation::Results::Results *rslts_;
  Settings::Optimizer *seto_;
  Model::Model *model_;

  vector<Augmented::Term*> terms_;
  vector<double> npv_coeffs_ = {0.0, 0.0, 0.0};
  double wcut_limit_ = 1.0;
  double wbreakthrough_ = 0.1;
  vector<double>  objf_scal_;

  string printTerm(Augmented::Term *term, double t_val, int lvl) {
    string pn, pt, ps, pc, im;
    int pad = 20;

    pn = term->prop_name_str;
    pt = rslts_->GetPropertyKey(term->prop_type);
    ps = rslts_->GetPropertyKey(term->prop_spec);
    pc = num2str(term->coeff, 5, 1, 12);

    pad_text(pn, pad);
    pad_text(pt, pad);
    pad_text(ps, pad);
    pad_text(pc, pad);

    im = "prop_name: "+ pn;
    if(lvl == 0) {
      im += "prop_coeff: "+ pc;
      im += "prop_spec: "+ ps + "|";

    } else if(lvl == 1) {
      im += "prop_active: " + to_string(term->active) + "    ";
      im += "prop_type: "+ pt;
      im += "prop_spec: "+ ps;
      im += "term value: " + num2str(t_val, 5, 1, 12) + "|";
    }
    return im;
  }

  string im_ = "", wm_ = "", em_ = "";
  string cl_ = "Augmented";
  string md_ = "Optimization::Objective";

  VectorXd vd_ = VectorXd::Zero(0);
  string fl_ = "dbg_aug_obj.py";
};

}
}

#endif //FIELDOPT_OPTIMIZATION_OBJECTIVE_AUGMENTED_H_
