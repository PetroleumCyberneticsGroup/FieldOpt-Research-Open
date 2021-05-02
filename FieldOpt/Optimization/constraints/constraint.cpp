/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2021 Mathias Bellout
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

#include <iostream>
#include "constraint.h"

namespace Optimization {
namespace Constraints {

Constraint::Constraint(SO& seto, VPC *vars, SV vp) {
  seto_ = seto;
  vars_ = vars;

  logging_enabled_ = false;
  penalty_weight_ = 0.0;
  vp_ = vp;
}

void Constraint::DBG_SnapCase(int cs, string s0, string s1) {
  string tm;
  if (vp_.vOPT >= 4 && cs == 1) {
    tm = "=> " + name() + " enabled: " + num2str(isEnabled(), 0);
    tm += " -- Check bounds: [" + DBG_prntDbl(min_) + DBG_prntDbl(max_) + "] ";
    tm += "for case: " + s0;
    pad_text(tm, vp_.lnw);
    tm += "x: " + s1;
    ext_info(tm, md_, cl_, vp_.lnw);

  } else if (vp_.vOPT >= 4 && cs == 2) {
    pad_text(s0, vp_.lnw );
    tm += "x: " + s1;
    ext_info(tm, md_, cl_, vp_.lnw);
  }
}

void Constraint::EnableLogging(QString output_directory_path) {
  logging_enabled_ = true;
  constraint_log_path_ = output_directory_path + "/log_constraints.txt";
}

Eigen::VectorXd Constraint::GetLowerBounds(QList<QUuid> id_vector) const {
  throw std::runtime_error("Attempted to get bounds from a non-bound constraint.");
}

Eigen::VectorXd Constraint::GetUpperBounds(QList<QUuid> id_vector) const {
  throw std::runtime_error("Attempted to get bounds from a non-bound constraint.");
}

double Constraint::Penalty(Case *c) { return 0.0; }

long double Constraint::PenaltyNormalized(Case *c) {
  return normalizer_.normalize(Penalty(c));
}

void Constraint::InitializeNormalizer(QList<Case *> cases) {
  if (!normalizer_.is_ready()) {
    cout << "WARNING: using default normalization parameter values" << endl;
    normalizer_.set_midpoint(0.0L);
    normalizer_.set_max(1.0L);
    normalizer_.set_steepness(1.0L);
  }
}

void Constraint::PrntWellInfo(string wi, int cs) {

  if (cs == 1) {
    string ws;
    if (seto_.well.isEmpty()) {
      for (const auto& wn : seto_.wells) { ws += wn.toStdString() + "; "; }
    } else { ws = seto_.well.toStdString(); }

    if (vp_.vOPT >= 3) {
      im_ = "Adding " + wi + " constraint for " + seto_.well.toStdString();
      ext_info(im_, md_, cl_, vp_.lnw);

      im_ = "Adding " + wi + " constraint with [min, max] = [";
      im_ += num2str(min_, 5) + ", " + num2str(max_, 5);
      im_ += "] for well " + ws + " with variables: ";
    }
  }

}

}
}

