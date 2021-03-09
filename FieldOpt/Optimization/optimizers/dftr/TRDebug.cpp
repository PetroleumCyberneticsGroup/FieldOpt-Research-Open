/***********************************************************
Created by bellout on 21.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
TRDebug.cpp built from TrustRegionModel.cpp
created by thiagols on 27.11.18
Copyright (C) 2018
Thiago Lima Silva <thiagolims@gmail.com>
Caio Giuliani <caiogiuliani@gmail.com>
Mathias Bellout <chakibbb.pcg@gmail.com>
--

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

#include "TRDebug.h"

namespace Optimization {
namespace Optimizers {

TRDebug::TRDebug(const string& pn, TRFrame* trm) {
  trm_ = trm;
  string dn =  "tr-dfo-dbg/";
  fn_pivp_ = dn + "DFTR_PivotPolyns_" + pn + ".txt";
  fn_mdat_ = dn + "DFTR_ModelData_" + pn + ".txt";
  fn_sdat_ = dn + "DFTR_SettingsData_" + pn + ".txt";
  fn_xchp_ = dn + "DFTR_ExchPoint_" + pn + ".txt";
  fn_co2m_ = dn + "DFTR_Coeffs2Mat_" + pn + ".txt";
  fn_pcfs_ = dn + "DFTR_PolyCoeffs_" + pn + ".txt";

  string cmdstr = "mkdir " + dn;
  std::system(cmdstr.c_str());

  cmdstr = "rm " + dn + "*" + pn + ".txt";
  std::system(cmdstr.c_str());

  prntToFl(fn_pivp_, fn_pivp_ + "\n");
  prntToFl(fn_mdat_, fn_mdat_ + "\n");
  prntToFl(fn_sdat_, fn_sdat_ + "\n");
  prntToFl(fn_xchp_, fn_xchp_ + "\n");
  prntToFl(fn_co2m_, fn_co2m_ + "\n");
  prntToFl(fn_pcfs_, fn_pcfs_ + "\n");
}

string TRDebug::prntDbl(double out, string fd, string fn) {
  stringstream ss;
  char buffer [100];
  sprintf(buffer, fd.c_str(), out);
  ss << buffer;
  if (fn != "") { prntToFl(fn, ss.str() + "\n"); }
  return ss.str();
}

string TRDebug::prntVecXd(VectorXd vec, string mv,
                          string fv, string fn) {
  stringstream ss;
  if ( vec.size() > 0 ) {
    if (mv != "") { ss << mv; }
    ss << "[";
    for (int ii = 0; ii < vec.size() - 1; ii++) {
      ss << prntDbl(vec(ii), fv);
    }
    ss << prntDbl(vec(vec.size() - 1)) << "]";
  } else {
    ss << "[ Empty ]";
  }

  if (fn != "") { prntToFl(fn, ss.str() + "\n"); }
  return ss.str();
}

string TRDebug::prntMatXd(MatrixXd mat, string mm, string fm) {
  stringstream ss;
  for (int ii = 0; ii < mat.cols(); ii++) {
    if (mm != "") ss << "\n" << mm << "#" << ii;
    ss << prntVecXd(mat.col(ii), mm, fm);
  }
  return ss.str();
}

void TRDebug::prntToFl(string fn, string sout) {
  FILE * pFile;
  pFile = std::fopen(fn.c_str(), "a");
  std::fprintf(pFile, "%s", sout.c_str());
  fclose (pFile);
}

void TRDebug::prntPolys(string msg, poly p) {
  stringstream ss;
  if (msg != "") { prntHeader(ss, msg, 1); }
  ss << prntVecXd(p.coeffs) << endl;
  prntToFl(fn_pcfs_, ss.str());
}

void TRDebug::prntPivotVals(string msg) {
  stringstream ss;
  if (msg != "") { prntHeader(ss, msg, 0); }
  ss << prntVecXd(trm_->pivot_values_) << endl;
  ss << prntVecXd(trm_->pivot_order_) << endl;
  prntToFl(fn_pivp_, ss.str());
}

void TRDebug::prntPivotPolys(string msg) {
  stringstream ss;
  prntHeader(ss, msg);

  for (int kk = 0; kk < trm_->pivot_polys_.size(); kk++) {
    auto vec = trm_->pivot_polys_[kk].coeffs;
    stringstream si; si << "piv.p#" << kk << " ";
    ss << prntVecXd(vec, si.str()) << endl;
  }

  if (trm_->pivot_values_.size() > 0) {
    ss << prntVecXd(trm_->pivot_values_, "piv_vls ") << endl;
  } else {
    ss << "pivot_values_.size(): " << trm_->pivot_values_.size() << endl;
  }
  prntToFl(fn_pivp_, ss.str());
}

void TRDebug::prntHeader(stringstream &ss, string msg, int htype) {
  string ENDSTRN = string(100, '=');

  if (msg != "") {
    if ( htype == 0 ) {
      ss << ENDSTRN << endl;
      ss << "[" << msg << "]" << endl;

    } else if ( htype == 1 ) {
      Printer::pad_text(msg, 28);
      ss << "[" << msg << "]";
    }
  }
}

void TRDebug::prntModelData(string msg) {
  stringstream ss;
  prntHeader(ss, msg);

  ss << "cached_fvalues_: " << prntVecXd(trm_->cached_fvals_) << "\n";
  ss << "cached_points_: "  << prntMatXd(trm_->cached_pts_, " ") << "\n";
  ss << "fvalues_: "        << prntVecXd(trm_->fvals_) << "\n";

  ss << "all_points_: "     << prntMatXd(trm_->all_pts_, " ") << "\n";
  ss << "points_abs_: "     << prntMatXd(trm_->pts_abs_, " ") << "\n";
  ss << "points_shifted_: " << prntMatXd(trm_->pts_shftd_, " ") << "\n";

  ss << "repl_new_points_shifted_: " << prntMatXd(trm_->repl_new_pts_shftd_, " ") << "\n";
  ss << "repl_new_point_shifted_: " << prntMatXd(trm_->repl_new_pt_shftd_, " ") << "\n";
  ss << "repl_new_pivots_: " << prntMatXd(trm_->repl_new_pivots_, " ") << "\n";

  ss << "repl_new_point_abs_: " << prntMatXd(trm_->repl_new_pt_abs_, " ") << "\n";
  ss << "repl_new_fvalues_: "   << prntVecXd(trm_->repl_new_fvals_) << "\n";

  ss << "lb_: " << prntVecXd(trm_->lb_) << "\n";
  ss << "ub_: " << prntVecXd(trm_->ub_) << "\n";

  prntToFl(fn_mdat_, ss.str());
}

string TRDebug::prntSettingsData(string msg) {
  stringstream ss;
  prntHeader(ss, msg);

  ss << "tr parameters      " << " & "  << "value" << "\\\\ \\toprule \n";
  // ss << "tr\\_prob\\_nm: " << " &  " << trm_->tr_prob_name_           << " \\\\ \n";
  ss << "tr\\_init\\_rad:   " << " & $" << prntDbl(trm_->tr_init_rad_)   << " $ \\\\ \n";
  ss << "tr\\_tol\\_f:      " << " & $" << prntDbl(trm_->tr_tol_f_)      << " $ \\\\ \n";
  ss << "tr\\_eps\\_c:      " << " & $" << prntDbl(trm_->tr_eps_c_)      << " $ \\\\ \n";
  ss << "tr\\_eta\\_0:      " << " & $" << prntDbl(trm_->tr_eta_0_)      << " $ \\\\ \n";
  ss << "tr\\_eta\\_1:      " << " & $" << prntDbl(trm_->tr_eta_1_)      << " $ \\\\ \\midrule \n";
  ss << "tr\\_piv\\_thld:   " << " & $" << prntDbl(trm_->tr_pivot_thld_) << " $ \\\\ \n";
  ss << "tr\\_add\\_thld:   " << " & $" << prntDbl(trm_->tr_add_thld_)   << " $ \\\\ \n";
  ss << "tr\\_exch\\_thld:  " << " & $" << prntDbl(trm_->tr_exch_thld_)  << " $ \\\\ \n";
  ss << "tr\\_rad\\_max:    " << " & $" << prntDbl(trm_->tr_rad_max_)    << " $ \\\\ \n";
  ss << "tr\\_rad\\_fac:    " << " & $" << prntDbl(trm_->tr_rad_fac_)    << " $ \\\\ \\midrule \n";
  ss << "tr\\_tol\\_rad:    " << " & $" << prntDbl(trm_->tr_tol_rad_)    << " $ \\\\ \n";
  ss << "tr\\_gamma\\_inc:  " << " & $" << prntDbl(trm_->tr_gamma_inc_)  << " $ \\\\ \n";
  ss << "tr\\_gamma\\_dec:  " << " & $" << prntDbl(trm_->tr_gamma_dec_)  << " $ \\\\ \n";
  ss << "tr\\_crit\\_mu:    " << " & $" << prntDbl(trm_->tr_crit_mu_)    << " $ \\\\ \n";
  ss << "tr\\_crit\\_omega: " << " & $" << prntDbl(trm_->tr_crit_omega_) << " $ \\\\ \\midrule \n";
  ss << "tr\\_crit\\_beta:  " << " & $" << prntDbl(trm_->tr_crit_beta_)  << " $ \\\\ \n";
  ss << "tr\\_lower\\_b:    " << " & $" << prntDbl(trm_->tr_lower_bnd_)  << " $ \\\\ \n";
  ss << "tr\\_upper\\_b:    " << " & $" << prntDbl(trm_->tr_upper_bnd_)  << " $ \\\\ \n";
  ss << "tr\\_iter\\_max:   " << " & $" << prntDbl(trm_->tr_iter_max_)   << " $ \\\\ \n";
  ss << "tr\\_num\\_initx:  " << " & $" << prntDbl(trm_->tr_num_init_x_) << " $ \\\\ \n";
  ss << "tr\\_basis\\_:     " << " & $" << trm_->tr_basis_               << " $ \\\\ \n";
  ss << "tr\\_init\\_smpln: " << " & $" << trm_->tr_init_smpln_          << " $ \\\\ \n";
  ss << "tr\\_rng\\_seed:   " << " & $" << prntDbl(trm_->tr_rng_seed_)   << " $ \\\\ \\bottomrule \n";

  prntToFl(fn_sdat_, ss.str());
  return ss.str();
}

void TRDebug::prntFunctionData(string fnm, string msg,
                               VectorXd v0, VectorXd v1, VectorXd v2,
                               double d0, double d1, double d2) {

  stringstream ss;
  prntHeader(ss, "FOEx");
  string fn;

  if (fnm == "exchangePoint") {
    ss << "[ p_abs_.rows(): " << prntDbl(trm_->pts_abs_.rows(),   "% 2.0f") << " ]";
    ss << "[ last_p: "        << prntDbl(trm_->pts_abs_.cols()-1, "% 2.0f") << " ]";
    ss << "[ center_i: "      << prntDbl((double)trm_->tr_center_,         "% 2.0f") << " ]";
    ss << "[ radius_: "       << prntDbl((double)trm_->radius_,"% 6.3e") << " ]\n";
    ss << "new_point: "       << prntVecXd(v0);
    ss << "shift_center: "    << prntVecXd(v1) << "\n";
    ss << "pivot_values_: "   << prntVecXd(trm_->pivot_values_) << "\n";

    ss << "cached_fvalues_: " << prntVecXd(trm_->cached_fvals_) << "\n";
    ss << "cached_points_: "  << prntMatXd(trm_->cached_pts_, " ") << "\n";
    ss << "fvalues_: "        << prntVecXd(trm_->fvals_) << "\n";
    ss << "points_abs_: "     << prntMatXd(trm_->pts_abs_, " ") << "\n";
    ss << "points_shifted_: " << prntMatXd(trm_->pts_shftd_, " ") << "\n";
    fn = fn_xchp_;

  } else if (fnm == "none") {
    cout << "Specify function name" << "\n";
  }

  prntToFl(fn, ss.str());
}

}
}