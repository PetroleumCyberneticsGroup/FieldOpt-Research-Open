/***********************************************************
Created by Mathias C. Bellout on 27.11.18
Copyright (C) 2018-2019 Mathias Bellout
<chakibbb-pcg@gmail.com>

Modified 2018-2019 Thiago Lima Silva
<thiagolims@gmail.com>

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

#ifndef FIELDOPT_TEST_TR_SUPPORT_HPP
#define FIELDOPT_TEST_TR_SUPPORT_HPP

using std::vector;

#include "Optimization/optimizers/ensemble_exp_value.h"

namespace TestResources {

void CheckSamePolynomials(
    vector<Polynomial> xp, Eigen::MatrixXd pcm,
    double tol, std::string msg) {

  for (int ii=0; ii<xp.size(); ii++) {

    ASSERT_TRUE(xp[ii].coefficients.isApprox(pcm.col(ii),tol));

    if( ~msg.empty() ) {
      stringstream ss; ss << "[          ] " << FMAGENTA;
      cout << ss.str() << msg.c_str() << ii << " -> ok" << END << endl;
    }

  }

}

void CheckSameVector(VectorXd xa, VectorXd xb,
                     double tol, std::string msg, bool dbg=false) {

  stringstream ss; ss << "[          ] " << FMAGENTA;

  if(dbg) {
    cout << ss.str() << msg.c_str() << END << endl;
    cout << ss.str() << "xa: [" << xa.transpose() << "]" << END << endl;
    cout << ss.str() << "xb: [" << xb.transpose() << "]" << END << endl;
  }

  ASSERT_TRUE(xa.isApprox(xb,tol));

  if(dbg && ~msg.empty() ) {
    cout << ss.str() << msg.c_str() << " -> ok" << END << endl;
  }
}

void CheckSameX(VectorXd xa, VectorXd xb,
                vector<int> idx, double tol,
                std::string msg) {

  VectorXd xbr(xb.size());
  for (int ii=0; ii<xb.rows(); ii++) {
    xbr.row(ii) << xb(idx[ii]);
  }
  // cout << "xbr: [" << xbr.transpose() << "]" << endl; // FO tr model data (reordered)
  // cout << "xa: [" << xa.transpose() << "]" << endl; // Matlab tr model data
  // cout << "xb: [" << xb.transpose() << "]" << endl; // FO tr model data

  ASSERT_TRUE(xa.isApprox(xbr,tol));

  if( ~msg.empty() ) {
    stringstream ss; ss << "[          ] " << FMAGENTA;
    cout << ss.str() << msg.c_str() << " -> ok" << END << endl;
  }
}

void FindVarSequence(TestResources::TrustRegionModelData::prob &prob,
                     Optimization::Case &base_case_tr_dfo_probs) {

  // First vector (Eigen + std formats)
  // va -> VectorXd
  VectorXd va = prob.xm.col(0);
  vector<double> v1;
  // v1 -> std
  v1.resize(va.size());
  VectorXd::Map(&v1[0], va.size()) = va;

  // Second vector (Eigen + std formats)
  Eigen::VectorXd vb = base_case_tr_dfo_probs.GetRealVarVector();
  vector<double> v2(vb.size());
  VectorXd::Map(&v2[0], vb.size()) = vb; // <- v2: std

  // --------------
  // Eigen approach
  double tol = 1e-6;
  VectorXi idxa(vb.rows(),1);
  // cout << "va: " << va.transpose() << endl;
  // cout << "vb: " << vb.transpose() << endl << endl;

  for (int ii=0; ii<vb.size(); ii++) {

    // cout << "ii: " << ii << endl << endl;

    for (int jj=0; jj<va.size(); jj++) {

      // cout << "jj: " << jj << endl;
      // cout << "va(jj=" << jj << ") : " << va(jj) << endl;
      // cout << "vb(ii=" << ii << ") : " << vb(ii) << endl;

      if ( norm( va(jj) - vb(ii) ) < tol ) {
        idxa.row(ii) << jj;
        va(jj) = std::numeric_limits<double>::lowest();

        // cout << "----" << endl << "va: " << va.transpose() << endl;
        // cout << "idxa: [ " << idxa.transpose() << " ]" << endl << endl;
        break;

      }

    }
  }
  // dbg
  // cout << "idxa-end: [ " << idxa.transpose() << " ]" << endl;

  // -----------
  // std approah
  vector<int> idx1;

  for (int ii=0; ii<v2.size(); ii++) {

    auto it =
        find_if(v1.begin(), v1.end(), [&](const double& x) {
          return norm( x - v2[ii] ) < tol;
        } );

    idx1.push_back( (int)distance(v1.begin(), it) );
    v1[idx1.back()] = std::numeric_limits<double>::lowest();

  }
  prob.idx = idx1; // <- set to function param

  // dbg
  // cout << "idx1: [ ";
  // for (int ii=0; ii<idx1.size(); ii++) { cout << idx1[ii] << " "; }
  // cout << "]" << endl;

  // --------------------
  // Check b/e approaches
  Eigen::VectorXi idx1e = Eigen::Map<VectorXi>(idx1.data(),
                                               idx1.size());
  if (!idxa.isApprox(idx1e)) { throw runtime_error(""); };

} // End of void FindVarSequence

void PrintCaseData(Optimization::Case &c) {

  // PRINT CASE DATA (ID, X, F)
  stringstream ss; ss << "[          ] " << FMAGENTA;
  cout << ss.str()
       << "---------------------------------------------------------" << END << endl;
  cout << ss.str() << "Case id: " << c.GetId().toString().toStdString() << END << endl;
  cout << ss.str() << "x: [" << c.GetRealVarVector().transpose() << "]" << END << endl;
  cout << ss.str() << "f: [" << c.objf_value() << "]" << END << endl;

}

void OverrideSecondPoint(TestResources::TrustRegionModelData::prob &prob,
                         Optimization::Case &c) {

  VectorXd xbr(prob.xm.rows(),1);
  for (int ii=0; ii<prob.xm.rows(); ii++) {
    xbr.row(ii) << prob.xm(ii,1); // a proper order is ensured previously.
  }

  // dbg
  // cout << "prob.xm: " << endl << prob.xm << endl;
  // cout << "xbr: " << endl << xbr << endl;
  c.SetRealVarValues(xbr.col(0));
  c.set_objf_value(prob.fm(0, 1));
}

void OverrideSecondPointEn(TestResources::TrustRegionModelData::prob &prob,
                           Optimization::Case &c,
                           Optimization::Optimizers::EnsembleExpValue *tr_en_) {

  VectorXd xbr(prob.xm.rows(),1);
  for (int ii=0; ii<prob.xm.rows(); ii++) {
    xbr.row(ii) << prob.xm(ii,1); // a proper order is ensured previously.
  }

  // dbg
  // cout << "prob.xm: " << endl << prob.xm << endl;
  // cout << "xbr: " << endl << xbr << endl;
  c.SetRealVarValues(xbr.col(0));
  c.set_objf_value(tr_en_->initialValue(xbr.col(0)));
}

}
#endif //FIELDOPT_TEST_TR_SUPPORT_HPP
