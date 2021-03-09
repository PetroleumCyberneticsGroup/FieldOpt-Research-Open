/***********************************************************
Created by bellout on 21.02.21
Copyright (C) 2021
Mathias Bellout <chakibbb.pcg@gmail.com>

--
TRDebug.h built from TrustRegionModel.h
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

#ifndef FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRDEBUG_H_
#define FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRDEBUG_H_

#include "TRFrame.h"

struct poly;
struct modMatrix;

namespace Optimization {
namespace Optimizers {

class TRFrame;

class TRDebug {

 public:
  TRDebug(const string& pn, TRFrame* trm);

  string prntDbl(double out, string fd="% 10.3e ", string fn="");

  string prntVecXd(VectorXd vec, string mv="",
                   string fv="% 10.3e ",
                   string fn="");

  string prntMatXd(MatrixXd mat, string mm="",
                   string fm="% 10.3e ");

  void prntToFl(string fn, string sout);


  void prntPolys(string msg, poly p);
  void prntPivotVals(string msg="");
  void prntPivotPolys(string msg="");

  void prntHeader(stringstream &ss, string msg="", int htype=0);
  void prntModelData(string msg);

  string prntSettingsData(string msg);

  void prntFunctionData(
    string fnm = "none", string msg = "",
    VectorXd v0 = VectorXd::Zero(0),
    VectorXd v1 = VectorXd::Zero(0),
    VectorXd v2 = VectorXd::Zero(0),
    double d0 = -9.99, double d1 = -9.99,
    double d2 = -9.99);

  string fn_pivp_, fn_mdat_, fn_sdat_;
  string fn_xchp_, fn_co2m_, fn_pcfs_;

  TRFrame* trm_;
};

}
}
#endif //FIELDOPT_OPTIMIZATION_OPTIMIZERS_DFTR_TRDEBUG_H_