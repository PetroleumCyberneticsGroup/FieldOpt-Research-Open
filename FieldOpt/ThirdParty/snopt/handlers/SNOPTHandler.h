/***********************************************************
 Created by bellout on 2/10/18 from SNOPT template

 Copyright (C) 2017-2019 Mathias Bellout
 <mathias.bellout@petroleumcyberneticsgroup.no>
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

#ifndef FIELDOPT_SNOPTHANDLER_H
#define FIELDOPT_SNOPTHANDLER_H

// ---------------------------------------------------------
#include "../../snopt/include/snopt/snopt.hpp"
// #include "uc_1ThirdParty/snopt/include/snopt/snopt.hpp"
// #include "uc_1ThirdParty/snopt/include/snopt-7.2.12/snopt.hpp"

#include <cstdio>
#include <vector>
#include <iostream>
#include <string>

// ---------------------------------------------------------
using std::cout;
using std::endl;
using std::vector;
using std::string;

namespace Optimization {
namespace Optimizers {

// =========================================================
class SNOPTHandler {

 public:
  SNOPTHandler(const char *prntname,
               const char *summname,
               const char *specname);
  ~SNOPTHandler();

 protected:
  // -------------------------------------------------------
  char Prob[200];
  integer exitCode_;

  integer n, neF;
  integer ObjRow;
  double ObjAdd;

  // -------------------------------------------------------
  integer lenA, neA;
  integer *iAfun, *jAvar;
  double *A;
  integer lenG, neG;
  integer *iGfun, *jGvar;

  // -------------------------------------------------------
  double *x, *xlow, *xupp, *xmul;
  double *F, *Flow, *Fupp, *Fmul;

  integer *xstate, *Fstate;

  // -------------------------------------------------------
  char *xnames, *Fnames;
  integer nxnames, nFnames;

  My_fp usrfun;

  // -------------------------------------------------------
  // SNOPT Files needed
  FILE *specsFile;
  FILE *printFile;

 private:
  // -------------------------------------------------------
  void increment_();
  void decrement_();
  void errorMessageOnExit_(const char *var);
  void init2zero_();
  void areAllUserDataSet_();
  void init_();

  // -------------------------------------------------------
  void setMemory_();
  void alloc_(integer lencw,
              integer leniw,
              integer lenrw);

  void realloc_(integer lencw,
                integer leniw,
                integer lenrw);

  // -------------------------------------------------------
  void memcpyIn_(char *tcw, integer *tiw, double *trw,
                 integer tlencw,
                 integer tleniw,
                 integer tlenrw);

  void memcpyOut(char *tcw, integer *tiw, double *trw,
                 integer tlencw,
                 integer tleniw,
                 integer tlenrw);

  int memoryGuess_(integer &mincw,
                   integer &miniw,
                   integer &minrw);

  // -------------------------------------------------------
  integer inform, fortranStyleObj, fortranStyleAG;
  integer initCalled;
  integer minrw, miniw, mincw;
  integer lenrw, leniw, lencw;

  double *rw;
  integer *iw;
  char *cw;

  // -------------------------------------------------------
  integer iSpecs, iSumm, iPrint;
  char specname[200], printname[200], summname[200];
  integer spec_len, prnt_len, summ_len;

 public:
  // -------------------------------------------------------
  string getErrorMsg();
  int getExitCode();
  void computeJac();

  // -------------------------------------------------------
  void solve(integer starttype,
             vector<double> &xsol,
             vector<double> &fsol);

  void setPrintFile(const char *printname);
  void setSummaryFile(const char *summname);
  void setSpecFile(const char *specname);

  bool has_snopt_option_file;

  // -------------------------------------------------------
  // Functions that set up the problem data:
  void setProblemSize(integer n,
                      integer neF);

  void setObjective(integer ObjRow,
                    double ObjAdd = 0.0);

  void setA(integer lenA,
            integer *iAfun,
            integer *jAvar,
            double *A);

  void setNeA(integer neA);

  void setG(integer lenG,
            integer *iGfun,
            integer *jGvar);

  void setNeG(integer neG);

  // -------------------------------------------------------
  void setX(double *x, double *xlow, double *xupp,
            double *xmul, integer *xstate);

  void setF(double *F, double *Flow, double *Fupp,
            double *Fmul, integer *Fstate);

  // -------------------------------------------------------
  void setXNames(char *xnames, integer nxnames);
  void setFNames(char *Fnames, integer nFnames);
  void setProbName(const char *Prob);
  void setUserFun(My_fp usrfun);

  // -------------------------------------------------------
  void setProblem(integer nvariables,
                  integer nproblemFunctions,
                  integer objectiveRow,
                  integer nA,
                  integer *iAfun,
                  integer *jAvar,
                  double *A,
                  integer nG,
                  integer *iGfun,
                  integer *jGvar,
                  double *x,
                  double *xlow,
                  double *xupp,
                  double *xmul,
                  integer *xstate,
                  double *F,
                  double *Flow,
                  double *Fupp,
                  double *Fmul,
                  integer *Fstate,
                  integer nonzerosOfA,
                  integer nonzerosOfG,
                  My_fp userfun,
                  bool haveDerivatives,
                  integer maxiterations,
                  char *xnames,
                  integer nxnames,
                  char *Fnames,
                  integer nFnames,
                  const char *problemname,
                  const char *printfile,
                  const char *specfile = 0);

  // -------------------------------------------------------
  // Functions for settings/retrieving parameters
  void setParameter(char *stroptin);
  void getParameter(char *stroptin, char *stroptout);

  void setIntParameter(const char *stropt, integer opt);
  void getIntParameter(const char *stropt, integer &opt);
  void setRealParameter(const char *stropt, double opt);
  void getRealParameter(const char *stropt, double &opt);
};

}
}

#endif //FIELDOPT_SNOPTHANDLER_H
