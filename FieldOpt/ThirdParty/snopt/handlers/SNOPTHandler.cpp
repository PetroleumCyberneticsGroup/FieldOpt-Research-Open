/***********************************************************
Created by bellout on 2/10/18 from sqp_snopt template

Copyright (C) 2017-2019 Mathias Bellout
<mathias.bellout@petroleumcyberneticsgroup.no>
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

#include "SNOPTHandler.h"
#include "../../../Utilities/filehandling.hpp"
#include "../../../Utilities/colors.hpp"

using Utilities::FileHandling::FileExists;

namespace Optimization {
namespace Optimizers {

SNOPTHandler::SNOPTHandler(const char *prntname,
                           const char *summname,
                           const char *specname) {
  vp_ = {};
  vp_.vUTI = 1;

  initCalled = 0;
  init2zero_();

  // Nothing incremented yet
  fortranStyleObj = 0;
  fortranStyleAG = 0;

  // Create temporary memory for the call to sninit_.
  // Lengths must all be >= 500.
  lencw = 500;
  leniw = 500;
  lenrw = 500;
  this->alloc_(500, 500, 500);

  // sninit_ "undefines" the optional parameters

  iPrint = 12;
  setPrintFile(prntname);
  iSumm = 29;  // Console/terminal: iSumm = 6;
  // setSummaryFile( summname );
  this->init_();

  this->has_snopt_option_file = true;
  if (!FileExists(QString(specname), vp_, md_, cl_)) {
    this->has_snopt_option_file = false;
    cout << "-- sqp_snopt: Couldn't open file: " << string(specname) << "\n";
    cout << "-- sqp_snopt: Using default options for sqp_snopt\n";
  }

  if (this->has_snopt_option_file) {
    iSpecs = 15;
    setSpecFile(specname);
  } else {
    iSpecs = 0;
  }
}


int SNOPTHandler::getExitCode() {
  return exitCode_;
}


// Destructor
SNOPTHandler::~SNOPTHandler() {

  // Close print and spec files if necessary.
  if (iPrint != 0) {
    closefile_(&iPrint);
  }


  if ((iSumm != 0) && (iSumm != 6)) {
    closefile_(&iSumm);
  }


  if (iSpecs != 0) {
    closefile_(&iSpecs);
  }

  // Delete work arrays.
    delete[] rw;
    delete[] iw;
    delete[] cw;
}


void SNOPTHandler::init2zero_() {
  // Data that must be set by user.

  n = 0;
  neF = 0;

  ObjRow = 0;
  ObjAdd = 0;

  iAfun = 0;
  jAvar = 0;
  A = 0;
  iGfun = 0;
  jGvar = 0;

  x = 0;
  xlow = 0;
  xupp = 0;
  xmul = 0;
  F = 0;
  Flow = 0;
  Fupp = 0;
  Fmul = 0;

  xstate = 0;
  Fstate = 0;
  nxnames = 0;
  nFnames = 0;

  usrfun = 0;

  neA = -1;  // Indicate that neA is yet to be assigned
  neG = -1;
  lenA = -1;
  lenG = -1;
}


void SNOPTHandler::areAllUserDataSet_() {

  if (n == 0) errorMessageOnExit_("n");
  if (neF == 0) errorMessageOnExit_("neF");

  if (x == 0) errorMessageOnExit_("x");
  if (xlow == 0) errorMessageOnExit_("xlow");
  if (xupp == 0) errorMessageOnExit_("xupp");
  if (xmul == 0) errorMessageOnExit_("xmul");

  if (F == 0) errorMessageOnExit_("F");
  if (Flow == 0) errorMessageOnExit_("Flow");
  if (Fupp == 0) errorMessageOnExit_("Fupp");
  if (Fmul == 0) errorMessageOnExit_("Fmul");

  if (xnames == 0) errorMessageOnExit_("xnames");
  if (Fnames == 0) errorMessageOnExit_("Fnames");
  if (nxnames == 0) errorMessageOnExit_("nxnames");
  if (nFnames == 0) errorMessageOnExit_("nFnames");

  if (usrfun == 0) errorMessageOnExit_("usrfun");
  if (lenA == -1) errorMessageOnExit_("lenA");
  if (lenG == -1) errorMessageOnExit_("lenG");

  if (neA > 0 && iAfun == 0) errorMessageOnExit_("iAfun");
  if (neA > 0 && jAvar == 0) errorMessageOnExit_("jAvar");
  if (neA > 0 && A == 0) errorMessageOnExit_("A");

  if (neG > 0 && iGfun == 0) errorMessageOnExit_("iGfun");
  if (neG > 0 && jGvar == 0) errorMessageOnExit_("jGvar");

}


void SNOPTHandler::errorMessageOnExit_(const char *var) {
  throw(string(var) + " must be set prior to call to \n" +
      "SNOPTOptimizer::solve() or SNOPTOptimizer::computeJac()");
}


void SNOPTHandler::setMemory_() {

  int memoryGuess;
  memoryGuess = this->memoryGuess_(mincw, miniw, minrw);


  if (mincw > lencw || miniw > leniw || minrw > lenrw) {

    // Reallocate memory while retaining the values set in sninit_
    this->realloc_(mincw, miniw, minrw);

    // Save the lengths of the new work arrays.
    this->setIntParameter("Total real workspace   ", lenrw);
    this->setIntParameter("Total integer workspace", leniw);

    // -----------------------------------------------------
    // Did we have to guess values of neA and neG for snmema()
    if (memoryGuess == 1) {
      this->computeJac();

      memoryGuess = this->memoryGuess_(mincw, miniw, minrw);
      assert(memoryGuess == 0);

      this->realloc_(mincw, miniw, minrw);
      this->setIntParameter("Total real workspace   ", lenrw);
      this->setIntParameter("Total integer workspace", leniw);
    }
  }
}


void SNOPTHandler::alloc_(integer alencw,
                          integer aleniw,
                          integer alenrw) {

  // Reset work array lengths.
  lencw = alencw;
  leniw = aleniw;
  lenrw = alenrw;

  // Allocate new memory for work arrays.
  cw = new char[8 * lencw];
  iw = new integer[leniw];
  rw = new doublereal[lenrw];
}


void SNOPTHandler::realloc_(integer alencw,
                            integer aleniw,
                            integer alenrw) {


  // Call to this->alloc will overwrite these values => must save.
  integer tlencw = lencw;
  integer tleniw = leniw;
  integer tlenrw = lenrw;


  // Call to this->alloc will create new values for cw, iw, rw => must save.
  char *tcw = cw;
  integer *tiw = iw;
  double *trw = rw;


  // Allocate new memory
  this->alloc_(alencw, aleniw, alenrw);

  // Copy in old values, previously set
  this->memcpyIn_(tcw, tiw, trw, tlencw, tleniw, tlenrw);

  // Delete temporary work arrays
  delete[] trw;
  delete[] tiw;
  delete[] tcw;

}


void SNOPTHandler::memcpyIn_(char *tcw,
                             integer *tiw, double *trw,
                             integer tlencw,
                             integer tleniw,
                             integer tlenrw) {

  integer mlencw = lencw < tlencw ? lencw : tlencw;
  integer mleniw = leniw < tleniw ? leniw : tleniw;
  integer mlenrw = lenrw < tlenrw ? lenrw : tlenrw;

  memcpy(cw, tcw, 8 * mlencw * sizeof(char));
  memcpy(iw, tiw, mleniw * sizeof(integer));
  memcpy(rw, trw, mlenrw * sizeof(double));
}


void SNOPTHandler::memcpyOut(char *tcw,
                             integer *tiw, double *trw,
                             integer tlencw,
                             integer tleniw,
                             integer tlenrw) {

  integer mlencw = lencw < tlencw ? lencw : tlencw;
  integer mleniw = leniw < tleniw ? leniw : tleniw;
  integer mlenrw = lenrw < tlenrw ? lenrw : tlenrw;

  memcpy(tcw, cw, 8 * mlencw * sizeof(char));
  memcpy(tiw, iw, mleniw * sizeof(integer));
  memcpy(trw, rw, mlenrw * sizeof(double));
}


void SNOPTHandler::increment_() {


  if (!fortranStyleObj) {
    // Increment row indicator.
    ObjRow++;
    fortranStyleObj = 1;
  }


  if (!fortranStyleAG) {
    // Increment A indices.
    int k;
    for (k = 0; k < neA; k++) {
      iAfun[k]++;
      jAvar[k]++;
    }

    // -----------------------------------------------------
    // Increment G indices.
    for (k = 0; k < neG; k++) {
      iGfun[k]++;
      jGvar[k]++;
    }
    fortranStyleAG = 1;
  }
}


void SNOPTHandler::decrement_() {


  if (fortranStyleObj) {
    // Decrement row indicator.
    ObjRow--;
    fortranStyleObj = 0;
  }


  if (fortranStyleAG) {
    int k;

    // Decrement A indices.
    for (k = 0; k < neA; k++) {
      iAfun[k]--;
      jAvar[k]--;
    }

    // Decrement G indices.
    for (k = 0; k < neG; k++) {
      iGfun[k]--;
      jGvar[k]--;
    }
    fortranStyleAG = 0;
  }
}


void SNOPTHandler::computeJac() {

  cout << " ---------- ---------- ---------- ---------- " << endl;
  cout << " ---------- ---------- ---------- ---------- " << endl;


  // Ensures all user data has been initialized.
  // areAllUserDataSet_();
  // this->memoryGuess_(mincw, miniw, minrw);

  // if (mincw > lencw || miniw > leniw || minrw > lenrw) {
  //   // Reallocate memory while retaining the values set in sninit_
  //   this->realloc_(mincw, miniw, minrw);
  //   cout << "ait -----------------" << endl;
  //   // Save the lengths of the new work arrays.
  //   this->setIntParameter("Total real workspace   ", lenrw);
  //   this->setIntParameter("Total integer workspace", leniw);
  // }


  // cout << mincw << "\t" << miniw << "\t" << minrw << endl;
  // iAfun = 0;jAvar = 0;lenA=0;neA=0;A=0;
  // cout << iAfun  << "\t" << jAvar << "\t" << lenA << "\t" << neA  << "\t" << A  << endl;
  // cout << leniw << "\t" << lencw << "\t" << lenrw << endl;
  cout << iGfun << "\t" << jGvar << "\t"
       << lenG  << "\t" << neG   << "\t" << endl;


  char *c_cw;
  integer *c_iw;
  double *c_rw;
  c_cw = new char[8 * lencw];
  c_iw = new integer[leniw];
  c_rw = new doublereal[lenrw];


  for (int i = 0; i < lencw; i++){
      c_cw[i]  = cw[i];
  }

  for (int i = 0; i < lenrw; i++){
      c_rw[i]  = rw[i];
  }

  for (int i = 0; i < leniw; i++){
      c_iw[i]  = iw[i];
  }

  for(int i = 0; i < n; i++){
      cout << x[i] << "\t";
  }


  integer c_inform = inform;
  integer c_neF = neF;
  integer c_n = n;
  integer c_lenA = lenA;
  integer c_neA = neA;
  integer c_lenG = lenG;
  integer c_neG = neG;
  integer c_lencw = lencw;
  integer c_leniw = leniw;
  integer c_lenrw = lenrw;
  integer c_mincw = mincw;
  integer c_miniw = miniw;
  integer c_minrw = minrw;


  integer *c_iAfun, *c_jAvar;
  double *c_A;
  integer *c_iGfun, *c_jGvar;

  c_iAfun = new integer[lenA];
  c_jAvar = new integer[lenA];
  c_A     = new double[lenA];

  c_iGfun = new integer[lenG];
  c_jGvar = new integer[lenG];

  double *c_x, *c_xlow, *c_xupp;


  c_x = new double[n];
  c_xlow = new double[n];
  c_xupp = new double[n];

  for (int i = 0; i < n; i++){
    c_x[i] = x[i];
    c_xlow[i] = xlow[i];
    c_xupp[i] = xupp[i];
  }

  for (int i = 0; i < lenA; i++){
    c_A[i] = A[i];
    c_iAfun[i] = iAfun[i];
    c_jAvar[i] = jAvar[i];
  }

  for (int i = 0; i < lenG; i++){
    c_iGfun[i] = iGfun[i];
    c_jGvar[i] = jGvar[i];
  }


  snjac_(&inform, &neF, &n, usrfun,
       iAfun, jAvar, &lenA, &neA, A,
       iGfun, jGvar, &lenG, &neG,
       x, xlow, xupp, &mincw, &miniw, &minrw,
       cw, &lencw, iw, &leniw, rw, &lenrw,
       cw, &lencw, iw, &leniw, rw, &lenrw,
       8 * 500, 8 * 500);

  //cout << iGfun  << "\t" << jGvar << "\t" << lenG << "\t" << neG  << "\t" << endl;
  //cout << leniw << "\t" << lencw << "\t" << lenrw << endl;
  //cout << iAfun  << "\t" << jAvar << "\t" << lenA << "\t" << neA  << "\t" << A  << endl;
  //cout << mincw << "\t" << miniw << "\t" << minrw << endl;


  for(int i = 0; i < n; i++){
      if(c_x[i] != x[i]) {
        cout << "x: " << i << "\t"
             << c_x[i] << "  !=  " << x[i] << endl;
      }

      if(c_xlow[i] != xlow[i]) {
        cout << "xlow: " << i << "\t"
             << c_xlow[i] << "  !=  " << xlow[i] << endl;
      }

      if(c_xupp[i] != xupp[i]) {
        cout << "xupp: " << i << "\t"
             << c_xupp[i] << "  !=  " << xupp[i] << endl;
      }
  }


  for (int i = 0; i < lenA; i++){
      if(c_A[i] != A[i]) {
        cout << "A: " << i << "\t"
             << c_A[i] << "  !=  " << A[i] << endl;
      }

      if(c_iAfun[i] != iAfun[i]) {
        cout << "iAfun: " << i << "\t"
             << c_iAfun[i] << "  !=  " << iAfun[i] << endl;
      }

      if(c_jAvar[i] != jAvar[i]) {
        cout << "jAvar: " << i << "\t"
             << c_jAvar[i] << "  !=  " << jAvar[i] << endl;
      }
  }


  for (int i = 0; i < lenG; i++){

      cout << "iGfun: " << i << "\t"
           << c_iGfun[i] << "\t\t\t   and  \t\t\t"
           << iGfun[i] << endl;

      cout << "jGvar: " << i << "\t"
           << c_jGvar[i] << "\t\t\t   and  \t\t\t"
          << jGvar[i] << endl;

      // if(c_iGfun[i] != iGfun[i])
      //     cout << "iGfun: " << i << "\t"
      //          << c_iGfun[i] << "\t\t\t  !=  \t\t\t"
      //          << iGfun[i] << endl;

      // if(c_jGvar[i] != jGvar[i])
      //     cout << "jGvar: " << i << "\t"
      //          << c_jGvar[i] << "\t\t\t  !=  \t\t\t"
      //          << jGvar[i] << endl;

  }


  for(int i = 0; i < lencw; i++){
      if (c_cw[i]  != cw[i]){
          cout << " ---  CW --- i = " << i << "\t "
               << c_cw[i] << " != " << cw[i] << endl;
      }
  }


  for(int i = 0; i < lenrw; i++){
      //rw[i] = c_rw[i];
      if (c_rw[i]  != rw[i]){
          cout << " ---  RW --- i = " << i << "\t "
               << c_rw[i] << " != " << rw[i] << endl;

      }
  }


  for(int i = 0; i < leniw; i++){
      //iw[i] = c_iw[i];
      if (c_iw[i]  != iw[i]){
          cout << " ---  IW --- i = " << i << "\t "
               << c_iw[i] << " != " << iw[i] << endl;

      }
  }


  //  inform = c_inform;
  if (c_inform != inform){
      cout << "inform: " << c_inform
           << "  !=  " << inform << endl;
  }

  if (c_neF != neF){
      cout << "neF: " << c_neF
           << "  !=  " << neF << endl;
  }

  if (c_n != n){
      cout << "n: " << c_n
          << "  !=  " << n << endl;
  }

  if (c_lenA != lenA){
      cout << "lenA: " << c_lenA
           << "  !=  " << lenA << endl;
  }

  if (c_lenG != lenG){
      cout << "lenG: " << c_lenG
           << "  !=  " << lenG << endl;
  }

  if (c_neA != neA){
      cout << "neA: " << c_neA
           << "  !=  " << neA << endl;
  }

  if (c_neG != neG){
      cout << "neG: " << c_neG
           << "  !=  " << neG << endl;
  }

  if (c_lencw != lencw){
      cout << "lencw: " << c_lencw
           << "  !=  " << lencw << endl;
  }

  if (c_leniw != leniw){
      cout << "leniw: " << c_leniw
           << "  !=  " << leniw << endl;
  }

  if (c_lenrw != lenrw){
      cout << "lenrw: " << c_lenrw
           << "  !=  " << lenrw << endl;
  }


  // miniw = c_miniw;
  // mincw = c_mincw;
  // minrw = c_minrw;

  if (c_mincw != mincw){
      cout << "mincw: " << c_mincw
           << "  !=  " << mincw << endl;
  }

  if (c_miniw != miniw){
      cout << "miniw: " << c_miniw
           << "  !=  " << miniw << endl;
  }

  if (c_minrw != minrw){
      cout << "minrw: " << c_minrw
           << "  !=  " << minrw << endl;
  }


  // mincw=500;
  // miniw=503;
  // minrw=521;
  cout << " ---------- ---------- ---------- ---------- " << endl;
  cout << " ---------- ---------- ---------- ---------- " << endl;

  // snjac_ will generate fortran style arrays.
  fortranStyleAG = 1;
}


int SNOPTHandler::memoryGuess_(integer &amincw,
                               integer &aminiw,
                               integer &aminrw) {
  int memoryGuess = 0;
  integer nxname = 1;
  integer nfname = 1;


  if (neA < 0) {
    neA = n * neF;
    memoryGuess = 1;
  }

  if (neG < 0) {
    neG = n * neF;
    memoryGuess = 1;
  }


  snmema_(&inform, &neF, &n, &nxname, &nfname, &neA, &neG,
          &amincw, &aminiw, &aminrw, cw, &lencw, iw, &leniw,
          rw, &lenrw, 8 * 500);

  return memoryGuess;
}


void SNOPTHandler::init_() {
  initCalled = 1;
  sninit_(&iPrint, &iSumm, cw,
          &lencw, iw, &leniw,
          rw, &lenrw, 8 * 500);
}


void SNOPTHandler::setParameter(char *stropt) {
  assert(initCalled == 1);

  integer iPrt = 0; // suppresses printing
  integer iSum = 0;
  ftnlen stropt_len = strlen(stropt);


  snset_(stropt, &iPrt, &iSum, &inform,
         cw, &lencw, iw, &leniw, rw,
         &lenrw, stropt_len, 8 * 500);
}


void SNOPTHandler::getParameter(char *stroptin,
                                char *stroptout) {
  assert(initCalled == 1);

  integer iPrt = 0;
  integer iSum = 0;
  ftnlen stroptin_len = strlen(stroptin);
  ftnlen stroptout_len = strlen(stroptout);


  sngetc_(stroptin, stroptout, &inform,  cw,
          &lencw, iw, &leniw, rw, &lenrw,
          stroptin_len, stroptout_len, 8 * 500);
}


void SNOPTHandler::setIntParameter(const char *stropt,
                                   integer opt) {
  assert(initCalled == 1);

  integer iPrt = 0; // suppresses printing
  integer iSum = 0;
  ftnlen stropt_len = strlen(stropt);

  char *cstr = const_cast<char *>(stropt);


  snseti_(cstr, &opt, &iPrt, &iSum, &inform,
          cw, &lencw, iw, &leniw, rw, &lenrw,
          stropt_len, 8 * 500);
}


void SNOPTHandler::getIntParameter(const char *stropt,
                                   integer &opt) {

  assert(initCalled == 1);

  integer iPrt = 0;
  integer iSum = 0;
  ftnlen stropt_len = strlen(stropt);

  char *cstr = const_cast<char *>(stropt);


  sngeti_(cstr, &opt, &inform, cw, &lencw,
          iw, &leniw, rw, &lenrw, stropt_len,
          8 * 500);
}


void SNOPTHandler::setRealParameter(const char *stropt,
                                    double opt) {
  assert(initCalled == 1);

  integer iPrt = 0; // suppresses printing
  integer iSum = 0;
  ftnlen stropt_len = strlen(stropt);

  char *cstr = const_cast<char *>(stropt);


  snsetr_(cstr, &opt, &iPrt, &iSum, &inform,
          cw, &lencw, iw, &leniw, rw, &lenrw,
          stropt_len, 8 * 500);
}


void SNOPTHandler::getRealParameter(const char *stropt,
                                    double &opt) {
  assert(initCalled == 1);
  integer iPrt = 0; // suppresses printing
  integer iSum = 0;
  ftnlen stropt_len = strlen(stropt);

  char *cstr = const_cast<char *>(stropt);


  sngetr_(cstr, &opt, &inform, cw, &lencw, iw,
          &leniw, rw, &lenrw, stropt_len, 8 * 500);
}


void SNOPTHandler::solve(integer starttype,
                         vector<double> &xsol,
                         vector<double> &fsol) {

  assert(initCalled == 1);
  // Ensures all user data initialized.
  areAllUserDataSet_();


  // Unlike snjac_ we also need neA and neG to be set.
  if (neA == -1 || neG == -1) {
    throw("neA and neG must be set before"
          "calling SNOPTOptimizer::solve()");
  }


  ftnlen npname = strlen(Prob);
  integer nS, nInf;
  double sInf;
  this->increment_(); // Convert array entries to Fortran style
  this->setMemory_();


  // this->computeJac();
  // neG = 0;
  // lenG = 0;
  snopta_(&starttype,
          &neF, &n,
          &nxnames, &nFnames,
          &ObjAdd, &ObjRow,
          Prob,
          usrfun, iAfun, jAvar, &lenA,
          &neA, A,
          iGfun, jGvar, &lenG, &neG,
          xlow, xupp, xnames, Flow,
          Fupp, Fnames, x, xstate,
          xmul, F, Fstate, Fmul,
          &inform, &mincw, &miniw,
          &minrw, &nS, &nInf, &sInf,
          cw, &lencw, iw, &leniw, rw, &lenrw,
          cw, &lencw, iw, &leniw, rw, &lenrw,
          npname,
          8 * short(nxnames),
          8 * short(nFnames),
          8 * 500,
          8 * 500);

  exitCode_ = iw[424-1];


  xsol.resize(n);
  memcpy(&xsol[0], x, n * sizeof(double));

  fsol.resize(neF);
  memcpy(&fsol[0], F, neF * sizeof(double));

  this->decrement_();  // Convert array entries to C style
}


void SNOPTHandler::setPrintFile(const char *aprintname) {
  // assert( initCalled = 1 );


// #ifdef UNIX
  strcpy(printname, aprintname);
// #else
//  strcpy_s(printname, aprintname);
// #endif
  prnt_len = strlen(printname);
  // snopenappend_( &iPrint, printname,   &inform, prnt_len );


  integer mode = 3; // read only

  openfile_(&iPrint, printname, &mode, prnt_len);
  // this->setIntParameter("Print file", iPrint);
}


void SNOPTHandler::setSummaryFile(const char *asummname) {
  //assert( initCalled = 1 );


// #ifdef UNIX
  strcpy(summname, asummname);
// #else
//  strcpy_s(summname, asummname);
// #endif
  summ_len = strlen(summname);


  integer mode = 3; // read only

  openfile_(&iSumm, summname, &mode, summ_len);
  // this->setIntParameter("Summary file", iSumm);
}


void SNOPTHandler::setSpecFile(const char *aspecname) {
  assert(initCalled == 1);


  // iSpecs = 4;
//#ifdef UNIX
  strcpy(specname, aspecname);
//#else
//  strcpy_s(specname, aspecname);
//#endif
  spec_len = strlen(specname);

  //  snfilewrapper_( specname, &iSpecs, &inform, cw, &lencw,
  //                  iw, &leniw, rw, &lenrw, spec_len, 8*lencw);


  integer mode = 1; // read only

  openfile_(&iSpecs, specname, &mode, spec_len);

  snspec_(&iSpecs, &inform, cw, &lencw, iw, &leniw, rw, &lenrw,
          (ftnlen) 8);


  switch (inform) {
    case 101:
    cout << "-- sqp_snopt: Specs file read successfully ! \n";
      break;

    case 131:
    cout << "-- sqp_snopt: Specs file specified\n";
      break;

    case 132:
    cout << "-- sqp_snopt: End-of-file encountered while looking for Specs file\n";
      break;

    case 133:
    cout << "-- sqp_snopt: End-of-file encountered before finding End\n";
      break;

    case 134:
    cout << "-- sqp_snopt: Endrun found before any valid sets of options\n";
      break;

    default:
    cout << "-- sqp_snopt: There were " << inform - 134
         << " errors while reading Specs\n";
      break;
  }


  // int snfilewrapper_(char *name__, integer *ispec, integer *
  // inform__, char *cw, integer *lencw, integer *iw, integer *leniw,
  // doublereal *rw, integer *lenrw, ftnlen name_len, ftnlen cw_len)

  if (inform != 101)
    cout << "-- sqp_snopt: Warning: unable to find specs file "
         << specname << "\n" << "\n";
}


void SNOPTHandler::setProblemSize(integer an, integer aneF) {
  n = an;
  neF = aneF;
}


void SNOPTHandler::setObjective(integer aObjRow,
                                double aObjAdd) {
  // checkSet = checkSet+2;
  ObjRow = aObjRow;
  ObjAdd = aObjAdd;
}


void SNOPTHandler::setA(integer alenA, integer *aiAfun,
                        integer *ajAvar, double *aA) {
  lenA = alenA;
  iAfun = aiAfun;
  jAvar = ajAvar;
  A = aA;
}


void SNOPTHandler::setNeA(integer aneA) {
  neA = aneA;
}


void SNOPTHandler::setG(integer alenG,
                        integer *aiGfun,
                        integer *ajGvar) {
  lenG = alenG;
  iGfun = aiGfun;
  jGvar = ajGvar;
}


void SNOPTHandler::setNeG(integer aneG) {
  neG = aneG;
}


void SNOPTHandler::setX(double *ax, double *axlow,
                        double *axupp, double *axmul,
                        integer *axstate) {
  x = ax;
  xlow = axlow;
  xupp = axupp;
  xmul = axmul;
  xstate = axstate;
}


void SNOPTHandler::setF(double *aF, double *aFlow,
                        double *aFupp, double *aFmul,
                        integer *aFstate) {
  F = aF;
  Flow = aFlow;
  Fupp = aFupp;
  Fmul = aFmul;
  Fstate = aFstate;
}


void SNOPTHandler::setXNames(char *axnames,
                             integer anxnames) {
  xnames = axnames;
  nxnames = anxnames;
}


void SNOPTHandler::setFNames(char *aFnames,
                             integer anFnames) {
  Fnames = aFnames;
  nFnames = anFnames;
}


void SNOPTHandler::setProbName(const char *aProb) {
  sprintf(Prob, "%8s", aProb );
}


void SNOPTHandler::setUserFun(My_fp ausrfun) {
  usrfun = ausrfun;
}


void
SNOPTHandler::setProblem(integer nvariables,
                         integer nproblemFunctions,
                         integer objectiveRow,
                         integer nA,
                         integer *iAfun,
                         integer *jAvar,
                         double *A,
                         integer lenG,
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
                         integer neA,
                         integer neG,
                         My_fp userfun,
                         bool haveDerivatives,
                         integer maxiterations,
                         char *xnames,
                         integer nxnames,
                         char *Fnames,
                         integer nFnames,
                         const char *problemname,
                         const char *printfile,
                         const char *specfile) {


  setProblemSize(nvariables, nproblemFunctions);
  setObjective(objectiveRow);

  setA(nA, iAfun, jAvar, A);
  setG(lenG, iGfun, jGvar);
  setX(x, xlow, xupp, xmul, xstate);
  setF(F, Flow, Fupp, Fmul, Fstate);


  setXNames(xnames, nxnames);
  setFNames(Fnames, nFnames);

  setNeA(neA);
  setNeG(neG);
  setUserFun(userfun);


  setProbName(problemname);
  setPrintFile(printfile);

  //  setSpecFile( specfile );
  //  setIntParameter( "Derivative option", integer(haveDerivatives) );
  //  setIntParameter( "Major Iteration limit", maxiterations );
}

}  // namespace Optimization
}  // namespace Optimizers