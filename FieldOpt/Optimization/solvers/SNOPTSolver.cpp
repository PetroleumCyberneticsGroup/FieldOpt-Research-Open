/***********************************************************
Created by bellout on 11.03.19
Copyright (C) 2019- Mathias Bellout
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

#include "SNOPTSolver.h"

namespace Optimization {
namespace Optimizers {

static Optimization::Optimizer::EmbeddedProblem *prob_ = nullptr;

// Rosenbrock cost function def
// ╦═╗  ╔═╗  ╔═╗  ╔═╗  ╔╗╔  ╔╗   ╦═╗  ╔═╗  ╔═╗  ╦╔═
// ╠╦╝  ║ ║  ╚═╗  ║╣   ║║║  ╠╩╗  ╠╦╝  ║ ║  ║    ╠╩╗
// ╩╚═  ╚═╝  ╚═╝  ╚═╝  ╝╚╝  ╚═╝  ╩╚═  ╚═╝  ╚═╝  ╩ ╩
#ifdef __cplusplus
extern "C" {
#endif
int Rosenbrock_(integer *Status, integer *n, double x[],
                integer *needF,  integer *neF, double F[],
                integer *needG,  integer *neG, double G[],
                char *cu, integer *lencu,
                integer iu[], integer *leniu,
                double ru[], integer *lenru);
#ifdef __cplusplus
}
#endif

int Rosenbrock_(integer  *Status, integer *n,    double x[],
                integer  *needF,  integer *neF,  double F[],
                integer  *needG,  integer *neG,  double G[],
                char     *cu,     integer *lencu,
                integer  iu[],    integer *leniu,
                double   ru[],    integer *lenru ) {


  int m = (int)*neF - 1;  // # of constraints

  if (*needF > 0) {
    // Rosenbrock testproblem
    F[0] = (1.0 - x[0]) * (1.0 - x[0])
      + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);

    if (m) {
      // Nonlinear constraint: unit disk
      F[1] = x[0] * x[0] + x[1] * x[1];
    }
  }

  if (*needG > 0) {
    // dF/dx0, dF/dx1
    G[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
    G[1] = 200.0 * (x[1] - x[0] * x[0]);

    if (m) {
      // dC/dx0, dC/dx1
      G[2] = 2.0 * x[0];
      G[3] = 2.0 * x[1];
    }
  }
  return 0;
}

// Cost function def
// ╔╦╗  ╦═╗  ╔╦╗  ╔═╗  ╔╦╗  ╔═╗
//  ║   ╠╦╝  ║║║  ║ ║   ║║  ╠═╣
//  ╩   ╩╚═  ╩ ╩  ╚═╝  ═╩╝  ╩ ╩
#ifdef __cplusplus
extern "C" {
#endif
int TRMod_A_(integer *Status, integer *n, doublereal x[],
             integer *needF, integer *neF, doublereal F[],
             integer *needG, integer *neG, doublereal G[],
             char *cu, integer *lencu,
             integer iu[], integer *leniu,
             doublereal ru[], integer *lenru);
#ifdef __cplusplus
}
#endif


int TRMod_A_(integer *Status, integer *n, doublereal x[],
             integer *needF, integer *neF, doublereal F[],
             integer *needG, integer *neG, doublereal G[],
             char *cu, integer *lencu,
             integer iu[], integer *leniu,
             doublereal ru[], integer *lenru) {

  bool dbg = false;

  double c = prob_->getTrC();
  Eigen::VectorXd g = prob_->getTrG();
  Eigen::MatrixXd H = prob_->getTrH();

  // x -> eigen format
  Eigen::VectorXd xv(*n);
  for (int ii = 0; ii < *n; ii++) {
    xv(ii) = x[ii];
  }

  if (dbg) {
    cout << "Debug c, g, gT, H:";
    cout << "c = \n" << c << endl;
    cout << "g = \n" << g << endl;
    cout << "H = \n" << H << endl;
    cout << "xv = \n" << xv << endl;
  }

  int m = (int)*neF - 1;  // # of constraints

  if (*needF > 0) {
    F[0] = c + g.transpose()*xv + 0.5*xv.transpose()*H*xv;


    if (m) {
      // Add nonlinear constraints to tr prob
      // F[1] = ...
      // F[2] = ...
    }
  }

  if (*needG > 0) {
    Eigen::VectorXd gv = g.transpose() + xv.transpose()*H;
    if (dbg) {
      cout << "gv = \n" << gv << endl;
    }

    for (int ii = 0; ii < *n; ii++) {
      G[ii] = gv(ii);
    }

    if (m) {
      // Standalone nonlinear constraint
      int n_nlc = 1; // # of standalone nlc
      // G[(*n) + n_nlc] = ...  // standalone constraint #1
      n_nlc++;
      // G[(*n) + n_nlc] = ...  // standalone constraint #2

      // Matrix of nl constraints gradients (1 constraint / row)
      // # of constraints = # of rows
      Eigen::MatrixXd gcm;  // = ...
      // Derivatives of constraints
      for (int ii = 0; ii < gcm.rows(); ii++) { // # of c
        for (int jj = 0; jj < gcm.cols(); jj++) { // # of vars
          // G[(*n) + i + j] = ...
        }
      }
    }
  }
  return 0;
}

string prntVecXd(VectorXd v, string lbl) {
  stringstream ss; char buffer [100];
  ss <<"                  " << lbl <<  ": [ ";
  for (int ii = 0; ii < v.rows(); ii++) {
    if (ii == v.rows() - 1) {
      sprintf(buffer, "% 10.3e ]\n", v(ii));
    } else {
      sprintf(buffer, "% 10.3e ", v(ii));
    }
    ss << buffer;
  }
  return ss.str();
}

// Dbg constructs
string md_ = "Optimization/solvers";
string cl_ = "SNOPTSolver";
string im_ = "", wm_ = "", em_ = "";

// SQP cost function def
// ╔═╗  ╔═╗   ╔═╗
// ╚═╗  ║═╬╗  ╠═╝
// ╚═╝  ╚═╝╚  ╩
#ifdef __cplusplus
extern "C" {
#endif
int SQP_(integer *Status, integer *n, double x[],
         integer *needF,  integer *neF, double F[],
         integer *needG,  integer *neG, double G[],
         char *cu, integer *lencu,
         integer iu[], integer *leniu,
         double ru[], integer *lenru);
#ifdef __cplusplus
}
#endif

int SQP_(integer  *Status, integer *n,    double x[],
         integer  *needF,  integer *neF,  double F[],
         integer  *needG,  integer *neG,  double G[],
         char     *cu,     integer *lencu,
         integer  iu[],    integer *leniu,
         double   ru[],    integer *lenru ) {

  bool dbg = true;

  // x -> eigen format
  Eigen::VectorXd xv(*n);
  for (int ii = 0; ii < *n; ii++) {
    xv(ii) = x[ii];
  }

  // if (dbg) {
  //   cout << "Debug c, g, gT, H:";
  //   cout << "c = \n" << c << endl;
  //   cout << "g = \n" << g << endl;
  //   cout << "H = \n" << H << endl;
  //   cout << "xv = \n" << xv << endl;
  // }

  int m = (int)*neF - 1;  // # of constraints

  // double c = 0.0;
  double f = 0.0;
  Eigen::VectorXd g00 = Eigen::VectorXd::Zero(*n);
  cout << "xv: " << xv.transpose() <<  endl;

  prob_->tcase_->SetRealVarValues(xv, g00);
  prob_->model_->ApplyCase(prob_->tcase_);
  prob_->simulator_->Evaluate();

  // Rerun to populate var container with gradients (after eval)
  g00 = prob_->objf_->grad(); // save unscaled g for dbg
  if (g00.size() == 0) {
    em_ = "No grad read from sim. Probably grad assembly failed";
    em_ += "due to well shut down during simulation.";
    ext_warn(em_, md_, cl_);
  }
  cout << "g00: " << g00.transpose() <<  endl;
  prob_->tcase_->SetRealVarValues(xv, g00);

  // Switch b/e FO.computed obj or read-in value from ECL GRD file
  prob_->tcase_->set_objf_value(prob_->objf_->value()); // FO.obj
  // prob_->tcase_->set_objf_value(prob_->objf_->fval()); // DBG: ECL GRD file

  f = prob_->tcase_->objf_value();
  // Get scaled gradient
  auto grads = prob_->tcase_->GetRealVarGradx(prob_->model_->variables());
  // grads[0]: scaled g (if scaling=true)
  // grads[1]: scaled + type.normalized g

  // dbg: f
  stringstream ss; char buffer [100];
  sprintf(buffer, "f: [ % 9.4e ] \n", f);
  ss << buffer;
  ss << prntVecXd(xv, "x  ");
  ss << prntVecXd(g00, "g00");
  ss << prntVecXd(grads[0], "g0s");
  ss << prntVecXd(grads[1], "gsn");

  // print to file
  FILE * pFile;
  pFile = std::fopen("sqp-snopt.out", "a");
  std::fprintf(pFile, "%s", ss.str().c_str());
  fclose (pFile);


  if (*needF > 0) {
    if (dbg) {
      cout << "f = " << f << endl;
      cout << "xv = [ " << xv.transpose() << " ] " << endl;
    }

    // ┌─┐  ┌─┐  ┌┬┐  ┌─┐  ┬ ┬  ┌┬┐  ┌─┐    ╔═╗
    // │    │ │  │││  ├─┘  │ │   │   ├┤     ╠╣
    // └─┘  └─┘  ┴ ┴  ┴    └─┘   ┴   └─┘    ╚
    // F[0] = c + g.transpose()*xv + 0.5*xv.transpose()*H*xv;
    F[0] = f;

    if (m) {
      // Add nonlinear constraints to tr prob
      // F[1] = ...
      // F[2] = ...
    }
  }

  if (*needG > 0) {
    if (dbg) {
      cout << "g00 = [ " << g00.transpose() << " ] " << endl;
      cout << "g0s = [ " << grads[0].transpose() << " ] " << endl;
      cout << "gsn = [ " << grads[1].transpose() << " ] " << endl;
    }

    // ┌─┐  ┌─┐  ┌┬┐    ╔═╗
    // │ ┬  ├┤    │     ║ ╦
    // └─┘  └─┘   ┴     ╚═╝
    // Eigen::VectorXd gv = g.transpose() + xv.transpose()*H;
    // Eigen::VectorXd gv = g;

    VectorXd g;
    if (prob_->getSeto()->scaleVars() == 0) {
      g = g00;
    } else if (prob_->getSeto()->scaleVars() == 1) {
      g = grads[0];
    } else if (prob_->getSeto()->scaleVars() == 2) {
      g = grads[1];
    }

    for (int ii = 0; ii < *n; ii++) {
      G[ii] = g(ii);
    }

    if (m) {
      // Standalone nonlinear constraint
      int n_nlc = 1; // # of standalone nlc
      // G[(*n) + n_nlc] = ...  // standalone constraint #1
      n_nlc++;
      // G[(*n) + n_nlc] = ...  // standalone constraint #2

      // Matrix of nl constraints gradients (1 constraint / row)
      // # of constraints = # of rows
      Eigen::MatrixXd gcm;  // = ...
      // Derivatives of constraints
      for (int ii = 0; ii < gcm.rows(); ii++) { // # of c
        for (int jj = 0; jj < gcm.cols(); jj++) { // # of vars
          // G[(*n) + i + j] = ...
        }
      }
    }
  }
  return 0;
}

SNOPTSolver::SNOPTSolver() {
  loadSNOPT();
}

SNOPTSolver::~SNOPTSolver() {
  delete[] iGfun_;
  delete[] jGvar_;
  delete[] x_;
  delete[] xlow_;
  delete[] xupp_;
  delete[] xmul_;
  delete[] xstate_;
  delete[] F_;
  delete[] Flow_;
  delete[] Fupp_;
  delete[] Fmul_;
  delete[] Fstate_;
  delete[] xnames_;
  delete[] Fnames_;
}

// ╔═╗  ╔═╗  ╔╦╗  ╦ ╦  ╔═╗    ╔═╗  ╔╗╔  ╔═╗  ╔═╗  ╔╦╗    ╔═╗  ╔═╗  ╦    ╦  ╦  ╔═╗  ╦═╗
// ╚═╗  ║╣    ║   ║ ║  ╠═╝    ╚═╗  ║║║  ║ ║  ╠═╝   ║     ╚═╗  ║ ║  ║    ╚╗╔╝  ║╣   ╠╦╝
// ╚═╝  ╚═╝   ╩   ╚═╝  ╩      ╚═╝  ╝╚╝  ╚═╝  ╩     ╩     ╚═╝  ╚═╝  ╩═╝   ╚╝   ╚═╝  ╩╚═
void SNOPTSolver::setUpSNOPTSolver(Optimization::Optimizer::EmbeddedProblem &prob) {

  if (prob.getProbName() == "TRMod_A") {
    subprobTRModA(prob);
    SNOPTRun_ = "TRMod_A";

  } else  if (prob.getProbName() == "TRMod_B") {
    SNOPTRun_ = "TRMod_B";

  } else  if (prob.getProbName() == "SQP") {
    subprobSQP(prob);
    SNOPTRun_ = "SQP";

  } else if (prob.getProbName() == "Rosenbrock") {
    subprobRosenbrock(prob);
    SNOPTRun_ = "Rosenbrock";
  }

  // For Reference: Initialization of SNOPTHandler pointer
  // Necessary to call destructor explicitly (see below),
  // otherwise we get: "Fortran runtime error: Cannot change
  // STATUS parameter in OPEN statement" [Issue avoided if
  // setting up SNOPTHandler within initSNOPTHandler(), then
  // passing SNOPTHandler object as reference]
  // SNOPTHandler_->~SNOPTHandler();
  // SNOPTHandler_ = nullptr;

  prnt_file_ = SNOPTRun_ + ".opt.prnt";
  smry_file_ = SNOPTRun_ + ".opt.smry";
  optn_file_ = SNOPTRun_ + ".opt.optn";

  auto d = prob.getSeto()->paths()->GetPath(Paths::SIM_EXEC_DIR);
  optn_file_ = d + "/" + optn_file_;

  SNOPTHandler snoptHandler = initSNOPTHandler();
  setSNOPTOptions(snoptHandler, prob);

  // Pointer to access prob data in SNOPTusrFG function
  prob_ = &prob;

  integer Cold = 0, Basis = 1, Warm = 2;

  vector<double> xsol;
  vector<double> fsol;

  snoptHandler.solve(Cold, xsol, fsol);

  prob.setXSol(xsol);
  prob.setFSol(fsol);

  if (prob.getSeto()->verbParams().vOPT >= 1) {
    cout << "sqp_snopt Exit code: " + std::to_string(snoptHandler.getExitCode()) << endl;
  }
  prob.setSNOPTExitCode(snoptHandler.getExitCode());
}

// ╔═╗  ╦ ╦  ╔╗   ╔═╗  ╦═╗  ╔═╗  ╔╗     ╔═╗  ╔═╗   ╔═╗
// ╚═╗  ║ ║  ╠╩╗  ╠═╝  ╠╦╝  ║ ║  ╠╩╗    ╚═╗  ║═╬╗  ╠═╝
// ╚═╝  ╚═╝  ╚═╝  ╩    ╩╚═  ╚═╝  ╚═╝    ╚═╝  ╚═╝╚  ╩
void SNOPTSolver::subprobSQP(Optimization::Optimizer::EmbeddedProblem &prob) {

//  cout << "Specs for TRMod subproblem A" << endl;

  n_ = prob.getNunVars(); // # of variables
  m_ = prob.getNunNnlConst(); // # of nonlinear c

  // Total # of constraints + objective (length of F vector)
  neF_ = m_ + 1; // l = 2
  // # of linear constraints
  lenA_ = prob.getNumLinConst();

  // Length of grad vector (dim of iGfun jGvar arrays)
  // --> n: (length of obj.grad = # of vars) +
  // --> m * n: (# of constraints) * (# of vars)
  lenG_ = n_ + m_ * n_;

  // The objective could be any component of F; here the
  // objective is specified as the first component of F
  objRow_ = 0;

  // Add nothing to objective for output purposes
  objAdd_ = 0.0;

  // Allocate memory

  // Objective
  F_ = new double[neF_];
  Flow_ = new double[neF_];
  Fupp_ = new double[neF_];
  Fmul_ = new double[neF_];
  Fstate_ = new integer[neF_];

  // Variables
  x_ = new double[n_];
  xlow_ = new double[n_];
  xupp_ = new double[n_];
  xmul_ = new double[n_];  // Lagrange mulpliers
  xstate_ = new integer[n_];  // state of variables

  iGfun_ = new integer[lenG_];
  jGvar_ = new integer[lenG_];

  // Linear constraints

  // iAfun_ = new integer[lenA_];
  // jAvar_ = new integer[lenA_];
  // A_ = new double[lenA_];

  // No linear constraints
  iAfun_ = nullptr;
  jAvar_ = nullptr;
  A_     = nullptr;

  xnames_ = new char[nxnames_ * 8];
  Fnames_ = new char[nxnames_ * 8];

  // Bounds objective
  Flow_[0] = prob.getFLb()(0);
  Fupp_[0] = prob.getFUb()(0);
  // Flow_[0] = -infinity_;
  // Fupp_[0] = infinity_;

  // Bounds nonlinear constraints
  for (int ii = 1; ii < neF_; ii++) {
    Flow_[ii] = prob.getFLb()(ii);
    Fupp_[ii] = prob.getFUb()(ii);
    // Flow_[ii] = -infinity_;
    // Fupp_[ii] = infinity_;
  }

  // Bounds variables
  for (int ii = 0; ii < n_; ii++) {
    xlow_[ii] = prob.getXLb()(ii);
    xupp_[ii] = prob.getXUb()(ii);
    // xlow_[ii] = -infinity_;
    // xupp_[ii] = infinity_;
  }

  // Initial point
  for (int ii = 0; ii < n_; ii++) {
    x_[ii] =  prob.getXInit()(ii);
    // x_[ii] = 0.0;
  }

  // Gradient indices of objective
  for (int ii = 0; ii < n_; ii++) {
    iGfun_[ii] = 0;
    jGvar_[ii] = ii;

    // ilc: loop over # of (linear) constraints
    for (int ilc = 0; ilc < lenA_; ilc++) {
      iAfun_[ii + ilc * n_] = m_ + 1 + ilc; // <- sure about this m_?
      jAvar_[ii + ilc * n_] = ii;
      A_[ii + ilc * n_] = 0; // <- coeff. of linear constraints
      //TODO If linear c: introduce coefficients here
    }
  }

  if (m_ != 0) {  // if having nonlinear constraints
    // Fill in constraint indices
    for (int jj = 1; jj <= m_; jj++) {
      for (int ii = 0; ii < n_; ii++) {
        iGfun_[ii + jj * n_] = jj;
        jGvar_[ii + jj * n_] = ii;
      }
    }
  }

  // Nonzero structure of the Jacobian
  neG_ = lenG_;
  neA_ = lenA_;

  // If initial guess provided then states should be zero
  for (int ii = 0; ii < n_; ii++) {
    xstate_[ii] = 0;
  }

  // Fmul is the vector with the estimation of the Lagrange
  // Multipliers; it will be always zero except in very rare
  // cases of benchmarking performance with them set to some
  // initial guess
  for (int ii = 0; ii < neF_; ii++) {
    Fmul_[ii] = 0;
  }
  // TODO Extract value from sqp_snopt operation: Fmul, dual, ++
}

// ╔═╗  ╦ ╦  ╔╗   ╔═╗  ╦═╗  ╔═╗  ╔╗     ╦═╗  ╔═╗  ╔═╗  ╔═╗  ╔╗╔  ╔╗   ╦═╗  ╔═╗  ╔═╗  ╦╔═
// ╚═╗  ║ ║  ╠╩╗  ╠═╝  ╠╦╝  ║ ║  ╠╩╗    ╠╦╝  ║ ║  ╚═╗  ║╣   ║║║  ╠╩╗  ╠╦╝  ║ ║  ║    ╠╩╗
// ╚═╝  ╚═╝  ╚═╝  ╩    ╩╚═  ╚═╝  ╚═╝    ╩╚═  ╚═╝  ╚═╝  ╚═╝  ╝╚╝  ╚═╝  ╩╚═  ╚═╝  ╚═╝  ╩ ╩
void SNOPTSolver::subprobRosenbrock(Optimization::Optimizer::EmbeddedProblem &prob) {

  // cout << "Specs for Rosenbrock subproblem" << endl;

  n_ = 2; // # of variables
  m_ = 1; // # of nonlinear c

  // Total # of constraints + objective (length of F vector)
  neF_ = m_ + 1; // l = 2
  // # of linear constraints
  lenA_ = 0;

  // Length of grad vector (dim of iGfun jGvar arrays)
  // --> n: (length of obj.grad = # of vars) +
  // --> m * n: (# of constraints) * (# of vars)
  lenG_ = n_ + m_ * n_;

  // The objective could be any component of F; here the
  // objective is specified as the first component of F
  objRow_ = 0;

  // Add nothing to objective for output purposes
  objAdd_ = 0.0;

  // Allocate memory

  // Objective
  F_ = new double[neF_];
  Flow_ = new double[neF_];
  Fupp_ = new double[neF_];
  Fmul_ = new double[neF_];
  Fstate_ = new integer[neF_];

  // Variables
  x_ = new double[n_];
  xlow_ = new double[n_];
  xupp_ = new double[n_];
  xmul_ = new double[n_];  // Lagrange mulpliers
  xstate_ = new integer[n_];  // state of variables

  iGfun_ = new integer[lenG_];
  jGvar_ = new integer[lenG_];

  // Linear constraints

  // iAfun_ = new integer[lenA_];
  // jAvar_ = new integer[lenA_];
  // A_ = new double[lenA_];

  // No linear constraints
  iAfun_ = nullptr;
  jAvar_ = nullptr;
  A_     = nullptr;

  xnames_ = new char[nxnames_ * 8];
  Fnames_ = new char[nxnames_ * 8];

  // Bounds objective
  Flow_[0] = -infinity_;
  Fupp_[0] = infinity_;

  // Bounds nonlinear constraints
  for (int ii = 1; ii < neF_; ii++) {
    Flow_[ii] = -infinity_;
    Fupp_[ii] = infinity_;
  }

  // Bounds variables
  for (int ii = 0; ii < n_; ii++) {
    xlow_[ii] = -infinity_;
    xupp_[ii] = infinity_;
  }

  // Initial point
  for (int ii = 0; ii < n_; ii++) {
    x_[ii] = 0.0;
  }

  // Gradient indices of objective
  for (int ii = 0; ii < n_; ii++) {
    iGfun_[ii] = 0;
    jGvar_[ii] = ii;

    // ilc: loop over # of (linear) constraints
    for (int ilc = 0; ilc < lenA_; ilc++) {
      iAfun_[ii + ilc * n_] = m_ + 1 + ilc; // <- sure about this m_?
      jAvar_[ii + ilc * n_] = ii;
      A_[ii + ilc * n_] = 0; // <- coeff. of linear constraints
      //TODO If linear c: introduce coefficients here
    }
  }

  if (m_ != 0) {  // if having nonlinear constraints
    // Fill in constraint indices
    for (int jj = 1; jj <= m_; jj++) {
      for (int ii = 0; ii < n_; ii++) {
        iGfun_[ii + jj * n_] = jj;
        jGvar_[ii + jj * n_] = ii;
      }
    }
  }

  // Nonzero structure of the Jacobian
  neG_ = lenG_;
  neA_ = lenA_;

  // If initial guess provided then states should be zero
  for (int ii = 0; ii < n_; ii++) {
    xstate_[ii] = 0;
  }

  // Fmul: estimation of Lagrange Multipliers
  for (int ii = 0; ii < neF_; ii++) {
    Fmul_[ii] = 0;
  }
}

// ╔═╗  ╦ ╦  ╔╗   ╔═╗  ╦═╗  ╔═╗  ╔╗     ╔╦╗  ╦═╗  ╔╦╗  ╔═╗  ╔╦╗  ╔═╗
// ╚═╗  ║ ║  ╠╩╗  ╠═╝  ╠╦╝  ║ ║  ╠╩╗     ║   ╠╦╝  ║║║  ║ ║   ║║  ╠═╣
// ╚═╝  ╚═╝  ╚═╝  ╩    ╩╚═  ╚═╝  ╚═╝     ╩   ╩╚═  ╩ ╩  ╚═╝  ═╩╝  ╩ ╩
void SNOPTSolver::subprobTRModA(Optimization::Optimizer::EmbeddedProblem &prob) {

//  cout << "Specs for TRMod subproblem A" << endl;

  n_ = prob.getNunVars(); // # of variables
  m_ = prob.getNunNnlConst(); // # of nonlinear c

  // Total # of constraints + objective (length of F vector)
  neF_ = m_ + 1; // l = 2
  // # of linear constraints
  lenA_ = prob.getNumLinConst();

  // Length of grad vector (dim of iGfun jGvar arrays)
  // --> n: (length of obj.grad = # of vars) +
  // --> m * n: (# of constraints) * (# of vars)
  lenG_ = n_ + m_ * n_;

  // The objective could be any component of F; here the
  // objective is specified as the first component of F
  objRow_ = 0;

  // Add nothing to objective for output purposes
  objAdd_ = 0.0;

  // Allocate memory

  // Objective
  F_ = new double[neF_];
  Flow_ = new double[neF_];
  Fupp_ = new double[neF_];
  Fmul_ = new double[neF_];
  Fstate_ = new integer[neF_];

  // Variables
  x_ = new double[n_];
  xlow_ = new double[n_];
  xupp_ = new double[n_];
  xmul_ = new double[n_];  // Lagrange mulpliers
  xstate_ = new integer[n_];  // state of variables

  iGfun_ = new integer[lenG_];
  jGvar_ = new integer[lenG_];

  // Linear constraints

  // iAfun_ = new integer[lenA_];
  // jAvar_ = new integer[lenA_];
  // A_ = new double[lenA_];

  // No linear constraints
  iAfun_ = nullptr;
  jAvar_ = nullptr;
  A_     = nullptr;

  xnames_ = new char[nxnames_ * 8];
  Fnames_ = new char[nxnames_ * 8];

  // Bounds objective
  Flow_[0] = prob.getFLb()(0);
  Fupp_[0] = prob.getFUb()(0);
  // Flow_[0] = -infinity_;
  // Fupp_[0] = infinity_;

  // Bounds nonlinear constraints
  for (int ii = 1; ii < neF_; ii++) {
    Flow_[ii] = prob.getFLb()(ii);
    Fupp_[ii] = prob.getFUb()(ii);
    // Flow_[ii] = -infinity_;
    // Fupp_[ii] = infinity_;
  }

  // Bounds variables
  for (int ii = 0; ii < n_; ii++) {
    xlow_[ii] = prob.getXLb()(ii);
    xupp_[ii] = prob.getXUb()(ii);
    // xlow_[ii] = -infinity_;
    // xupp_[ii] = infinity_;
  }

  // Initial point
  for (int ii = 0; ii < n_; ii++) {
    x_[ii] =  prob.getXInit()(ii);
    // x_[ii] = 0.0;
  }

  // Gradient indices of objective
  for (int ii = 0; ii < n_; ii++) {
    iGfun_[ii] = 0;
    jGvar_[ii] = ii;

    // ilc: loop over # of (linear) constraints
    for (int ilc = 0; ilc < lenA_; ilc++) {
      iAfun_[ii + ilc * n_] = m_ + 1 + ilc; // <- sure about this m_?
      jAvar_[ii + ilc * n_] = ii;
      A_[ii + ilc * n_] = 0; // <- coeff. of linear constraints
      //TODO If linear c: introduce coefficients here
    }
  }

  if (m_ != 0) {  // if having nonlinear constraints
    // Fill in constraint indices
    for (int jj = 1; jj <= m_; jj++) {
      for (int ii = 0; ii < n_; ii++) {
        iGfun_[ii + jj * n_] = jj;
        jGvar_[ii + jj * n_] = ii;
      }
    }
  }

  // Nonzero structure of the Jacobian
  neG_ = lenG_;
  neA_ = lenA_;

  // If initial guess provided then states should be zero
  for (int ii = 0; ii < n_; ii++) {
    xstate_[ii] = 0;
  }

  // Fmul is the vector with the estimation of the Lagrange
  // Multipliers; it will be always zero except in very rare
  // cases of benchmarking performance with them set to some
  // initial guess
  for (int ii = 0; ii < neF_; ii++) {
    Fmul_[ii] = 0;
  }
  // TODO Extract value from sqp_snopt operation: Fmul, dual, ++
}

// ╔═╗  ╔═╗  ╔╦╗    ╔═╗  ╔╗╔  ╔═╗  ╔═╗  ╔╦╗    ╔═╗  ╔═╗  ╔╦╗  ╦  ╔═╗  ╔╗╔  ╔═╗
// ╚═╗  ║╣    ║     ╚═╗  ║║║  ║ ║  ╠═╝   ║     ║ ║  ╠═╝   ║   ║  ║ ║  ║║║  ╚═╗
// ╚═╝  ╚═╝   ╩     ╚═╝  ╝╚╝  ╚═╝  ╩     ╩     ╚═╝  ╩     ╩   ╩  ╚═╝  ╝╚╝  ╚═╝
void SNOPTSolver::setSNOPTOptions(SNOPTHandler &H, Optimization::Optimizer::EmbeddedProblem &prob) {

  H.setProbName(prob.getProbName().c_str());
  H.setParameter((char*)"Minimize");

  H.setProblemSize( n_, neF_ );
  H.setObjective  ( objRow_ );
  H.setA          ( lenA_, iAfun_, jAvar_, A_ );
  H.setG          ( lenG_, iGfun_, jGvar_ );
  H.setX          ( x_, xlow_, xupp_, xmul_, xstate_ );
  H.setF          ( F_, Flow_, Fupp_, Fmul_, Fstate_ );
  H.setXNames     ( xnames_, nxnames_ );
  H.setFNames     ( Fnames_, nFnames_ );
  H.setNeA        ( neA_ );
  H.setNeG        ( neG_ );

  auto p = prob.getSeto()->parameters();
  p.sqp_ftol;

  // ┌┬┐  ┬─┐  ┌┬┐  ┌─┐  ┌┬┐    ┌─┐
  //  │   ├┬┘  │││  │ │   ││    ├─┤
  //  ┴   ┴└─  ┴ ┴  └─┘  ─┴┘    ┴ ┴
  if (prob.getProbName() == "TRMod_A") {
    H.setUserFun(TRMod_A_);
    H.setParameter((char *) "Scale Print");
    H.setParameter((char *) "Solution  Yes");
    H.setIntParameter("Minor print level", 10);
    H.setIntParameter("Scale option", 1);
    H.setParameter((char *) "Scale Print");
    H.setParameter((char *) "Solution  Yes");

    // ┌┬┐  ┬─┐  ┌┬┐  ┌─┐  ┌┬┐    ┌┐
    //  │   ├┬┘  │││  │ │   ││    ├┴┐
    //  ┴   ┴└─  ┴ ┴  └─┘  ─┴┘    └─┘
  } else  if (prob.getProbName() == "TRMod_B") {

    // ┌─┐  ┌─┐   ┌─┐
    // └─┐  │─┼┐  ├─┘
    // └─┘  └─┘└  ┴
  } else  if (prob.getProbName() == "SQP") {
    H.setUserFun(SQP_);

    // 5spot BC01/BC02
    H.setParameter((char*) "Maximize");
    H.setProbName((char*) "SQP");
    H.setParameter((char*) "System information Yes");
    H.setParameter((char *) "Solution  Yes");

    // H.setRealParameter("Linesearch tolerance",        0.9);
    // H.setRealParameter("Major feasibility tolerance", 1e-12);
    // H.setIntParameter("Major Iterations Limit", 100);

    H.setIntParameter("Verify level", -1);

    // H.setIntParameter("Scale option", 1); // not working
    H.setParameter((char *) "Scale Print");
    H.setParameter((char *) "Solution Yes");

    // H.setRealParameter("Major step limit", 100);
    H.setRealParameter("Major optimality tolerance", 1e-12);

    double f_prec = 1e-16;
    H.setRealParameter("Function precision", f_prec);
    H.setRealParameter("Difference interval",std::pow(f_prec, 1.0/2.0));
    H.setRealParameter("Central difference interval",std::pow(f_prec, 1.0/3.0));

    H.setIntParameter("Major print level", 1);
    H.setIntParameter("Minor print level", 1);

    H.setParameter((char*) "Print frequency 1");
    H.setParameter((char*) "Summary frequency 1");

    // ┬─┐  ┌─┐  ┌─┐  ┌─┐  ┌┐┌  ┌┐   ┬─┐  ┌─┐  ┌─┐  ┬┌─
    // ├┬┘  │ │  └─┐  ├┤   │││  ├┴┐  ├┬┘  │ │  │    ├┴┐
    // ┴└─  └─┘  └─┘  └─┘  ┘└┘  └─┘  ┴└─  └─┘  └─┘  ┴ ┴
  } else if (prob.getProbName() == "Rosenbrock") {
    H.setUserFun(Rosenbrock_);
    H.setParameter((char *) "Scale Print");
    H.setParameter((char *) "Solution  Yes");
    H.setIntParameter("Minor print level", 10);

    H.setIntParameter("Scale option", 1);
    H.setParameter((char *) "Scale Print");
    H.setParameter((char *) "Solution  Yes");
  }

  // ┌┬┐  ┌─┐  ┌─┐  ┌─┐  ┬ ┬  ┬    ┌┬┐
  //  ││  ├┤   ├┤   ├─┤  │ │  │     │
  // ─┴┘  └─┘  └    ┴ ┴  └─┘  ┴─┘   ┴
  // H.setParameter(                 "Backup basis file 0");
  // H.setRealParameter(      "Central difference interval",
  //                    2 * derivativeRelativePerturbation);

  // H.setIntParameter(              "Check frequency", 60);
  // H.setParameter("Cold Start                      Cold");

  // H.setParameter("Crash option                       3");
  // H.setParameter("Crash tolerance                  0.1");
  // H.setParameter("Derivative level                   3");

  // H.setParameter(     (char*)"Nonderivative linesearch");
  // H.setParameter(        (char*)"Derivative linesearch");
  // H.setIntParameter(             "Derivative option", 0);

  // H.setRealParameter("Difference interval",        ????);


  // H.setParameter("Dump file                          0");
  // H.setParameter("Elastic weight                1.0e+4");
  // H.setParameter("Expand frequency               10000");
  // H.setParameter("Factorization frequency           50");
  // H.setRealParameter("Function precision",        ?????);

  // H.setParameter("Hessian full memory");
  // H.setParameter("Hessian limited memory");

  // H.setIntParameter("Hessian frequency",           ????);
  // H.setIntParameter("Hessian updates",                0);

  // Does NOT work in the current version of sqp_snopt!!!
  // H.setIntParameter("Hessian flush",    1);

  // H.setParameter("Insert file                        0");
  // H.setRealParameter("Infinite bound",            ?????);

  // H.setParameter("Iterations limit");
  // H.setRealParameter("Linesearch tolerance",        0.9);
  // H.setParameter("Load file                          0");
  // H.setParameter("Log frequency                    100");
  // H.setParameter("LU factor tolerance             3.99");
  // H.setParameter("LU update tolerance             3.99");
  // H.setRealParameter("LU factor tolerance",        3.99);
  // H.setRealParameter("LU update tolerance",        3.99);
  // H.setParameter("LU partial pivoting");
  // H.setParameter("LU density tolerance             0.6");
  // H.setParameter("LU singularity tolerance     3.2e-11");

  // Target nonlinear constraint violation
  // H.setRealParameter("Major feasibility tolerance", 0.000001);
  // H.setIntParameter("Major Iterations Limit", 1000);

  // Target complementarity gap
  // H.setRealParameter("Major optimality tolerance", 0.0001);

  // H.setParameter("Major Print level  11111");  // 000001
  // H.setRealParameter("Major step limit", 0.2);
  // H.setIntParameter("Minor iterations limit", 200);


  // For satisfying the QP bounds
  // H.setRealParameter("Minor feasibility tolerance", ???);
  // H.setIntParameter("Minor print level", 10);

  // H.setParameter("New basis file                 0");
  // H.setParameter("New superbasics limit          99");
  // H.setParameter("Objective Row");
  // H.setParameter("Old basis file                 0");
  // H.setParameter("Partial price                  1");
  // H.setParameter("Pivot tolerance                3.7e-11");
  // H.setParameter("Print frequency                100");
  // H.setParameter("Proximal point method          1");
  // H.setParameter("QPSolver Cholesky");
  // H.setParameter("Reduced Hessian dimension");
  // H.setParameter("Save frequency                 100");

  // H.setIntParameter("Scale option", 1);
  // H.setParameter("Scale tolerance                0.9");
  // H.setParameter((char *) "Scale Print");
  // H.setParameter((char *) "Solution  Yes");

  // H.setParameter(   "Start Objective Check at Column 1");
  // H.setParameter(  "Start Constraint Check at Column 1");
  // H.setParameter(      "Stop Objective Check at Column");
  // H.setParameter(     "Stop Constraint Check at Column");
  // H.setParameter(                "Sticky parameters No");
  // H.setParameter(               "Summary frequency 100");
  // H.setParameter(                   "Superbasics limit");
  // H.setParameter(                 "Suppress parameters");
  // H.setParameter(        (char*)"System information No");
  // H.setParameter(                      "Timing level 3");
  // H.setRealParameter("Unbounded objective value 1.0e+15");
  // H.setParameter(         "Unbounded step size 1.0e+18");
  // H.setIntParameter(                 "Verify level", -1); // -1
  // H.setRealParameter(           "Violation limit", 1e-8); // 1e-8

}

// ╦  ╔╗╔  ╦  ╔╦╗    ╔═╗  ╔╗╔  ╔═╗  ╔═╗  ╔╦╗  ╦ ╦  ╔═╗  ╔╗╔  ╔╦╗  ╦    ╔═╗  ╦═╗
// ║  ║║║  ║   ║     ╚═╗  ║║║  ║ ║  ╠═╝   ║   ╠═╣  ╠═╣  ║║║   ║║  ║    ║╣   ╠╦╝
// ╩  ╝╚╝  ╩   ╩     ╚═╝  ╝╚╝  ╚═╝  ╩     ╩   ╩ ╩  ╩ ╩  ╝╚╝  ═╩╝  ╩═╝  ╚═╝  ╩╚═
SNOPTHandler SNOPTSolver::initSNOPTHandler() {

  if (false) {
    cout << "Initiate SNOPTHandler" << endl;
    cout << "prnt_file: " << prnt_file_ << endl;
    cout << "smry_file: " << smry_file_ << endl;
    cout << "optn_file: " << optn_file_ << endl;
  }

  SNOPTHandler snoptHandler(prnt_file_.c_str(),
                            smry_file_.c_str(),
                            optn_file_.c_str());

  return snoptHandler;
}


// ╦═╗  ╔═╗  ╔═╗  ╔═╗  ╔╦╗    ╔═╗  ╔╗╔  ╔═╗  ╔═╗  ╔╦╗    ╔═╗  ╔═╗  ╦    ╦  ╦  ╔═╗  ╦═╗
// ╠╦╝  ║╣   ╚═╗  ║╣    ║     ╚═╗  ║║║  ║ ║  ╠═╝   ║     ╚═╗  ║ ║  ║    ╚╗╔╝  ║╣   ╠╦╝
// ╩╚═  ╚═╝  ╚═╝  ╚═╝   ╩     ╚═╝  ╝╚╝  ╚═╝  ╩     ╩     ╚═╝  ╚═╝  ╩═╝   ╚╝   ╚═╝  ╩╚═
void SNOPTSolver::resetSNOPTSolver() {

  for (int i = 0; i < n_; i++) {
    Fstate_[i] = 0;
    xstate_[i] = 0;
    x_[i] = 0.0;
    xmul_[i] = 0;
  }

  for (int h = 0; h < neF_; h++) {
    F_[h] = 0.0;
    Fmul_[h] = 0.0;
  }
}

// ╦    ╔═╗  ╔═╗  ╔╦╗    ╔═╗  ╔╗╔  ╔═╗  ╔═╗  ╔╦╗
// ║    ║ ║  ╠═╣   ║║    ╚═╗  ║║║  ║ ║  ╠═╝   ║
// ╩═╝  ╚═╝  ╩ ╩  ═╩╝    ╚═╝  ╝╚╝  ╚═╝  ╩     ╩
bool SNOPTSolver::loadSNOPT(const string lib_name) {

// #ifdef NDEBUG

  if (LSL_isSNOPTLoaded()) {
    printf("\x1b[33mSNOPT is already loaded.\n\x1b[0m");
    return true;
  }


  char buf[256];
  int rc;
  if (lib_name.empty()) {
    rc = LSL_loadSNOPTLib(nullptr, buf, 255);
  } else {
    rc = LSL_loadSNOPTLib(lib_name.c_str(), buf, 255);
  }
  // printf("\x1b[33mSNOPT load exit # = %d.\n\x1b[0m", rc);


  if (rc) {
    string errmsg;
    errmsg = "Selected NLP solver sqp_snopt not available.\n"
             "Tried to obtain sqp_snopt from shared library \"";
    errmsg += LSL_SNOPTLibraryName();
    errmsg += "\", but the following error occured:\n";
    errmsg += buf;
    cout << errmsg << endl;
    return false;
  }
// #endif

  return true;
}

}  // namespace Optimizers
}  // namespace Optimization
