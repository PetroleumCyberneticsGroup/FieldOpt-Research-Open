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

#ifndef FOEX_DC_5SNOPTSOLVER_H
#define FOEX_DC_5SNOPTSOLVER_H

#include <Optimization/optimizer.h>

#include <ThirdParty/snopt/handlers/SNOPTHandler.h>
#include <ThirdParty/snopt/handlers/SNOPTLoader.h>

#include <Utilities/printer.hpp>

#include <Eigen/Core>

namespace Optimization {
namespace Optimizers {

class SNOPTSolver {

 public:
  SNOPTSolver();
  ~SNOPTSolver();

  void setUpSNOPTSolver(Optimization::Optimizer::EmbeddedProblem &prob);

 private:

  int n_;  // # of variables
  int m_;  // # of nonlinear constraints

  integer neF_;  // # of elements in F
  integer neG_;
  integer lenG_;
  integer objRow_;
  double objAdd_;

  integer *iAfun_ = nullptr;
  integer *jAvar_ = nullptr;
  double *A_ = nullptr;
  integer lenA_;
  integer neA_;

  integer *iGfun_ = nullptr;
  integer *jGvar_ = nullptr;

  double *x_ = nullptr;
  double *F_ = nullptr;
  double *Flow_ = nullptr;
  double *Fupp_ = nullptr;
  double *Fmul_ = nullptr;
  integer *Fstate_ = nullptr;
  char *xnames_ = nullptr;
  char *Fnames_ = nullptr;

  integer nxnames_ = 1;
  integer nFnames_ = 1;

  // Lower, upper bounds
  double *xlow_ = nullptr;
  double *xupp_ = nullptr;

  // Initial guess Lagrange multipliers
  double *xmul_ = nullptr;;

  // State of variables (whether the optimal
  // is likely to be on the boundary or not)
  integer *xstate_ = nullptr;

  // Value sqp_snopt considers as infinity
  const double infinity_ = 1e20;

  // Setup for SNOPTHandler
  // optn_file, smry_file and prnt_file file names should
  // be read from JSON and stored in settings_
  Settings::Optimizer *settings_;

  string SNOPTRun_ = "DEF";
  string prnt_file_ = SNOPTRun_ + ".opt.prnt";
  string smry_file_ = SNOPTRun_ + ".opt.smry";
  string optn_file_ = SNOPTRun_ + ".opt.optn";

  SNOPTHandler initSNOPTHandler();

  // Setup for subproblems
  void subprobTRModA(Optimization::Optimizer::EmbeddedProblem &prob);
  void subprobRosenbrock(Optimization::Optimizer::EmbeddedProblem &prob);
  void subprobSQP(Optimization::Optimizer::EmbeddedProblem &prob);

  // sqp_snopt functions / set sqp_snopt params
  bool loadSNOPT(string lib_name = "libsnopt-7.2.12.2.so");

  void setSNOPTOptions(
      SNOPTHandler &snoptHandler,
      Optimization::Optimizer::EmbeddedProblem &prob);

 public:
  void resetSNOPTSolver(); // Needed?
};

}  // namespace Optimizers
}  // namespace Optimization

#endif //FOEX_DC_5SNOPTSOLVER_H

// This class will find _one_ maximum of a quadratic
// function (specified by c_, g_ and H_) subject to
// some specified constraints.

// The constraints must be specified by whomever uses
// this function and they must be specified by editing
// the code explicitly (i.e, they cannot be set through
// function calls).

// The objective function and the constraints (except the
// simple/basic bounds) are put together into one set of
// equations:

// Let cl_i represent a linear constraint,
// and let cn_i represent a nonlinear one.

// Let n_l and m be the numbers of linear constraints
// and nonlinear constraints, respectively.

// F = [f_obj, cl_0, cl_1, ..., cl_n_l, cn_0, cn_1, ..., cn_m]^T;

// The inequalities for the constraints and
// the objective function are given by:
// Flow <= F <= Fupp

// If you don't have a limit, use the infinity_ *(+-1) value.

// maximize    f_obj = c_ + g_^T x + x^T H_ x

// subject to     Flow <= F <= Fupp
// 				        xlow <= x <= xupp

// Linear and nonlinear constraints are specified differently.
//
// Note!
// Because F contains both the objective functions and
// the constraints, the row-indices will start at row 1
// and not row 0.

// Linear
//  The linear constraints are specified through lenA,
//  iAfun, jAvar and A.

//  The first 2 are used because the matrix hatA
//  ( "lower <= hatA*x <= upper") might be very
//  sparse. An example will illustrate the usage:
//
// 	A = [ 0 1 0
// 		    5 0 6
// 		    0 0 8 ].
//
// 	lenA = 4; // There are four nonzero elements in A

//  iAfun[0] = 1;
// 	jAvar[0] = 1;
//  These two indices belong to the element (0,1) in A
//  (namely, the value 1). Now A[0] must be set to 1.
// 	A[0] = 1;

//  These two indices belong to the element (1,0) in A
//  (namely, the value 5). Now A[1] must be set to 5.
//	iAfun[1] = 2;
// 	jAvar[1] = 0;
// 	A[1] = 5;
//
//  These two indices belongs to the element (1,2) in A
//  (namely, the value 6). Now A[2] must be set to 6.
// 	iAfun[2] = 2;
// 	jAvar[2] = 2;
// 	A[2] = 6;
//
//  These two indices belongs to the element (2,2) in A
//  (namely, the value 8). Now A[3] must be set to 8.
//  iAfun[3] = 3;
// 	jAvar[3] = 2;
// 	A[3] = 8;
//
//
// Nonlinear
// 	The nonlinear constraints are specified through lenG,
//  neG, iGfun and jGvar. The actual G matrix is specified
//  in the userfunc.

// 	The G matrix contains both the derivative of the
//  objective function and the derivative of the nonlinear
//  constraints.

// 	The nonlinear constraints must also be specified by the
//  userfunc, and put into appropriate place in F:
// 	F[4]=x_1^2 + x_2^2 + x_3^2;

// 	Let's say that we now, in addition to the linear
//  constraints above, also has 1 nonlinear constraint.
//  ("lower <= x_1^2 + x_2^2 + x_3^2  <= upper").

// 	The partial derivatives with respect
//  to the variables will then be:
// 	G[y] = 2*x_1;
// 	G[y+1] = 2*x_2;
// 	G[y+2] = 2*x_3;

// 	Where y is the number of partial derivatives
//  of the objective functions.
//
//  If f_obj = x_1 + x_2 + x_3, then y = 3; and
//
// 	F[0] = x_1 + x_2 + x_3;
// 	G[0] = 1;
// 	G[1] = 1;
// 	G[2] = 1;
//
// 	Now we must specify iGfun and jGvar.
//
// 	//From the objective function:
// 	iGfun[0] = 0;
// 	jGvar[0] = 0;
// 	iGfun[1] = 0;
// 	jGvar[1] = 1;
// 	iGfun[2] = 0;
// 	jGvar[2] = 2;
//
//
//	//From the nonlinear constraints:
// 	iGfun[3] = 4; NOTE NOTE! The reason why this value is 4
//  is because the first row is for the objective function,
//  then  we have 3 linear constraints
// 	jGvar[3] = 0;
//
//	iGfun[4] = 4;
// 	jGvar[4] = 1;
//
//	iGfun[5] = 4;
// 	jGvar[5] = 2;
//
//	neG = lenG = 6;





