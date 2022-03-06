/***********************************************************
Sun May 02 2021 01:07:08 week 18 CET+0200
Copyright (C) 2021-
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

#ifndef FIELDOPT_OPTIMIZATION_OPTIMIZERS_SQP_SNOPT_SQP_SNOPT_H_
#define FIELDOPT_OPTIMIZATION_OPTIMIZERS_SQP_SNOPT_SQP_SNOPT_H_

#include "Optimization/optimizer.h"
#include <Optimization/solvers/SNOPTSolver.h>

using namespace Eigen;
using namespace std;

namespace Optimization {
namespace Optimizers {

using namespace Eigen;
using Model::Properties::VarPropContainer;
using Optimization::Constraints::ConstraintHandler;

using TC = Optimization::Optimizer::TerminationCondition;

class SQP_SNOPT : public Optimizer {

 public:
  SQP_SNOPT(Settings::Optimizer *seto,
            Case *bscs,
            Model::Model *model,
            Simulation::Simulator *simulator,
            Optimization::Objective::Objective *objf,
            Logger *logger, CaseHandler *case_handler = nullptr,
            Constraints::ConstraintHandler *constraint_handler = nullptr);

  TC IsFinished() override;

 protected:
  void iterate() override;
  void handleEvaluatedCase(Case *c) override {};

 private:
  // FO constructs
  Settings::Optimizer *seto_;
  Model::Properties::VarPropContainer *vars_;

  Case *bscs_;
  Model::Model *model_;
  Simulation::Simulator *simulator_;
  Optimization::Objective::Objective *objf_;

  // SNOPT solver
  SNOPTSolver *SNOPTSolver_;

  VectorXd lb_, ub_;

  tuple<VectorXd, double, int> solver();

  VectorXd x_;
  int getXDim() { return (int)x_.rows(); }

  // Dbg constructs
  string md_ = "Optimization/optimizers/sqp_snopt";
  string cl_ = "SQP_SNOPT";
  string im_ = "", wm_ = "", em_ = "";
  Settings::VerbParams vp_;
  QString sqp_snopt_log_;
};

}
}

#endif //FIELDOPT_OPTIMIZATION_OPTIMIZERS_SQP_SNOPT_SQP_SNOPT_H_
