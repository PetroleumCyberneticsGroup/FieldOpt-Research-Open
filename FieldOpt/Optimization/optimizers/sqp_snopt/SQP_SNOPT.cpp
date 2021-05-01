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

#include "SQP_SNOPT.h"

namespace Optimization {
namespace Optimizers {

SQP_SNOPT::SQP_SNOPT(Settings::Optimizer *seto, Case *bscs,
  VarPropContainer *vars, Reservoir::Grid::Grid *grid,
  Logger *logger, CaseHandler *case_handler,
  ConstraintHandler *constraint_handler)
  : Optimizer(seto, bscs, vars, grid, logger,
    case_handler, constraint_handler) {

  vars_ = vars;
  bscs_ = bscs;

  if (enable_logging_) {
    logger_->AddEntry(this);
  }

}

}
}
