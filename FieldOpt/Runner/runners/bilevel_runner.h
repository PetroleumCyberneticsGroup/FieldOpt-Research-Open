/***********************************************************
Copyright (C) 2021
Mathias Bellout <chakibbb-pcg@gmail.com>

bellout - Sat Aug 28 2021 14:16:50 week 34 CET+0200

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

#ifndef FIELDOPT_RUNNER_RUNNERS_BILEVEL_RUNNER_H_
#define FIELDOPT_RUNNER_RUNNERS_BILEVEL_RUNNER_H_

#include "abstract_runner.h"

namespace Runner {

class MainRunner;

class BilevelRunner : public AbstractRunner
{
  friend class MainRunner;
 private:
  explicit BilevelRunner(RuntimeSettings *runtime_settings);

  string md_ = "Runner";
  string cl_ = "BilevelRunner";

 private:
  void Execute() override;

  void prntDbg(int lvl,
               Optimization::Optimizer *optmzr,
               Optimization::Case *cs);
};

}
#endif //FIELDOPT_RUNNER_RUNNERS_BILEVEL_RUNNER_H_
