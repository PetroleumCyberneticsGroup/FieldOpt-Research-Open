/***********************************************************
Copyright (C) 2020-
Thiago L. Silva <thiagolims@gmail.com>
Mathias Bellout <chakibbb-pcg@gmail.com>

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

#ifndef FIELDOPT_DRILLING_RUNNER_H
#define FIELDOPT_DRILLING_RUNNER_H

#include "abstract_runner.h"
#include "Model/model.h"
#include "Model/drilling/drilling.h"
#include "Model/drilling/drilling_schedule.h"

namespace Runner {

class MainRunner;

class DrillingRunner : public AbstractRunner {

 friend class MainRunner;

 private:
  Model::Drilling::Drilling* drilling_;

  void InitializeDrillingWorkflow();

 // AbstractRunner interface
 public:
  explicit DrillingRunner(RuntimeSettings *runtime_settings);
  void Execute() override;




};

}
#endif //FIELDOPT_DRILLING_RUNNER_H
