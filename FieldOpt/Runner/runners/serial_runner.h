/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2021 Mathias Bellout
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

#ifndef SERIALRUNNER_H
#define SERIALRUNNER_H

#include "abstract_runner.h"

namespace Runner {

using TC = Optimization::Optimizer::TerminationCondition;
using ES = Optimization::Case::CaseState::EvalStatus;
using EM = Optimization::Case::CaseState::ErrorMessage;
using OT = Settings::Optimizer::ObjectiveType;

class MainRunner;

/*!
 * \brief The SerialRunner class is a simple implementation
 * of a serial runner. It interacts with the optimizer and
 * starts simulations in series (i.e. it is _very_ slow).
 */
class SerialRunner : public AbstractRunner
{
  friend class MainRunner;
 private:
  explicit SerialRunner(RuntimeSettings *runtime_settings);
  explicit SerialRunner(RuntimeSettings *runtime_settings,
                        Optimization::Case* base_case,
                        Model::ModelSynchronizationObject* mso);

  string md_ = "Runner";
  string cl_ = "SerialRunner";

  // AbstractRunner interface
 private:
  void Execute() override;
};

}

#endif // SERIALRUNNER_H
