/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020- Mathias Bellout
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

#ifndef RUNNER_H
#define RUNNER_H

#include <stdexcept>
#include "runtime_settings.h"
#include "abstract_runner.h"

namespace Runner {

/*!
 * \brief MainRunner class takes care of initializing and
 * starting the actual runner indicated in runtime settings.
 */
class MainRunner
{
 public:
  explicit MainRunner(RuntimeSettings *rts);
  explicit MainRunner(RuntimeSettings *rts,
                      Optimization::Case* base_case,
                      Model::ModelSynchronizationObject* mso);

  /*!
   * \brief Execute Start the optimization run by
   * calling the Execute function in the simulator.
   */
  void Execute();

  Model::Model* getModel() { return runner_->getModel(); }
  Settings::Settings* getSettings() { return runner_->getSettings(); }
  Optimization::Optimizer* getOptimizer() { return runner_->getOptimizer(); }
  Simulation::Simulator* getSimulator() { return runner_->getSimulator(); }

  /*!
   * @brief Initializes runner modules.
   */
   void InitializeModules() { runner_->InitializeModules(); }

   /*!
    * @brief Replace the optimizer.
    * @param opt New optimizer setttings.
    */
    void ReplaceOptimizer(Settings::Optimizer *opt) {
      runner_->ReplaceOptimizer(opt);
    }

   /*!
    * @brief Evaluates the base case.
    * @return Base case with objective value.
    */
    Optimization::Case* evaluateBaseCase();
  
 private:
  RuntimeSettings *runtime_settings_;
  AbstractRunner *runner_;

  // void E(string m) const {
  //   m = "[mod: " + md_ + "] [cls: " + cl_ + "] " + m;
  //   throw runtime_error(m);
  // };

  string md_ = "Runner";
  string cl_ = "MainRunner";
};

}

#endif // RUNNER_H
