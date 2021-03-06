/***********************************************************
Created: 16.12.2015 2015 by einar

Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#ifndef ABSTRACTRUNNER_H
#define ABSTRACTRUNNER_H

#include "runtime_settings.h"
#include "Model/model.h"
#include "Optimization/optimizer.h"
#include "Optimization/case.h"
#include "Simulation/simulator_interfaces/simulator.h"
#include "Settings/settings.h"
#include "bookkeeper.h"
#include "Runner/logger.h"
#include "ensemble_helper.h"
#include <vector>

#include "Optimization/objective/objective.h"
#include "Optimization/objective/NPV.h"
#include "Optimization/objective/weightedsum.h"
#include "Optimization/objective/augmented.h"

namespace Runner {

using Printer::info;
using Printer::ext_info;
using Printer::num2str;
using Printer::E;
using std::runtime_error;

using TC = Optimization::Optimizer::TerminationCondition;

class MainRunner;

/*!
 * \brief The AbstractRunner class is the abstract parent
 * class for all runners. It should only be constructed by
 * the MainRunner class.
 *
 * This class initializes the primary objects needed and
 * provides some utility functions for logging.
 *
 * It also defines the purely virtual Execute() method which
 * should be implemented by all concrete runners.
 *
 * todo: Create a method to get the timeout it seconds that
 * uses the recorded simulation times and the timeout argument.
 */
class AbstractRunner
{
  friend class MainRunner;
 private:

  /*!
   * \brief Execute starts the actual optimization run
   * and should not return until the optimization is done.
   */
  virtual void Execute() = 0;

  //!< Value to be used as a sentinel value for the
  //!< objective function of cases that cannot be evaluated.
  const double sentinel_value_ = 0.0001;

  string md_ = "Runner";
  string cl_ = "AbstractRunner";

 protected:
  explicit AbstractRunner(RuntimeSettings *runtime_settings);

  Bookkeeper *bookkeeper_;
  Model::Model *model_;
  Settings::Settings *settings_;
  RuntimeSettings *rts_;
  Optimization::Case *base_case_;
  Optimization::Case *optz_case_;

  Optimization::Optimizer *optimizer_;
  Optimization::Objective::Objective *objf_;
  Simulation::Simulator *simulator_;
  Logger *logger_;
  std::vector<int> sim_times_;
  bool is_ensemble_run_;
  EnsembleHelper ensemble_helper_;

  // void E(string m) const {
  //   m = "[mod: " + md_ + "] [cls: " + cl_ + "] " + m;
  //   throw runtime_error(m);
  // };

  string im_ = "", wm_ = "", em_ = "";
  Settings::VerbParams vp_;

  void PrintCompletionMessage();

  void ComputeOptmzdCase();

  /*!
   * \brief sentinelValue Get the sentinel value to be used
   * as objective function values for cases that fail.
   *
   * When maximizing, this will be 0.0001; when minimizing,
   * this will be -0.0001.
   * \return
   */
  double sentinelValue() const;

  /*!
   * @brief Get timeout value to be used when starting
   * simulations. It is calculated from the recorded
   * (successful) simulation times and the timeout value
   * provided as an argument when running the program.
   *
   * If there either have not been any recorded simulation
   * times or the timeout argument was not provided, 10,000
   * will be returned.
   * @return
   */
  int timeoutVal() const;

  void InitializeSettings(const QString& output_subdirectory="");
  void InitializeModel();
  void InitializeSimulator();
  void EvaluateBaseModel();
  void InitializeObjectiveFunction();
  void InitializeBaseCase();
  void InitializeOptimizer();
  void InitializeBookkeeper();

  //!< Write the pre-run summary
  void FinalizeInitialization(bool write_logs);

  //!< Finalize the run, writing data to the summary log.
  void FinalizeRun(bool write_logs);

  /*!
   * @brief Initialize the logger.
   * @param output_subdir Optional subdir in
   * the output dir to write the logs in.
   */
  void InitializeLogger(QString output_subdir="", bool write_logs=true);
};

}

#endif // ABSTRACTRUNNER_H
