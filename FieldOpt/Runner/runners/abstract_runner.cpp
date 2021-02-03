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

#include <Simulation/simulator_interfaces/flowsimulator.h>
#include <Optimization/optimizers/APPS.h>
#include <Optimization/optimizers/GeneticAlgorithm.h>
#include <Optimization/optimizers/RGARDD.h>
#include <Optimization/hybrid_optimizer.h>
#include <Optimization/optimizers/bayesian_optimization/EGO.h>
#include <Optimization/optimizers/trust_region/TrustRegionOptimization.h>
#include "Optimization/optimizers/PSO.h"
#include "Optimization/optimizers/CMA_ES.h"
#include "Optimization/optimizers/VFSA.h"
#include "Optimization/optimizers/SPSA.h"
#include "Simulation/simulator_interfaces/ix_simulator.h"
#include "abstract_runner.h"
#include "Optimization/optimizers/compass_search.h"
#include "Optimization/optimizers/ExhaustiveSearch2DVert.h"
#include "Simulation/simulator_interfaces/eclsimulator.h"
#include "Simulation/simulator_interfaces/adgprssimulator.h"
#include "Utilities/math.hpp"
#include "Utilities/printer.hpp"
#include "Utilities/verbosity.h"

namespace Runner {

using Printer::info;
using Printer::ext_info;
using Printer::num2str;
using std::runtime_error;

AbstractRunner::AbstractRunner(RuntimeSettings *runtime_settings) {
  runtime_settings_ = runtime_settings;

  settings_ = nullptr;
  model_ = nullptr;
  simulator_ = nullptr;
  objf_ = nullptr;
  base_case_ = nullptr;
  optimizer_ = nullptr;
  bookkeeper_ = nullptr;

  logger_ = nullptr;
  is_ensemble_run_ = false;
}

double AbstractRunner::sentinelValue() const {
  if (settings_->optimizer()->mode() == Settings::Optimizer::OptimizerMode::Minimize)
    return -1*sentinel_value_;
  return sentinel_value_;
}

void AbstractRunner::InitializeSettings(QString output_subdir) {
  QString output_dir = runtime_settings_->paths().GetPathQstr(Paths::OUTPUT_DIR);
  if (output_subdir.length() > 0) {
    output_dir.append(QString("/%1/").arg(output_subdir));
  }
  runtime_settings_->paths().SetPath(Paths::OUTPUT_DIR, output_dir.toStdString());
  settings_ = new Settings::Settings(runtime_settings_->paths());
  vp_ = settings_->global()->verbParams();
  // settings_->global()->showVerbParams();

  if (!DirExists(output_dir, vp_)) {
    CreateDir(output_dir, vp_);
  }

  if (settings_->simulator()->is_ensemble()) {
    is_ensemble_run_ = true;
    ensemble_helper_ = EnsembleHelper(settings_->simulator()->get_ensemble(),
                                      settings_->optimizer()->parameters().rng_seed);
  } else {
    is_ensemble_run_ = false;
  }
}

void AbstractRunner::InitializeModel() {
  if (settings_ == nullptr)
    throw std::runtime_error("Settings must be initialized before Model.");

  if (is_ensemble_run_) {
    settings_->paths().SetPath(Paths::GRID_FILE,
                               ensemble_helper_.GetBaseRealization().grid());
  }

  model_ = new Model::Model(*settings_, logger_);
}

void AbstractRunner::InitializeSimulator() {
  if (model_ == nullptr)
    throw std::runtime_error("Model must be initialized before Simulator.");

  switch (settings_->simulator()->type()) {
    case ::Settings::Simulator::SimulatorType::ECLIPSE: {
      string tm = "Using ECLIPSE simulator.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::ECLSimulator(settings_,
                                                model_);
      break;
    }

    case ::Settings::Simulator::SimulatorType::ADGPRS: {
      auto tm = "Using ADGPRS simulator.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::AdgprsSimulator(settings_,
                                                   model_);
      break;
    }

    case ::Settings::Simulator::SimulatorType::Flow: {
      auto tm = "Using Flow simulator.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::ECLSimulator(settings_,
                                                model_);
      break;
    }

    case ::Settings::Simulator::SimulatorType::INTERSECT: {
      auto tm = "Using INTERSECT simulator.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::IXSimulator(settings_,
                                               model_);
      break;
    }

    default: {
      string em = "Unable to initialize runner: simulator ";
      em += "type set in JSON driver file not recognized.";
      throw runtime_error(em);
    }
  }
}

void AbstractRunner::EvaluateBaseModel() {
  if (simulator_ == nullptr) {
    string em = "Simulator must be initialized before evaluating base model.";
    throw runtime_error(em);
  }

  if (is_ensemble_run_) {
    auto tm = "Simulating ensemble base case.";
    if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }

    auto base_rlz = ensemble_helper_.GetBaseRealization();
    simulator_->Evaluate(base_rlz, 10000, 4);

  } else if (!simulator_->results()->isAvailable()) {
    auto tm = "Simulating base case.";
    if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }

    simulator_->Evaluate();
  }
}

void AbstractRunner::InitializeObjectiveFunction() {
  if (simulator_ == nullptr || settings_ == nullptr) {
    string em = "Simulator and Settings must be initialized before ObjectiveFunction.";
    throw runtime_error(em);
  }

  switch (settings_->optimizer()->objective().type) {
    case Settings::Optimizer::ObjectiveType::WeightedSum: {
      auto tm = "Using WeightedSum-type objective function.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      objf_ =
        new Optimization::Objective::WeightedSum(settings_->optimizer(),
                                                 simulator_->results(),
                                                 model_);
      break;
    }

    case Settings::Optimizer::ObjectiveType::NPV: {
      auto tm = "Using NPV-type objective function.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      objf_ =
        new Optimization::Objective::NPV(settings_->optimizer(),
                                         simulator_->results(),
                                         model_);
      break;
    }

    case Settings::Optimizer::ObjectiveType::Augmented: {
      auto tm = "Using Augmented objective.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      objf_ =
        new Optimization::Objective::Augmented(settings_->optimizer(),
                                               simulator_->results(),
                                               model_);
      break;
    }

    default: {
      string em = "Unable to initialize Runner: ObjectiveFunction type not recognized.";
      throw runtime_error(em);
    }
  }
}

void AbstractRunner::InitializeBaseCase() {
  if (objf_ == nullptr || model_ == nullptr) {
    string em = "ObjectiveFunction and Model must be initialized before BaseCase.";
    throw runtime_error(em);
  }
  base_case_ = new Optimization::Case(model_->variables()->GetBinVarValues(),
                                      model_->variables()->GetDiscVarValues(),
                                      model_->variables()->GetContVarValues());

  if (!simulator_->results()->isAvailable()) {
    base_case_->set_objf_value(sentinelValue());

    string tm = "Simulation results are unavailable.";
    tm += "Base case objective set to sentinel value.";
    if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }

  } else {
    model_->wellCost(settings_->optimizer());
    if (settings_->optimizer()->objective().type == Settings::Optimizer::ObjectiveType::Augmented) {
      base_case_->set_objf_value(objf_->value(true));
    } else {
      base_case_->set_objf_value(objf_->value());
    }

    string tm = "Base case objective function value set to ";
    tm += num2str(base_case_->objf_value(), 8, 1);
    if (vp_.vRUN >= 1) { info(tm,  vp_.lnw); }
  }
}

namespace Optzr = Optimization::Optimizers;

void AbstractRunner::InitializeOptimizer() {
  if (base_case_ == nullptr || model_ == nullptr) {
    string em = "Base Case and Model must be initialized before the Optimizer";
    throw runtime_error(em);
  }

  switch (settings_->optimizer()->type()) {

    case Settings::Optimizer::OptimizerType::Compass: {
      if (vp_.vRUN >= 1) info("Using CompassSearch optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::CompassSearch(settings_->optimizer(),
                                            base_case_,
                                            model_->variables(),
                                            model_->grid(),
                                            logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::APPS: {
      if (vp_.vRUN >= 1) info("Using APPS optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::APPS(settings_->optimizer(),
                                   base_case_,
                                   model_->variables(),
                                   model_->grid(),
                                   logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::GeneticAlgorithm: {
      if (vp_.vRUN >= 1) info("Using GeneticAlgorithm optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::RGARDD(settings_->optimizer(),
                                     base_case_,
                                     model_->variables(),
                                     model_->grid(),
                                     logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::EGO: {
      if (vp_.vRUN >= 1) info("Using EGO optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::BayesianOptimization::EGO(settings_->optimizer(),
                                                        base_case_,
                                                        model_->variables(),
                                                        model_->grid(),
                                                        logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::ExhaustiveSearch2DVert: {
      if (vp_.vRUN >= 1) info("Using ExhaustiveSearch2DVert optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::ExhaustiveSearch2DVert(settings_->optimizer(),
                                                     base_case_,
                                                     model_->variables(),
                                                     model_->grid(),
                                                     logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::Hybrid: {
      if (vp_.vRUN >= 1) info("Using Hybrid optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::HybridOptimizer(settings_->optimizer(),
                                                     base_case_,
                                                     model_->variables(),
                                                     model_->grid(),
                                                     logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::TrustRegionOptimization: {
      if (VERB_RUN >= 1) info("Using Trust Region optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::TrustRegionOptimization(settings_->optimizer(),
                                                      base_case_,
                                                      model_->variables(),
                                                      model_->grid(),
                                                      logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::PSO: {
      if (VERB_RUN >= 1) info("Using PSO optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::PSO(settings_->optimizer(),
                                  base_case_,
                                  model_->variables(),
                                  model_->grid(),
                                  logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::CMA_ES: {
      if (VERB_RUN >= 1) info("Using CMA_ES optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::CMA_ES(settings_->optimizer(),
                                     base_case_,
                                     model_->variables(),
                                     model_->grid(),
                                     logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::VFSA: {
      if (VERB_RUN >= 1) info("Using VFSA optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::VFSA(settings_->optimizer(),
                                   base_case_,
                                   model_->variables(),
                                   model_->grid(),
                                   logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::SPSA: {
      if (VERB_RUN >= 1) info("Using SPSA optimization algorithm.", vp_.lnw);
      optimizer_ = new Optzr::SPSA(settings_->optimizer(),
                                   base_case_,
                                   model_->variables(),
                                   model_->grid(),
                                   logger_
      );
      break;
    }
    default: {
      string em = "Unable to initialize runner: ";
      em += "Optimization algorithm set in driver file not recognized.";
      throw runtime_error(em);
    }

  }
  optimizer_->EnableConstraintLogging(
    runtime_settings_->paths().GetPathQstr(Paths::OUTPUT_DIR));
}

void AbstractRunner::InitializeBookkeeper() {
  if (settings_ == nullptr || optimizer_ == nullptr) {
    string em = "The Settings and Optimizer must be initialized before the Bookkeeper.";
    throw runtime_error(em);
  }

  bookkeeper_ = new Bookkeeper(settings_,
                               optimizer_->case_handler());
}

void AbstractRunner::InitializeLogger(QString output_subdir, bool write_logs) {
  logger_ = new Logger(runtime_settings_, output_subdir, write_logs);
}

void AbstractRunner::PrintCompletionMessage() const {

  std::cout << "Optimization complete: ";
  switch (optimizer_->IsFinished()) {
    case Optimization::Optimizer::TerminationCondition::MAX_EVALS_REACHED:
      std::cout << "maximum number of evaluations reached (not converged)." << std::endl;
      break;
    case Optimization::Optimizer::TerminationCondition::MINIMUM_STEP_LENGTH_REACHED:
      std::cout << "minimum step length reached (converged)." << std::endl;
      break;
    default: std::cout << "Unknown termination reason." << std::endl;
  }

  std::cout << "Best case at termination:" << optimizer_->GetTentativeBestCase()->id().toString().toStdString() << std::endl;
  std::cout << "Variable values: " << std::endl;

  for (auto var : optimizer_->GetTentativeBestCase()->integer_variables()) {
    auto prop_name = model_->variables()->GetDiscreteVariable(var.first)->name();
    std::cout << "\t" << prop_name.toStdString() << "\t" << var.second << std::endl;
  }

  for (auto var : optimizer_->GetTentativeBestCase()->real_variables()) {
    auto prop_name = model_->variables()->GetContinuousVariable(var.first)->name();
    std::cout << "\t" << prop_name.toStdString() << "\t" << var.second << std::endl;
  }

  for (auto var : optimizer_->GetTentativeBestCase()->binary_variables()) {
    auto prop_name = model_->variables()->GetBinaryVariable(var.first)->name();
    std::cout << "\t" << prop_name.toStdString() << "\t" << var.second << std::endl;
  }
}

int AbstractRunner::timeoutValue() const {
  if (simulation_times_.size() == 0 || runtime_settings_->simulation_timeout() == 0)
    return 10000;
  else {
    return calc_median(simulation_times_) * runtime_settings_->simulation_timeout();
  }
}

void AbstractRunner::FinalizeInitialization(bool write_logs) {
  if (write_logs) {
    logger_->AddEntry(runtime_settings_);
    logger_->FinalizePrerunSummary();
  }
}

void AbstractRunner::FinalizeRun(bool write_logs) {
  if (optimizer_ != nullptr) { // This indicates whether or not we're on a worker process
    model_->ApplyCase(optimizer_->GetTentativeBestCase());
    simulator_->WriteDriverFilesOnly();
    PrintCompletionMessage();
  }
  model_->Finalize();
  if (write_logs)
    logger_->FinalizePostrunSummary();
}

}
