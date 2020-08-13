/***********************************************************
Created: 16.12.2015 2015 by einar

Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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
#include "Optimization/objective/weightedsum.h"
#include "Simulation/simulator_interfaces/eclsimulator.h"
#include "Simulation/simulator_interfaces/adgprssimulator.h"
#include "Utilities/math.hpp"
#include "Utilities/printer.hpp"
#include "Utilities/verbosity.h"

namespace Runner {

using Printer::info;
using Printer::ext_info;
using Printer::num2str;

AbstractRunner::AbstractRunner(RuntimeSettings *runtime_settings) {
  runtime_settings_ = runtime_settings;

  settings_ = nullptr;
  model_ = nullptr;
  simulator_ = nullptr;
  objective_function_ = nullptr;
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

void AbstractRunner::InitializeSettings(QString output_subdirectory) {
  QString output_directory = QString::fromStdString(runtime_settings_->paths().GetPath(Paths::OUTPUT_DIR));
  if (output_subdirectory.length() > 0) {
    output_directory.append(QString("/%1/").arg(output_subdirectory));
  }
  runtime_settings_->paths().SetPath(Paths::OUTPUT_DIR, output_directory.toStdString());
  settings_ = new Settings::Settings(runtime_settings_->paths());
  vp_ = settings_->global()->verbParams();
  // settings_->global()->showVerbParams();

  if (!DirExists(output_directory, vp_)) {
    CreateDir(output_directory, vp_);
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
    settings_->paths().SetPath(Paths::GRID_FILE,ensemble_helper_.GetBaseRealization().grid());
  }

  model_ = new Model::Model(*settings_, logger_);
}

void AbstractRunner::InitializeSimulator() {
  if (model_ == nullptr)
    throw std::runtime_error("Model must be initialized before Simulator.");

  switch (settings_->simulator()->type()) {
    case ::Settings::Simulator::SimulatorType::ECLIPSE: {
      string tm = "Using ECLIPSE simulator.";
      if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::ECLSimulator(settings_, model_);
      break;
    }

    case ::Settings::Simulator::SimulatorType::ADGPRS: {
      auto tm = "Using ADGPRS simulator.";
      if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::AdgprsSimulator(settings_, model_);
      break;
    }

    case ::Settings::Simulator::SimulatorType::Flow: {
      auto tm = "Using Flow simulator.";
      if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::ECLSimulator(settings_, model_);
      break;
    }

    case ::Settings::Simulator::SimulatorType::INTERSECT: {
      auto tm = "Using INTERSECT simulator.";
      if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::IXSimulator(settings_, model_);
      break;
    }

    default:
      throw std::runtime_error("Unable to initialize runner: simulator "
                               "type set in JSON driver file not recognized.");
  }
}

void AbstractRunner::EvaluateBaseModel() {
  if (simulator_ == nullptr)
    throw std::runtime_error("Simulator must be initialized before evaluating base model.");

  if (is_ensemble_run_) {
    auto tm = "Simulating ensemble base case.";
    if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
    simulator_->Evaluate(ensemble_helper_.GetBaseRealization(), 10000, 4);
  }
  else if (!simulator_->results()->isAvailable()) {
    auto tm = "Simulating base case.";
    if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
    simulator_->Evaluate();
  }
}

void AbstractRunner::InitializeObjectiveFunction() {
  if (simulator_ == nullptr || settings_ == nullptr)
    throw std::runtime_error("Simulator and the Settings must be initialized before ObjectiveFunction.");

  switch (settings_->optimizer()->objective().type) {
    case Settings::Optimizer::ObjectiveType::WeightedSum: {
      auto tm = "Using WeightedSum-type objective function.";
      if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      objective_function_ = new Optimization::Objective::WeightedSum(
        settings_->optimizer(), simulator_->results(), model_);
      break;
    }

    case Settings::Optimizer::ObjectiveType::NPV: {
      auto tm = "Using NPV-type objective function.";
      if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      objective_function_ = new Optimization::Objective::NPV(
        settings_->optimizer(), simulator_->results(), model_);
      break;
    }

    default:
      throw std::runtime_error("Unable to initialize Runner: "
                               "ObjectiveFunction type not recognized.");
  }
}

void AbstractRunner::InitializeBaseCase() {
  if (objective_function_ == nullptr || model_ == nullptr)
    throw std::runtime_error("ObjectiveFunction and Model must be initialized before BaseCase.");

  base_case_ = new Optimization::Case(model_->variables()->GetBinaryVariableValues(),
                                      model_->variables()->GetDiscreteVariableValues(),
                                      model_->variables()->GetContinousVariableValues());

  if (!simulator_->results()->isAvailable()) {
    auto tm = "Simulation results are unavailable. Setting base case objective function value to sentinel value.";
    if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm, vp_.lnw); }
    base_case_->set_objective_function_value(sentinelValue());
  } else {
    model_->wellCost(settings_->optimizer());
    base_case_->set_objective_function_value(objective_function_->value());
  }
  auto tm = "Base case objective function value set to " + num2str(base_case_->objective_function_value());
  if (VERB_RUN >= 1 || vp_.vRUN >= 1) { info(tm,  vp_.lnw); }
}

void AbstractRunner::InitializeOptimizer() {
  if (base_case_ == nullptr || model_ == nullptr)
    throw std::runtime_error("The Base Case and the Model must be initialized before the Optimizer");

  switch (settings_->optimizer()->type()) {

    case Settings::Optimizer::OptimizerType::Compass: {
      if (VERB_RUN >= 1) info("Using CompassSearch optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::CompassSearch(settings_->optimizer(),
                                                               base_case_,
                                                               model_->variables(),
                                                               model_->grid(),
                                                               logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::APPS: {
      if (VERB_RUN >= 1) info("Using APPS optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::APPS(settings_->optimizer(),
                                                      base_case_,
                                                      model_->variables(),
                                                      model_->grid(),
                                                      logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::GeneticAlgorithm: {
      if (VERB_RUN >= 1) info("Using GeneticAlgorithm optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::RGARDD(settings_->optimizer(),
                                                        base_case_,
                                                        model_->variables(),
                                                        model_->grid(),
                                                        logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::EGO: {
      if (VERB_RUN >= 1) info("Using EGO optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::BayesianOptimization::EGO(settings_->optimizer(),
                                                                           base_case_,
                                                                           model_->variables(),
                                                                           model_->grid(),
                                                                           logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::ExhaustiveSearch2DVert: {
      if (VERB_RUN >= 1) info("Using ExhaustiveSearch2DVert optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::ExhaustiveSearch2DVert(settings_->optimizer(),
                                                                        base_case_,
                                                                        model_->variables(),
                                                                        model_->grid(),
                                                                        logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::Hybrid: {
      if (VERB_RUN >= 1) info("Using Hybrid optimization algorithm.", vp_.lnw);
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
      optimizer_ = new Optimization::Optimizers::TrustRegionOptimization(settings_->optimizer(),
                                                                         base_case_,
                                                                         model_->variables(),
                                                                         model_->grid(),
                                                                         logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::PSO: {
      if (VERB_RUN >= 1) info("Using PSO optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::PSO(settings_->optimizer(),
                                                     base_case_,
                                                     model_->variables(),
                                                     model_->grid(),
                                                     logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::CMA_ES: {
      if (VERB_RUN >= 1) info("Using CMA_ES optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::CMA_ES(settings_->optimizer(),
                                                        base_case_,
                                                        model_->variables(),
                                                        model_->grid(),
                                                        logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::VFSA: {
      if (VERB_RUN >= 1) info("Using VFSA optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::VFSA(settings_->optimizer(),
                                                      base_case_,
                                                      model_->variables(),
                                                      model_->grid(),
                                                      logger_
      );
      break;
    }
    case Settings::Optimizer::OptimizerType::SPSA: {
      if (VERB_RUN >= 1) info("Using SPSA optimization algorithm.", vp_.lnw);
      optimizer_ = new Optimization::Optimizers::SPSA(settings_->optimizer(),
                                                      base_case_,
                                                      model_->variables(),
                                                      model_->grid(),
                                                      logger_
      );
      break;
    }
    default:
      throw std::runtime_error("Unable to initialize runner: optimization algorithm set in driver file not recognized.");
  }
  optimizer_->EnableConstraintLogging(QString::fromStdString(runtime_settings_->paths().GetPath(Paths::OUTPUT_DIR)));
}

void AbstractRunner::InitializeBookkeeper() {
  if (settings_ == nullptr || optimizer_ == nullptr)
    throw std::runtime_error("The Settings and the Optimizer must be initialized before the Bookkeeper.");

  bookkeeper_ = new Bookkeeper(settings_, optimizer_->case_handler());
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
  for (auto var : optimizer_->GetTentativeBestCase()->integer_variables().keys()) {
    auto prop_name = model_->variables()->GetDiscreteVariable(var)->name();
    auto prop_val = optimizer_->GetTentativeBestCase()->integer_variables()[var];
    std::cout << "\t" << prop_name.toStdString() << "\t" << prop_val << std::endl;
  }
  for (auto var : optimizer_->GetTentativeBestCase()->real_variables().keys()) {
    auto prop_name = model_->variables()->GetContinousVariable(var)->name();
    auto prop_val = optimizer_->GetTentativeBestCase()->real_variables()[var];
    std::cout << "\t" << prop_name.toStdString() << "\t" << prop_val << std::endl;
  }
  for (auto var : optimizer_->GetTentativeBestCase()->binary_variables().keys()) {
    auto prop_name = model_->variables()->GetBinaryVariable(var)->name();
    auto prop_val = optimizer_->GetTentativeBestCase()->binary_variables()[var];
    std::cout << "\t" << prop_name.toStdString() << "\t" << prop_val << std::endl;
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
