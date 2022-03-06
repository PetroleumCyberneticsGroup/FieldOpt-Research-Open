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
#include <APPS.h>
#include <PSO.h>
#include <trust_region/TrustRegionOptimization.h>
#include <dftr/DFTR.h>
#include <sqp_snopt/SQP_SNOPT.h>
#include <GeneticAlgorithm.h>
#include <RGARDD.h>
#include <bayesian_optimization/EGO.h>
#include <ExhaustiveSearch2DVert.h>
#include <compass_search.h>
#include <CMA_ES.h>
#include <VFSA.h>
#include <SPSA.h>

#include <Optimization/hybrid_optimizer.h>
#include "Simulation/simulator_interfaces/ix_simulator.h"
#include "abstract_runner.h"
#include "Simulation/simulator_interfaces/eclsimulator.h"
#include "Simulation/simulator_interfaces/adgprssimulator.h"
#include "Utilities/math.hpp"
#include "Utilities/printer.hpp"

namespace Runner {

using OptzrMod = Settings::Optimizer::OptimizerMode;
using OptzrTyp = Settings::Optimizer::OptimizerType;

using SimTyp = Settings::Simulator::SimulatorType;
using ObjTyp = Settings::Optimizer::ObjectiveType;

using Optimization::Objective::WeightedSum;
using Optimization::Objective::NPV;
using Optimization::Objective::Augmented;

namespace Optzr = Optimization::Optimizers;

// ┌─┐  ┌┐   ┌─┐  ┌┬┐  ┬─┐  ┌─┐  ┌─┐  ┌┬┐    ╦═╗  ╦ ╦  ╔╗╔  ╔╗╔  ╔═╗  ╦═╗
// ├─┤  ├┴┐  └─┐   │   ├┬┘  ├─┤  │     │     ╠╦╝  ║ ║  ║║║  ║║║  ║╣   ╠╦╝
// ┴ ┴  └─┘  └─┘   ┴   ┴└─  ┴ ┴  └─┘   ┴     ╩╚═  ╚═╝  ╝╚╝  ╝╚╝  ╚═╝  ╩╚═
AbstractRunner::AbstractRunner(RuntimeSettings *runtime_settings) {
  rts_ = runtime_settings;

  settings_ = nullptr;
  model_ = nullptr;
  simulator_ = nullptr;
  objf_ = nullptr;
  base_case_ = nullptr;
  optmzr_ = nullptr;
  bookkeeper_ = nullptr;

  logger_ = nullptr;
  is_ensemble_run_ = false;
}

// ┌─┐  ┌─┐  ┌┐┌  ┌┬┐  ┬  ┌┐┌  ┌─┐  ┬    ╦  ╦  ╔═╗  ╦
// └─┐  ├┤   │││   │   │  │││  ├┤   │    ╚╗╔╝  ╠═╣  ║
// └─┘  └─┘  ┘└┘   ┴   ┴  ┘└┘  └─┘  ┴─┘   ╚╝   ╩ ╩  ╩═╝
double AbstractRunner::sentinelValue() const {
  if (seto_->mode() == OptzrMod::Minimize) {
    return -1*sentinel_value_;
  }
  return sentinel_value_;
}

// ┬  ┌┐┌  ┬  ┌┬┐  ╔═╗  ╔═╗  ╔╦╗  ╔╦╗  ╦  ╔╗╔  ╔═╗  ╔═╗
// │  │││  │   │   ╚═╗  ║╣    ║    ║   ║  ║║║  ║ ╦  ╚═╗
// ┴  ┘└┘  ┴   ┴   ╚═╝  ╚═╝   ╩    ╩   ╩  ╝╚╝  ╚═╝  ╚═╝
void AbstractRunner::InitSettings(const QString& output_subdir) {

  // Output dir
  QString output_dir = rts_->paths().GetPathQstr(Paths::OUTPUT_DIR);
  if (output_subdir.length() > 0) {
    output_dir.append(QString("/%1/").arg(output_subdir));
  }
  if (!DirExists(output_dir, vp_)) { CreateDir(output_dir, vp_); }
  rts_->paths().SetPath(Paths::OUTPUT_DIR, output_dir.toStdString());

  // Optmzd dir
  string optz_dir = output_dir.toStdString();
  if(rts_->runner_type() == RuntimeSettings::SERIAL) {
    optz_dir += "/optcs";
  } else if(rts_->runner_type() == RuntimeSettings::MPISYNC) {
    optz_dir += "../optcs";
  }

  if (!DirExists(optz_dir, vp_)) { CreateDir(optz_dir, vp_); }
  rts_->paths().SetPath(Paths::OPTMZD_DIR, optz_dir);

  // Settings
  settings_ = new Settings::Settings(rts_->paths());
  vp_ = settings_->global()->verbParams();
  // settings_->global()->showVerbParams();
  seto_ = settings_->optimizer();

  if (settings_->simulator()->is_ensemble()) {
    is_ensemble_run_ = true;
    ensemble_helper_ = EnsembleHelper(
      settings_->simulator()->get_ensemble(),
      seto_->parameters().rng_seed);
  } else {
    is_ensemble_run_ = false;
  }
}

// ┬  ┌┐┌  ┬  ┌┬┐  ╔╦╗  ╔═╗  ╔╦╗  ╔═╗  ╦
// │  │││  │   │   ║║║  ║ ║   ║║  ║╣   ║
// ┴  ┘└┘  ┴   ┴   ╩ ╩  ╚═╝  ═╩╝  ╚═╝  ╩═╝
void AbstractRunner::InitModel() {
  if (settings_ == nullptr) { E("Settings must be initialized b/f Model.", md_, cl_); }

  if (is_ensemble_run_) {
    settings_->paths().SetPath(Paths::GRID_FILE,
                               ensemble_helper_.GetBaseRealization().grid());
  }

  model_ = new Model::Model(*settings_, logger_);
  vars_ = model_->variables();
}

// ┬  ┌┐┌  ┬  ┌┬┐  ╔═╗  ╦  ╔╦╗  ╦ ╦  ╦    ╔═╗  ╔╦╗  ╔═╗  ╦═╗
// │  │││  │   │   ╚═╗  ║  ║║║  ║ ║  ║    ╠═╣   ║   ║ ║  ╠╦╝
// ┴  ┘└┘  ┴   ┴   ╚═╝  ╩  ╩ ╩  ╚═╝  ╩═╝  ╩ ╩   ╩   ╚═╝  ╩╚═
void AbstractRunner::InitSimulator() {
  if (model_ == nullptr) { E("Model must be initialized b/f Simulator.", md_, cl_); }

  switch (settings_->simulator()->type()) {

    case SimTyp::ECLIPSE: {
      string tm = "Using ECLIPSE simulator.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::ECLSimulator(settings_, model_);
      break;
    }

    case SimTyp::ADGPRS: {
      auto tm = "Using ADGPRS simulator.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::AdgprsSimulator(settings_, model_);
      break;
    }

    case SimTyp::Flow: {
      auto tm = "Using Flow simulator.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::ECLSimulator(settings_, model_);
      break;
    }

    case SimTyp::INTERSECT: {
      auto tm = "Using INTERSECT simulator.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      simulator_ = new Simulation::IXSimulator(settings_, model_);
      break;
    }

    default: { E("Simulator not initialized: type set in JSON driver not recognized.", md_, cl_); }
  }
}

// ┌─┐  ┬  ┬  ┌─┐  ┬    ╔╗   ╔═╗  ╔═╗  ╔═╗  ╔╦╗  ╔═╗  ╔╦╗  ╔═╗  ╦
// ├┤   └┐┌┘  ├─┤  │    ╠╩╗  ╠═╣  ╚═╗  ║╣   ║║║  ║ ║   ║║  ║╣   ║
// └─┘   └┘   ┴ ┴  ┴─┘  ╚═╝  ╩ ╩  ╚═╝  ╚═╝  ╩ ╩  ╚═╝  ═╩╝  ╚═╝  ╩═╝
void AbstractRunner::EvalBaseModel() {
  if (simulator_ == nullptr) {
    E("Simulator must be initialized b/f evaluating base model.", md_, cl_);
  }

  if (is_ensemble_run_) {
    im_ = "Simulating ensemble base case.";
    if (vp_.vRUN >= 1) { ext_info(im_, md_, cl_, vp_.lnw); }

    auto base_rlz = ensemble_helper_.GetBaseRealization();
    simulator_->Evaluate(base_rlz, 10000, 4);

  } else if (!simulator_->results()->isAvailable()) {
    im_ = "Simulating base case.";
    if (vp_.vRUN >= 1) { ext_info(im_, md_, cl_, vp_.lnw); }
    simulator_->Evaluate();
  }
}

// ┬  ┌┐┌  ┬  ┌┬┐    ╔═╗  ╔╗    ╦  ╔═╗
// │  │││  │   │     ║ ║  ╠╩╗   ║  ╠╣
// ┴  ┘└┘  ┴   ┴     ╚═╝  ╚═╝  ╚╝  ╚
void AbstractRunner::InitObjF() {
  if (simulator_ == nullptr || settings_ == nullptr) {
    E("Simulator & Settings must be initialized b/f ObjectiveFunction.", md_, cl_);
  }

  switch (seto_->objective().type) {

    case ObjTyp::WeightedSum: {
      auto tm = "Using WeightedSum-type objective function.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      objf_ = new WeightedSum(seto_, simulator_->results(), model_);
      break;
    }

    case ObjTyp::NPV: {
      auto tm = "Using NPV-type objective function.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      objf_ = new NPV(seto_, simulator_->results(), model_);
      break;
    }

    case ObjTyp::Augmented: {
      auto tm = "Using Augmented objective.";
      if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }
      objf_ = new Augmented(seto_, simulator_->results(), model_);
      break;
    }

    default: { E("Runner not initialized: ObjectiveFunction type not recognized.", md_, cl_); }
  }
}

// ┌─┐  ┌─┐  ┌─┐  ┬    ┬ ┬    ╦═╗  ╔═╗  ╔═╗  ╔╦╗  ╔═╗  ╦═╗  ╔╦╗
// ├─┤  ├─┘  ├─┘  │    └┬┘    ╠╦╝  ║╣   ╚═╗   ║   ╠═╣  ╠╦╝   ║
// ┴ ┴  ┴    ┴    ┴─┘   ┴     ╩╚═  ╚═╝  ╚═╝   ╩   ╩ ╩  ╩╚═   ╩
void AbstractRunner::ApplyRestartBaseCase() {
  if (base_case_ == nullptr || settings_ == nullptr)
    E("base_case_ and settings_ not initialized", md_, cl_);

  QJsonArray rjson0, rjson1;

  cout << "BestPoint" << endl;
  rjson0 = seto_->restartJson()->value("BestPoint").toArray();
  for (const auto v: rjson0)
    if (v.toObject().contains("Fmax")) {
      base_case_->set_objf_value(v.toObject().value("Fmax").toDouble());
    } else if (v.toObject().contains("Variables")) {
      rjson1 = v.toObject().value("Variables").toArray();
    }

  for (auto var: *vars_->GetContinuousVariables()) {
    base_case_->set_real_variable_value(var.first, seto_->infd_n());

    bool found_var = false;
    for (const auto v: rjson1)
      if (v.toObject().keys()[0].contains(var.second->name())) {
        auto oval = v.toObject().value("Var#" + var.second->name()).toDouble();
        var.second->scaleOtherValue(oval);
        base_case_->set_real_variable_value(var.first, oval);
        found_var = true; // one of the rstrt vars should match the current var
        break;
      }

    if (!found_var)
      ext_warn("input var.rstrt not found in current var.container", md_, cl_);
  }

  // double check all vars in var.container have been updated
  for (auto var : vars_->GetContVarValues())
    if (var.second == seto_->infd_n()) ext_warn("var.container not updated", md_, cl_);
}

// ┬  ┌┐┌  ┬  ┌┬┐    ╔╗   ╔═╗  ╔═╗  ╔═╗    ╔═╗  ╔═╗  ╔═╗  ╔═╗
// │  │││  │   │     ╠╩╗  ╠═╣  ╚═╗  ║╣     ║    ╠═╣  ╚═╗  ║╣
// ┴  ┘└┘  ┴   ┴     ╚═╝  ╩ ╩  ╚═╝  ╚═╝    ╚═╝  ╩ ╩  ╚═╝  ╚═╝
void AbstractRunner::InitBaseCase() {
  if (objf_ == nullptr || model_ == nullptr) {
    E("ObjectiveFunction and Model must be initialized before BaseCase.", md_, cl_);
  }
  base_case_ = new Optimization::Case(vars_->GetBinVarValues(),
                                      vars_->GetDiscVarValues(),
                                      vars_->GetContVarValues());
  base_case_->SetVerbParams(vp_);
  if (seto_->restart())
    ApplyRestartBaseCase();

  if (!simulator_->results()->isAvailable()) {
    base_case_->set_objf_value(seto_->infd_n());

    string tm = "Simulation results are unavailable.";
    tm += "Base case objective set to sentinel value.";
    if (vp_.vRUN >= 1) { info(tm, vp_.lnw); }

  } else {
    model_->wellCost(seto_);
    if (seto_->objective().type == ObjTyp::Augmented) {
      base_case_->set_objf_value(objf_->value(true));
    } else {
      base_case_->set_objf_value(objf_->value());
    }

    string tm = "Base case objective function value set to ";
    tm += num2str(base_case_->objf_value(), 8, 1);
    if (vp_.vRUN >= 1) { info(tm,  vp_.lnw); }
  }
}

// ┬  ┌┐┌  ┬  ┌┬┐    ╔═╗  ╔═╗  ╔╦╗  ╦  ╔╦╗  ╦  ╔═╗  ╔═╗  ╦═╗
// │  │││  │   │     ║ ║  ╠═╝   ║   ║  ║║║  ║  ╔═╝  ║╣   ╠╦╝
// ┴  ┘└┘  ┴   ┴     ╚═╝  ╩     ╩   ╩  ╩ ╩  ╩  ╚═╝  ╚═╝  ╩╚═
void AbstractRunner::InitializeOptimizer() {
  if (base_case_ == nullptr || model_ == nullptr) {
    E("Base Case and Model must be initialized before the Optimizer", md_, cl_);
  }
  seto_->setPaths(&settings_->paths());

  switch (seto_->type()) {

    // ┌─┐  ┌─┐  ┌┬┐  ┌─┐  ┌─┐  ┌─┐  ┌─┐
    // │    │ │  │││  ├─┘  ├─┤  └─┐  └─┐
    // └─┘  └─┘  ┴ ┴  ┴    ┴ ┴  └─┘  └─┘
    case OptzrTyp::Compass: {
      if (vp_.vRUN >= 1) { ext_info("Using CompassSearch.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::CompassSearch(seto_,
                                         base_case_,
                                         vars_,
                                         model_->grid(),
                                         logger_,
                                         nullptr,
                                         model_->constraintHandler());
      break;
    }
      // ┌─┐  ┌─┐  ┌─┐  ┌─┐
      // ├─┤  ├─┘  ├─┘  └─┐
      // ┴ ┴  ┴    ┴    └─┘
    case OptzrTyp::APPS: {
      if (vp_.vRUN >= 1) { ext_info("Using APPS.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::APPS(seto_,
                                base_case_,
                                vars_,
                                model_->grid(),
                                logger_,
                                nullptr,
                                model_->constraintHandler());
      break;
    }
      // ┌─┐  ┌─┐  ┌─┐
      // ├─┘  └─┐  │ │
      // ┴    └─┘  └─┘
    case OptzrTyp::PSO: {
      if (vp_.vRUN >= 1) { ext_info("Using PSO.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::PSO(seto_,
                               base_case_,
                               vars_,
                               model_->grid(),
                               logger_,
                               nullptr,
                               model_->constraintHandler());
      break;
    }
      // ┌┬┐  ┬─┐  ┌┬┐  ┌─┐  ┌─┐
      //  │   ├┬┘   ││  ├┤   │ │
      //  ┴   ┴└─  ─┴┘  └    └─┘
    case OptzrTyp::TrustRegionOptimization: {
      if (vp_.vRUN >= 1) { ext_info("Using TR-DFO.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::TrustRegionOptimization(seto_,
                                                   base_case_,
                                                   vars_,
                                                   model_->grid(),
                                                   logger_,
                                                   nullptr,
                                                   model_->constraintHandler());
      break;
    }
      // ┌┬┐  ┌─┐  ┌┬┐  ┬─┐
      //  ││  ├┤    │   ├┬┘
      // ─┴┘  └     ┴   ┴└─
    case OptzrTyp::DFTR: {
      if (vp_.vRUN >= 1) { ext_info("Using DFTR [TR-DFO].", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::DFTR(seto_,
                                base_case_,
                                vars_,
                                model_->grid(),
                                logger_,
                                nullptr,
                                model_->constraintHandler());
      break;
    }
      // ┌─┐  ┌─┐   ┌─┐
      // └─┐  │─┼┐  ├─┘
      // └─┘  └─┘└  ┴
    case OptzrTyp::SQP: {
      if (vp_.vRUN >= 1) { ext_info("Using SQP [SNOPT].", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::SQP_SNOPT(seto_,
                                     base_case_,
                                     model_,
                                     simulator_,
                                     objf_,
                                     logger_,
                                     nullptr,
                                     model_->constraintHandler());
      break;
    }
      // ┌─┐  ┌─┐
      // │ ┬  ├─┤
      // └─┘  ┴ ┴
    case OptzrTyp::GeneticAlgorithm: {
      if (vp_.vRUN >= 1) { ext_info("Using GA.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::RGARDD(seto_,
                                  base_case_,
                                  vars_,
                                  model_->grid(),
                                  logger_);
      break;
    }
      // ┌─┐  ┌─┐  ┌─┐
      // ├┤   │ ┬  │ │
      // └─┘  └─┘  └─┘
    case OptzrTyp::EGO: {
      if (vp_.vRUN >= 1) { ext_info("Using EGO.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::BayesianOptimization::EGO(seto_,
                                                     base_case_,
                                                     vars_,
                                                     model_->grid(),
                                                     logger_);
      break;
    }
      // ┌─┐  ─┐ ┬  ┬ ┬  ┌─┐  ┬─┐  ┌─┐  ┬ ┬  ┬┬  ┌┬┐
      // ├┤   ┌┴┬┘  ├─┤  └─┐  ├┬┘  │    ├─┤  ││   ││
      // └─┘  ┴ └─  ┴ ┴  └─┘  ┴└─  └─┘  ┴ ┴  ┴┴  ─┴┘
    case OptzrTyp::ExhaustiveSearch2DVert: {
      if (vp_.vRUN >= 1) { ext_info("Using ExhaustiveSearch2DVert.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::ExhaustiveSearch2DVert(seto_,
                                                  base_case_,
                                                  vars_,
                                                  model_->grid(),
                                                  logger_);
      break;
    }
      // ┬ ┬  ┬ ┬  ┌┐   ┬─┐  ┬  ┌┬┐
      // ├─┤  └┬┘  ├┴┐  ├┬┘  │   ││
      // ┴ ┴   ┴   └─┘  ┴└─  ┴  ─┴┘
    case OptzrTyp::Hybrid: {
      if (vp_.vRUN >= 1) { ext_info("Using Hybrid optimization.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optimization::HybridOptimizer(seto_,
                                                  base_case_,
                                                  vars_,
                                                  model_->grid(),
                                                  logger_);
      break;
    }
      // ┌─┐  ┌┬┐  ┌─┐       ┌─┐  ┌─┐
      // │    │││  ├─┤  ───  ├┤   └─┐
      // └─┘  ┴ ┴  ┴ ┴       └─┘  └─┘
    case OptzrTyp::CMA_ES: {
      if (vp_.vRUN >= 1) { ext_info("Using CMA_ES.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::CMA_ES(seto_,
                                  base_case_,
                                  vars_,
                                  model_->grid(),
                                  logger_);
      break;
    }
      // ┬  ┬  ┌─┐  ┌─┐  ┌─┐
      // └┐┌┘  ├┤   └─┐  ├─┤
      //  └┘   └    └─┘  ┴ ┴
    case OptzrTyp::VFSA: {
      if (vp_.vRUN >= 1) { ext_info("Using VFSA.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::VFSA(seto_,
                                base_case_,
                                vars_,
                                model_->grid(),
                                logger_);
      break;
    }
      // ┌─┐  ┌─┐  ┌─┐  ┌─┐
      // └─┐  ├─┘  └─┐  ├─┤
      // └─┘  ┴    └─┘  ┴ ┴
    case OptzrTyp::SPSA: {
      if (vp_.vRUN >= 1) { ext_info("Using SPSA.", md_, cl_, vp_.lnw); }
      optmzr_ = new Optzr::SPSA(seto_,
                                base_case_,
                                vars_,
                                model_->grid(),
                                logger_);
      break;
    }
    default: {
      string em = "Unable to initialize runner: ";
      em += "Algorithm type in driver file not recognized.";
      E(em, md_, cl_);
    }

  }
  optmzr_->EnableConstraintLogging(
    rts_->paths().GetPathQstr(Paths::OUTPUT_DIR));
}

// ┬  ┌┐┌  ┬  ┌┬┐  ╔╗   ╔═╗  ╔═╗  ╦╔═  ╦╔═  ╔═╗  ╔═╗  ╔═╗  ╔═╗  ╦═╗
// │  │││  │   │   ╠╩╗  ║ ║  ║ ║  ╠╩╗  ╠╩╗  ║╣   ║╣   ╠═╝  ║╣   ╠╦╝
// ┴  ┘└┘  ┴   ┴   ╚═╝  ╚═╝  ╚═╝  ╩ ╩  ╩ ╩  ╚═╝  ╚═╝  ╩    ╚═╝  ╩╚═
void AbstractRunner::InitBookkeeper() {
  if (settings_ == nullptr || optmzr_ == nullptr) {
    E("Settings and Optimizer must be initialized before the Bookkeeper.", md_, cl_);
  }
  bookkeeper_ = new Bookkeeper(settings_, optmzr_->case_handler());
}

// ┬  ┌┐┌  ┬  ┌┬┐  ╦    ╔═╗  ╔═╗  ╔═╗  ╔═╗  ╦═╗
// │  │││  │   │   ║    ║ ║  ║ ╦  ║ ╦  ║╣   ╠╦╝
// ┴  ┘└┘  ┴   ┴   ╩═╝  ╚═╝  ╚═╝  ╚═╝  ╚═╝  ╩╚═
void AbstractRunner::InitLogger(QString output_subdir, bool write_logs) {
  logger_ = new Logger(rts_, output_subdir, write_logs);
}

// ┌─┐  ┬─┐  ┌┐┌  ┌┬┐  ╔═╗  ╔═╗  ╔╦╗  ╔═╗  ╦    ╔═╗  ╔╦╗  ╦  ╔═╗  ╔╗╔  ╔╦╗  ╔═╗  ╔═╗
// ├─┘  ├┬┘  │││   │   ║    ║ ║  ║║║  ╠═╝  ║    ║╣    ║   ║  ║ ║  ║║║  ║║║  ╚═╗  ║ ╦
// ┴    ┴└─  ┘└┘   ┴   ╚═╝  ╚═╝  ╩ ╩  ╩    ╩═╝  ╚═╝   ╩   ╩  ╚═╝  ╝╚╝  ╩ ╩  ╚═╝  ╚═╝
void AbstractRunner::PrintCompletionMessage() {
  cout << "Optimization complete: ";

  switch (optmzr_->IsFinished()) {

    case TC::MAX_ITERS_REACHED:
      cout << "Max # iterations reached." << endl;
      break;

    case TC::MAX_EVALS_REACHED:
      cout << "Max # fevals reached." << endl;
      break;

    case TC::MIN_STEP_LENGTH_REACHED:
      cout << "Min step length reached (converged)." << endl;
      break;

    case TC::OPT_CRITERIA_REACHED:
      cout << "Opt criteria reached" << endl;
      break;

    case TC::DFTR_CRIT_NORM_1ST_ORD_LT_TOLF:
      cout << "Crit.norm[1st-order] less than tolf." << endl;
      break;

    case TC::DFTR_MAX_NUM_RHO_INF_MET:
      cout << "Max # of rho=-inf met." << endl;
      break;

    default: cout << "Unknown termination reason." << endl;
  }

  cout << "Best case at termination:" << endl;
  cout << optmzr_->GetTentBestCase()->id().toString().toStdString() << endl;
  cout << "Variable values: " << endl;

  auto opt_vars_int = optmzr_->GetTentBestCase()->integer_variables();
  auto opt_vars_real = optmzr_->GetTentBestCase()->real_variables();
  auto opt_vars_bin = optmzr_->GetTentBestCase()->binary_variables();

  for (auto var : opt_vars_int) {
    auto prop_name = vars_->GetDiscreteVariable(var.first)->name();
    cout << "\t" << prop_name.toStdString() << "\t" << var.second << endl;
  }

  for (auto var : opt_vars_real) {
    auto prop_name = vars_->GetContinuousVariable(var.first)->name();
    auto vre = vars_->GetContinuousVariable(var.first)->value();
    auto vsc = vars_->GetContinuousVariable(var.first)->valueSc();
    cout << "\t" << prop_name.toStdString();
    cout << "\t" << num2str(vsc, 3, 1);
    cout << "\t" << num2str(vre, 3, 1) << endl;
  }

  for (auto var : opt_vars_bin) {
    auto prop_name = vars_->GetBinaryVariable(var.first)->name();
    cout << "\t" << prop_name.toStdString() << "\t" << var.second << endl;
  }

  ComputeOptmzdCase();
}

// ┌─┐  ┌─┐  ┌┬┐  ┌─┐  ┬ ┬  ┌┬┐  ┌─┐  ╔═╗  ╔═╗  ╔╦╗  ╔╦╗  ╔═╗  ╔╦╗  ╔═╗  ╔═╗
// │    │ │  │││  ├─┘  │ │   │   ├┤   ║ ║  ╠═╝   ║   ║║║  ╔═╝   ║║  ║    ╚═╗
// └─┘  └─┘  ┴ ┴  ┴    └─┘   ┴   └─┘  ╚═╝  ╩     ╩   ╩ ╩  ╚═╝  ═╩╝  ╚═╝  ╚═╝
void AbstractRunner::ComputeOptmzdCase() {
  cout << "Running simulation using optzd values" << endl;
  auto opt_vars_int = optmzr_->GetTentBestCase()->integer_variables();
  auto opt_vars_real = optmzr_->GetTentBestCase()->real_variables();
  auto opt_vars_bin = optmzr_->GetTentBestCase()->binary_variables();
  optz_case_ = new Optimization::Case(opt_vars_bin, opt_vars_int, opt_vars_real);

  settings_->paths().CopyPath(Paths::OUTPUT_DIR, Paths::OPTMZD_DIR);
  cout << "OUTPUT_DIR:" << settings_->paths().GetPath(Paths::OUTPUT_DIR) << endl;

  model_->ApplyCase(optz_case_);
  simulator_->UpdatePaths(settings_->paths());
  simulator_->Evaluate();

  if (seto_->objective().type == ObjTyp::Augmented) {
    optz_case_->set_objf_value(objf_->value(true));
  } else {
    optz_case_->set_objf_value(objf_->value());
  }
  cout << "objf_value: " << optz_case_->objf_value() << endl;
}

// ┌┬┐  ┬  ┌┬┐  ┌─┐  ┌─┐  ┬ ┬  ┌┬┐  ╦  ╦  ╔═╗  ╦
//  │   │  │││  ├┤   │ │  │ │   │   ╚╗╔╝  ╠═╣  ║
//  ┴   ┴  ┴ ┴  └─┘  └─┘  └─┘   ┴    ╚╝   ╩ ╩  ╩═╝
int AbstractRunner::timeoutVal() const {
  if (sim_times_.empty() || rts_->sim_timeout() == 0) {
    return 10000;
  } else {
    return calc_median(sim_times_) * rts_->sim_timeout();
  }
}

// ┌─┐  ┬  ┌┐┌  ┌─┐  ┬    ┬  ┌─┐  ┌─┐  ╦  ╔╗╔  ╦  ╔╦╗
// ├┤   │  │││  ├─┤  │    │  ┌─┘  ├┤   ║  ║║║  ║   ║
// └    ┴  ┘└┘  ┴ ┴  ┴─┘  ┴  └─┘  └─┘  ╩  ╝╚╝  ╩   ╩
void AbstractRunner::FinalizeInit(bool write_logs) {
  if (write_logs) {
    logger_->AddEntry(rts_);
    logger_->FinalizePrerunSummary();
  }
}

// ┌─┐  ┬  ┌┐┌  ┌─┐  ┬    ┬  ┌─┐  ┌─┐  ╦═╗  ╦ ╦  ╔╗╔
// ├┤   │  │││  ├─┤  │    │  ┌─┘  ├┤   ╠╦╝  ║ ║  ║║║
// └    ┴  ┘└┘  ┴ ┴  ┴─┘  ┴  └─┘  └─┘  ╩╚═  ╚═╝  ╝╚╝
void AbstractRunner::FinalizeRun(bool write_logs) {
  if (optmzr_ != nullptr) { // This indicates whether or not we're on a worker process
    model_->ApplyCase(optmzr_->GetTentBestCase());
    simulator_->WriteDriverFilesOnly();
    PrintCompletionMessage();
  }
  cout << "model_->Finalize() " << endl;
  model_->Finalize();
  if (write_logs)
    logger_->FinalizePostrunSummary();
}

void AbstractRunner::ReplaceOptimizer(Settings::Optimizer *opt) {
  settings_->setOptimizer(opt);
  InitializeOptimizer();
  InitializeBookkeeper();
  FinalizeInitialization(true);
}

}
