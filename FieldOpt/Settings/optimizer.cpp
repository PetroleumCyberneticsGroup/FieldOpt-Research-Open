/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb.pcg@gmail.com>

Modified 2019-2020 Brage Strand Kristoffersen
<brage.s.kristoffersen@ntnu.no>

Modified 2019-2020 Thiago Lima Silva
<thiagolims@gmail.com>

Modified 2019-2020 Caio Giuliani
<caiogiuliani@gmail.com>

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

#include "Utilities/printer.hpp"
#include "optimizer.h"

#include "settings_exceptions.h"
#include "Settings/helpers.hpp"

namespace Settings {

using std::runtime_error;
using Printer::ext_warn;
using Printer::E;

Optimizer::Optimizer(QJsonObject json_optimizer, VerbParams vp) {
  vp_=vp;

  // Get the root objects.
  QJsonObject json_parameters = json_optimizer["Parameters"].toObject();
  QJsonObject json_objective = json_optimizer["Objective"].toObject();
  QJsonArray json_constraints = json_optimizer["Constraints"].toArray();
  QString type = json_optimizer["Type"].toString();

  scale_vars_ = json_optimizer["ScaleVars"].toInt();
  scale_objf_ = json_optimizer["ScaleObjf"].toDouble();

  type_ = parseType(type);

  if (type_ != ExhaustiveSearch2DVert) {
    mode_ = parseMode(json_optimizer);
    if (type_ == Hybrid) {
      hybrid_components_ = parseHybridComponents(json_optimizer);
    }
    parameters_ = parseParameters(json_parameters);
  }
  objective_ = parseObjective(json_objective);


  // Optimizer constraints
  try {
    constraints_ = QList<Constraint>();
    // Iterate over all constraints
    for (int i = 0; i < json_constraints.size(); ++i) {
      QJsonObject json_constraint = json_constraints[i].toObject();
      constraints_.append(parseSingleConstraint(json_constraint));
    }
  } catch (std::exception const &ex) {
    throw UnableToParseOptimizerConstraintsSectionException(
        "Unable to parse optimizer constraints: "
            + std::string(ex.what()));
  }
}

Optimizer::Constraint
Optimizer::parseSingleConstraint(QJsonObject json_constraint) {

  Constraint optmzr_constraint = Constraint();

  optmzr_constraint.scaling_ = scale_vars_;

  if (json_constraint.contains("Well")) {
    optmzr_constraint.well = json_constraint["Well"].toString();
    optmzr_constraint.wells.append(optmzr_constraint.well);

  } else if (json_constraint.contains("Wells") &&
      json_constraint["Wells"].isArray()) {

    if (json_constraint["Wells"].toArray().size() == 1) {
      optmzr_constraint.well =
          json_constraint["Wells"].toArray()[0].toString();
    }

    for (auto wname : json_constraint["Wells"].toArray()) {
      optmzr_constraint.wells.append(wname.toString());
    }

  } else {
    em_ = "A constraint must always specify either the Well";
    em_ = " or the Wells property.";
    throw runtime_error(em_);
  }

  // Penalty function weight for the constraint
  if (json_constraint.contains("PenaltyWeight")) {
    optmzr_constraint.penalty_weight =
        json_constraint["PenaltyWeight"].toDouble();

  } else {
    optmzr_constraint.penalty_weight = 0.0;
  }

  // Constraint types BHP, Rate and Boundary2D
  QString constraint_type = json_constraint["Type"].toString();
  if (QString::compare(constraint_type, "BHP") == 0) {

    optmzr_constraint.type = ConstraintType::BHP;
    if (json_constraint.contains("Max")) {
      set_opt_prop_double(optmzr_constraint.max,
                          json_constraint, "Max", vp_);
    }
    if (json_constraint.contains("Min")) {
      set_opt_prop_double(optmzr_constraint.min,
                          json_constraint, "Min", vp_);
    }

    // Packer- and ICV Constraints
  } else if (QString::compare(constraint_type, "ICVConstraint") == 0) {

    optmzr_constraint.type = ConstraintType::ICVConstraint;
    set_opt_prop_double(optmzr_constraint.max,
                        json_constraint, "Max", vp_);

    if (optmzr_constraint.max >= 7.8540E-03) {
      string wm = "Maximum valve size is too big. Setting it to 7.8539-3.";
      ext_warn(wm, md_, cl_, vp_.lnw);
      optmzr_constraint.max = 7.8539E-3;
    }
    set_opt_prop_double(optmzr_constraint.min,
                        json_constraint, "Min", vp_);

  } else if (CnstrCmp(constraint_type, "PackerConstraint")) {
    optmzr_constraint.type = ConstraintType::PackerConstraint;

  } else if (QString::compare(constraint_type, "Rate") == 0) {
    optmzr_constraint.type = ConstraintType::Rate;
    if (json_constraint.contains("Max")) {
      set_opt_prop_double(optmzr_constraint.max,
                          json_constraint, "Max", vp_);
    }
    if (json_constraint.contains("Min")) {
      set_opt_prop_double(optmzr_constraint.min,
                          json_constraint, "Min", vp_);
    }

  } else if (CnstrCmp(constraint_type, "Boundary2D")) {

    optmzr_constraint.type = ConstraintType::PseudoContBoundary2D;
    set_opt_prop_double(optmzr_constraint.box_imin, json_constraint, "Imin", vp_);
    set_opt_prop_double(optmzr_constraint.box_imax, json_constraint, "Imax", vp_);
    set_opt_prop_double(optmzr_constraint.box_jmin, json_constraint, "Jmin", vp_);
    set_opt_prop_double(optmzr_constraint.box_jmax, json_constraint, "Jmax", vp_);

    // Constraint type Well Spline Points
  } else if (CnstrCmp(constraint_type, "WellSplinePoints")) {

    optmzr_constraint.type = ConstraintType::SplinePoints;

    // Spline points constraint input type
    QString optmzr_cnstrnt_spline_pnts_type =
        json_constraint["WellSplinePointsInputType"].toString();

    if (CnstrCmp(optmzr_cnstrnt_spline_pnts_type, "Function")) {

      optmzr_constraint.spline_points_type = ConstraintWellSplinePointsType::Function;
      json_constraint["Function"].toString();

    } else if (CnstrCmp(optmzr_cnstrnt_spline_pnts_type, "MaxMin")) {

      optmzr_constraint.spline_points_type =
          ConstraintWellSplinePointsType::MaxMin;

      optmzr_constraint.spline_points_limits =
          QList<Constraint::RealMaxMinLimit>();

      QJsonArray well_spline_point_limits =
          json_constraint["WellSplinePointLimits"].toArray();

      for (int i = 0; i < well_spline_point_limits.size(); ++i) {
        QJsonObject well_spline_point_limit = well_spline_point_limits[i].toObject();
        QJsonArray min_array = well_spline_point_limit["Min"].toArray();
        QJsonArray max_array = well_spline_point_limit["Max"].toArray();

        Constraint::RealCoordinate min{};
        min.x = min_array[0].toDouble();
        min.y = min_array[1].toDouble();
        min.z = min_array[2].toDouble();
        Constraint::RealCoordinate max{};
        max.x = max_array[0].toDouble();
        max.y = max_array[1].toDouble();
        max.z = max_array[2].toDouble();

        Constraint::RealMaxMinLimit limit{};
        limit.min = min;
        limit.max = max;
        optmzr_constraint.spline_points_limits.append(limit);
      }

    } else {
      throw UnableToParseOptimizerConstraintsSectionException(
          "Well spline constraint type not recognized.");
    }

  } else if (QString::compare(constraint_type, "WSplineLength") == 0
      || QString::compare(constraint_type, "PolarWellLength") == 0 ) {

    if (constraint_type == "WSplineLength"){
      optmzr_constraint.type = ConstraintType::WSplineLength;
    } else {
      optmzr_constraint.type = ConstraintType::PolarWellLength;
    }

    if (json_constraint.contains("Min")) {
      optmzr_constraint.min = json_constraint["Min"].toDouble();
      optmzr_constraint.min_length = json_constraint["Min"].toDouble();

    } else if (json_constraint.contains("MinLength")) {
      optmzr_constraint.min = json_constraint["MinLength"].toDouble();
      optmzr_constraint.min_length = json_constraint["MinLength"].toDouble();
    }
    else {
      throw std::runtime_error(
          "The MinLength field must be specified for well spline length constraints.");
    }

    if (json_constraint.contains("Max")) {
      optmzr_constraint.max = json_constraint["Max"].toDouble();
      optmzr_constraint.max_length = json_constraint["Max"].toDouble();

    } else if (json_constraint.contains("MaxLength")) {
      optmzr_constraint.max = json_constraint["MaxLength"].toDouble();
      optmzr_constraint.max_length = json_constraint["MaxLength"].toDouble();

    } else {
      throw std::runtime_error(
          "MaxLength must be specified for well length constraints.");
    }

  } else if (QString::compare(constraint_type, "WSplineInterwDist") == 0) {

    optmzr_constraint.type = ConstraintType::WSplineInterwDist;

    if (json_constraint.contains("Min")) {
      optmzr_constraint.min = json_constraint["Min"].toDouble();
      optmzr_constraint.min_distance = json_constraint["Min"].toDouble();

    } else if (json_constraint.contains("MinDistance")) {
      optmzr_constraint.min = json_constraint["MinDistance"].toDouble();
      optmzr_constraint.min_distance = json_constraint["MinDistance"].toDouble();
    }

    if (optmzr_constraint.wells.length() != 2) {
      throw UnableToParseOptimizerConstraintsSectionException(
          "WSplineInterwDist constraint needs"
          " a Wells array with exactly two well names specified.");
    }

  } else if (QString::compare(constraint_type, "PolarAzimuth") == 0
      || QString::compare(constraint_type, "PolarElevation") == 0) {

    if (constraint_type == "PolarAzimuth") {
      optmzr_constraint.type = ConstraintType::PolarAzimuth;
    } else {
      optmzr_constraint.type = ConstraintType::PolarElevation;
    }
    if (json_constraint.contains("Max")){
      optmzr_constraint.max = json_constraint["Max"].toDouble();
    }
    if (json_constraint.contains("Min")){
      optmzr_constraint.min = json_constraint["Min"].toDouble();
    }

  } else if (QString::compare(constraint_type, "ResBoundary") == 0
      || QString::compare(constraint_type, "PolarSplineBoundary") == 0
      || QString::compare(constraint_type, "ReservoirBoundaryToe") == 0) {

    if (QString::compare(constraint_type, "ResBoundary") == 0){
      optmzr_constraint.type = ConstraintType::ReservoirBoundary;

    } else if (QString::compare(constraint_type, "PolarSplineBoundary") == 0){
      optmzr_constraint.type = ConstraintType::PolarSplineBoundary;

    } else if (QString::compare(constraint_type, "ReservoirBoundaryToe") == 0){
      optmzr_constraint.type = ConstraintType::ReservoirBoundaryToe;
    }

    optmzr_constraint.box_imin = json_constraint["BoxImin"].toInt();
    optmzr_constraint.box_imax = json_constraint["BoxImax"].toInt();
    optmzr_constraint.box_jmin = json_constraint["BoxJmin"].toInt();
    optmzr_constraint.box_jmax = json_constraint["BoxJmax"].toInt();
    optmzr_constraint.box_kmin = json_constraint["BoxKmin"].toInt();
    optmzr_constraint.box_kmax = json_constraint["BoxKmax"].toInt();

  } else if (QString::compare(constraint_type, "PolarXYZBoundary") == 0
      || QString::compare(constraint_type, "ReservoirXYZBoundary") == 0
      || QString::compare(constraint_type, "WellXYZBox") == 0) {


    if(QString::compare(constraint_type, "PolarXYZBoundary") == 0) {
      optmzr_constraint.type = ConstraintType::PolarXYZBoundary;

    } else if(QString::compare(constraint_type, "ReservoirXYZBoundary") == 0) {
      optmzr_constraint.type = ConstraintType::ReservoirXYZBoundary;

    } else if(QString::compare(constraint_type, "WellXYZBox") == 0) {
      optmzr_constraint.type = ConstraintType::WellXYZBox;
    }

    optmzr_constraint.box_xyz_xmin = json_constraint["xMin"].toDouble();
    optmzr_constraint.box_xyz_ymin = json_constraint["yMin"].toDouble();
    optmzr_constraint.box_xyz_zmin = json_constraint["zMin"].toDouble();
    optmzr_constraint.box_xyz_xmax = json_constraint["xMax"].toDouble();
    optmzr_constraint.box_xyz_ymax = json_constraint["yMax"].toDouble();
    optmzr_constraint.box_xyz_zmax = json_constraint["zMax"].toDouble();

  } else if (QString::compare(constraint_type,
                              "MxWSplineLengthInterwDist") == 0) {

    optmzr_constraint.type = ConstraintType::MxWSplineLengthInterwDist;
    optmzr_constraint.min_length = json_constraint["MinLength"].toDouble();
    optmzr_constraint.max_length = json_constraint["MaxLength"].toDouble();
    optmzr_constraint.min_distance = json_constraint["MinDistance"].toDouble();
    optmzr_constraint.max_iterations = json_constraint["MaxIterations"].toInt();

    if (optmzr_constraint.wells.length() != 2) {
      throw UnableToParseOptimizerConstraintsSectionException(
          "WSplineInterwDist constraint"
          " needs a Wells array with exactly two well names specified.");
    }

  } else if (QString::compare(constraint_type,
                              "MxWSplineLengthInterwDistResBound") == 0) {

    optmzr_constraint.type = ConstraintType::MxWSplineLengthInterwDistResBound;
    optmzr_constraint.min_length = json_constraint["MinLength"].toDouble();
    optmzr_constraint.max_length = json_constraint["MaxLength"].toDouble();
    optmzr_constraint.min_distance = json_constraint["MinDistance"].toDouble();
    optmzr_constraint.max_iterations = json_constraint["MaxIterations"].toInt();

    optmzr_constraint.box_imin = json_constraint["BoxImin"].toInt();
    optmzr_constraint.box_imax = json_constraint["BoxImax"].toInt();
    optmzr_constraint.box_jmin = json_constraint["BoxJmin"].toInt();
    optmzr_constraint.box_jmax = json_constraint["BoxJmax"].toInt();
    optmzr_constraint.box_kmin = json_constraint["BoxKmin"].toInt();
    optmzr_constraint.box_kmax = json_constraint["BoxKmax"].toInt();

    if (optmzr_constraint.wells.length() != 2) {
      string em = "MxWSplineLengthInterwDistResBound constraint ";
      em += "needs a Wells array with exactly two well names specified.";
      throw UnableToParseOptimizerConstraintsSectionException(em);
    }

  } else {
    throw UnableToParseOptimizerConstraintsSectionException(
        "Constraint type " + constraint_type.toStdString()
            + " not recognized.");
  }
  return optmzr_constraint;
}

Optimizer::OptimizerMode Optimizer::parseMode(QJsonObject &json_optimizer) {
  OptimizerMode opt_mode;
  if (json_optimizer.contains("Mode")) {
    QString mode = json_optimizer["Mode"].toString();
    if (QString::compare(mode, "Minimize", Qt::CaseInsensitive) == 0) {
      opt_mode = OptimizerMode::Minimize;
    } else if (QString::compare(mode, "Maximize", Qt::CaseInsensitive) == 0) {
      opt_mode = OptimizerMode::Maximize;
    } else {
      throw UnableToParseOptimizerSectionException(
          "Optimizer Mode keyword must be specified.");
    }
  }
  return opt_mode;
}


Optimizer::Parameters Optimizer::parseParameters(QJsonObject &json_params) {
  Parameters params;
  string ind_val = "Invalid value for setting";
  string not_rec = "" + not_rec;

  try {

    // -----------------------------------------------------
    // GSS PARAMETERS
    // GSS parameters :: MaxEvaluations
    if (json_params.contains("MaxEvaluations"))
      params.max_evaluations = json_params["MaxEvaluations"].toInt();

    // GSS parameters :: AutoStepLengths
    if (json_params.contains("AutoStepLengths"))
      params.auto_step_lengths = json_params["AutoStepLengths"].toBool();

    // GSS parameters :: AutoStepInitScale
    if (json_params.contains("AutoStepInitScale"))
      params.auto_step_init_scale = json_params["AutoStepInitScale"].toDouble();

    // GSS parameters :: AutoStepConvScale
    if (json_params.contains("AutoStepConvScale"))
      params.auto_step_conv_scale = json_params["AutoStepConvScale"].toDouble();

    // GSS parameters :: InitialStepLength
    if (json_params.contains("InitialStepLength")) {
      set_req_prop_double(params.initial_step_length, json_params, "InitialStepLength", vp_);
    }

    // GSS parameters :: MinimumStepLength
    if (json_params.contains("MinimumStepLength")) {
      set_req_prop_double(params.minimum_step_length, json_params, "MinimumStepLength", vp_);
    }

    // GSS parameters :: ContractionFactor
    if (json_params.contains("ContractionFactor"))
      params.contraction_factor = json_params["ContractionFactor"].toDouble();
    else { params.contraction_factor = 0.5; }

    // GSS parameters :: ExpansionFactor
    if (json_params.contains("ExpansionFactor"))
      params.expansion_factor = json_params["ExpansionFactor"].toDouble();
    else { params.expansion_factor = 1.0; }

    // GSS parameters :: MaxQueueSize
    if (json_params.contains("MaxQueueSize"))
      params.max_queue_size = json_params["MaxQueueSize"].toInt();
    else { params.max_queue_size = 2; }

    // GSS parameters :: Pattern
    if (json_params.contains("Pattern"))
      params.pattern = json_params["Pattern"].toString();
    else params.pattern = "Compass";

    // -----------------------------------------------------
    // GA PARAMETERS
    // GA parameters :: MaxGenerations
    if (json_params.contains("MaxGenerations"))
      params.max_generations = json_params["MaxGenerations"].toInt();
    else { params.max_generations = 50; }

    // GA parameters :: PopulationSize
    if (json_params.contains("PopulationSize"))
      params.population_size = json_params["PopulationSize"].toInt();
    else { params.population_size = -1; } // Will be properly set in optimizer.

    // GA parameters :: CrossoverProbability
    if (json_params.contains("CrossoverProbability"))
      params.p_crossover = json_params["CrossoverProbability"].toDouble();
    else { params.p_crossover = 0.1; }

    // GA parameters :: DiscardParameter
    if (json_params.contains("DiscardParameter"))
      params.discard_parameter = json_params["DiscardParameter"].toDouble();
    else { params.discard_parameter = -1; } // Will be properly set in optimizer

    // GA parameters :: DecayRate
    if (json_params.contains("DecayRate"))
      params.decay_rate = json_params["DecayRate"].toDouble();
    else { params.decay_rate = 4.0; }

    // GA parameters :: MutationStrength
    if (json_params.contains("MutationStrength"))
      params.mutation_strength = json_params["MutationStrength"].toDouble();
    else { params.mutation_strength = 0.25; }

    // GA parameters :: StagnationLimit
    if (json_params.contains("StagnationLimit"))
      params.stagnation_limit = json_params["StagnationLimit"].toDouble();
    else { params.stagnation_limit = 1e-10; }

    // GA parameters :: LowerBound
    if (json_params.contains("LowerBound"))
      params.lower_bound = json_params["LowerBound"].toDouble();
    else { params.lower_bound = -10; }

    // GA parameters :: UpperBound
    if (json_params.contains("UpperBound"))
      params.upper_bound = json_params["UpperBound"].toDouble();
    else { params.upper_bound = 10; }

    // -----------------------------------------------------
    // PSO PARAMETERS
    // PSO parameters :: PSO-LearningFactor1
    if(json_params.contains("PSO-LearningFactor1")){
      params.pso_learning_factor_1 = json_params["PSO-LearningFactor1"].toDouble();
    } else { params.pso_learning_factor_1 = 2; }

    // PSO parameters :: PSO-LearningFactor2
    if(json_params.contains("PSO-LearningFactor2")){
      params.pso_learning_factor_2 = json_params["PSO-LearningFactor2"].toDouble();
    } else { params.pso_learning_factor_2 = 2; }

    // PSO parameters :: PSO-SwarmSize
    if(json_params.contains("PSO-SwarmSize")){
      params.pso_swarm_size = json_params["PSO-SwarmSize"].toDouble();
    } else { params.pso_swarm_size = 50; }

    // PSO parameters :: PSO-VelocityScale
    if(json_params.contains("PSO-VelocityScale")){
      params.pso_velocity_scale = json_params["PSO-VelocityScale"].toDouble();
    } else { params.pso_velocity_scale = 1.0; }

    // -----------------------------------------------------
    // EGO PARAMEETERS
    // EGO Parameters :: EGO-InitGuesses
    if (json_params.contains("EGO-InitGuesses")) {
      params.ego_init_guesses = json_params["EGO-InitGuesses"].toInt();
    }

    // EGO Parameters :: EGO-InitSamplingMethod
    if (json_params.contains("EGO-InitSamplingMethod")) {
      QString method = json_params["EGO-InitSamplingMethod"].toString();
      if (QString::compare(method, "Random") == 0 || QString::compare(method, "Uniform") == 0) {
        params.ego_init_sampling_method = json_params["EGO-InitSamplingMethod"].toString().toStdString();
      } else {
        Printer::error("EGO-InitSamplingMethod " + method.toStdString() + " not recognized.");
        throw std::runtime_error("Failed reading EGO settings.");
      }
    }

    // EGO Parameters :: EGO-Kernel
    if (json_params.contains("EGO-Kernel")) {
      QStringList available_kernels = { "CovLinearard", "CovLinearone", "CovMatern3iso",
                                        "CovMatern5iso", "CovNoise", "CovRQiso", "CovSEard",
                                        "CovSEiso", "CovPeriodicMatern3iso", "CovPeriodic"};
      if (available_kernels.contains(json_params["EGO-Kernel"].toString())) {
        params.ego_kernel = json_params["EGO-Kernel"].toString().toStdString();
      } else {
        Printer::error("EGO-Kernel " + json_params["EGO-Kernel"].toString().toStdString() + " not recognized.");
        Printer::info("Available kernels: " + available_kernels.join(", ").toStdString());
        throw std::runtime_error("Failed reading EGO settings.");
      }
    }

    // EGO Parameters :: EGO-AF
    if (json_params.contains("EGO-AF")) {
      QStringList available_afs = { "ExpectedImprovement", "ProbabilityOfImprovement" };
      if (available_afs.contains(json_params["EGO-AF"].toString())) {
        params.ego_af = json_params["EGO-AF"].toString().toStdString();

      } else {
        Printer::error("EGO-AF " + json_params["EGO-AF"].toString().toStdString() + " not recognized.");
        Printer::info("Available acquisition functions: " + available_afs.join(", ").toStdString());
        throw std::runtime_error("Failed reading EGO settings.");
      }
    }

    // -----------------------------------------------------
    // CMA-ES PARAMETERS
    if (json_params.contains("ImproveBaseCase")) {
      params.improve_base_case = json_params["ImproveBaseCase"].toBool();
    }

    // -----------------------------------------------------
    // VFSA PARAMETERS
    if (json_params.contains("VFSA-EvalsPrIteration")) {
      params.vfsa_evals_pr_iteration = json_params["VFSA-EvalsPrIteration"].toInt();
    }
    if (json_params.contains("VFSA-MaxIterations")) {
      params.vfsa_max_iterations = json_params["VFSA-MaxIterations"].toInt();
    }
    if (json_params.contains("VFSA-Parallel")) {
      params.vfsa_parallel = json_params["VFSA-Parallel"].toBool();
    }
    if (json_params.contains("VFSA-InitTemp")) {
      params.vfsa_init_temp = json_params["VFSA-InitTemp"].toDouble();
    }
    if (json_params.contains("VFSA-TempScale")) {
      params.vfsa_temp_scale = json_params["VFSA-TempScale"].toDouble();
    }

    // -----------------------------------------------------
    // SPSA PARAMETERS
    if (json_params.contains("SPSA-MaxIterations")) {
      params.spsa_max_iterations = json_params["SPSA-MaxIterations"].toInt();
    }
    if (json_params.contains("SPSA-gamma")) {
      params.spsa_gamma = json_params["SPSA-gamma"].toDouble();
    }
    if (json_params.contains("SPSA-alpha")) {
      params.spsa_alpha = json_params["SPSA-alpha"].toDouble();
    }
    if (json_params.contains("SPSA-c")) {
      params.spsa_c = json_params["SPSA-c"].toDouble();
    }
    if (json_params.contains("SPSA-A")) {
      params.spsa_A = json_params["SPSA-A"].toDouble();
    } else {
      params.spsa_A = 0.1 * params.spsa_max_iterations;
    }
    if (json_params.contains("SPSA_a")) {
      params.spsa_a = json_params["SPSA_a"].toDouble();
    }
    if (json_params.contains("SPSA-InitStepMagnitude")) {
      params.spsa_init_step_magnitude = json_params["SPSA-InitStepMagnitude"].toDouble();
    }

    // -----------------------------------------------------
    //  HYBRID PARAMETERS
    // Hybrid parameters :: HybridSwitchMode
    if (json_params.contains("HybridSwitchMode")) {
      if (json_params["HybridSwitchMode"].toString() == "OnConvergence") {
        params.hybrid_switch_mode = "OnConvergence";
      } else {
        throw std::runtime_error("HybridSwitchMode " + not_rec);
      }
    }

    // Hybrid parameters :: HybridTerminationCondition
    if (json_params.contains("HybridTerminationCondition")) {
      if (json_params["HybridTerminationCondition"].toString() == "NoImprovement") {
        params.hybrid_termination_condition = "NoImprovement";
      } else {
        throw std::runtime_error("HybridTerminationCondition " + not_rec);
      }
    }

    // Hybrid parameters :: HybridMaxIterations
    if (json_params.contains("HybridMaxIterations")) {
      if (json_params["HybridMaxIterations"].toInt() >= 1) {
        params.hybrid_max_iterations =
            json_params["HybridMaxIterations"].toInt();
      } else {
        throw std::runtime_error(ind_val + "HybridMaxIterations");
      }
    }

    // -----------------------------------------------------
    //  TRUST REGION PARAMETERS
    // Trust Region parameters :: Initial radius
    set_opt_prop_double(params.tr_init_rad, json_params, "TR-InitRad", vp_);

    set_opt_prop_double(params.tr_tol_f, json_params, "TR-TolF", vp_);
    set_opt_prop_double(params.tr_eps_c, json_params, "TR-EpsC", vp_);
    set_opt_prop_double(params.tr_eta_0, json_params, "TR-Eta0", vp_);
    set_opt_prop_double(params.tr_eta_1, json_params, "TR-Eta1", vp_);

    // Thesholds
    set_opt_prop_double(params.tr_piv_thld, json_params, "TR-PivThld", vp_);
    set_opt_prop_double(params.tr_add_thld, json_params, "TR-AddThld", vp_);
    set_opt_prop_double(params.tr_xch_thld, json_params, "TR-xchThld", vp_);

    // Radii
    set_opt_prop_double(params.tr_rad_max, json_params, "TR-RadMax", vp_);
    set_opt_prop_double(params.tr_rad_fac, json_params, "TR-RadFac", vp_);
    set_opt_prop_double(params.tr_rad_tol, json_params, "TR-RadTol", vp_);

    // Gamma factors
    set_opt_prop_double(params.tr_gamma_inc, json_params, "TR-GammaInc", vp_);
    set_opt_prop_double(params.tr_gamma_dec, json_params, "TR-GammaDec", vp_);

    // Criticality
    set_opt_prop_double(params.tr_crit_mu,    json_params, "TR-CritMu", vp_);
    set_opt_prop_double(params.tr_crit_omega, json_params, "TR-CritOmega", vp_);
    set_opt_prop_double(params.tr_crit_beta,  json_params, "TR-CritBeta", vp_);


    // Trust Region parameters :: Lower bound
    if (json_params.contains("TR-LowerBnd")) {
      if (json_params.contains("TR-UpperBnd")) {
        if (json_params["TR-LowerBnd"].toDouble() < json_params["TR-UpperBnd"].toDouble()) {
          params.tr_lower_bnd = json_params["TR-LowerBnd"].toDouble();
        } else { E(ind_val + "TR-LowerBnd", md_, cl_); }
      } else {
        params.tr_lower_bnd = json_params["TR-LowerBnd"].toDouble();
      }
    }

    // Trust Region parameters :: Upper bound
    if (json_params.contains("TR-UpperBnd")) {
      if (json_params.contains("TR-LowerBnd")) {
        if (json_params["TR-LowerBnd"].toDouble() < json_params["TR-UpperBnd"].toDouble()) {
          params.tr_upper_bnd = json_params["TR-UpperBnd"].toDouble();
        } else { E(ind_val + "TR-UpperBnd", md_, cl_); }
      } else {
        params.tr_upper_bnd = json_params["TR-UpperBnd"].toDouble();
      }
    }

    set_opt_prop_int(params.tr_num_init_x, json_params, "TR-NumInitX", vp_);
    set_opt_prop_int(params.tr_iter_max,   json_params, "TR-IterMax", vp_);

    set_opt_prop_string(params.tr_init_smpln, json_params, "TR-InitSmpln", vp_);
    set_opt_prop_string(params.tr_basis, json_params, "TR-Basis", vp_);
    set_opt_prop_string(params.tr_prob_name, json_params, "TR-ProbName", vp_);

    // Trust Region parameters :: RNG seed
    if (json_params.contains("RNGSeed")) {
      params.rng_seed = json_params["RNGSeed"].toInt();
    } else {
      params.rng_seed = std::time(0);
    }

    // -----------------------------------------------------
    //  SQP [SNOPT] PARAMS
    set_opt_prop_double(params.sqp_ftol, json_params, "SQP-FTol", vp_);
    set_opt_prop_double(params.sqp_upper_bnd, json_params, "SQP-UpperBnd", vp_);
    set_opt_prop_double(params.sqp_lower_bnd, json_params, "SQP-LowerBnd", vp_);

    set_opt_prop_double(params.sqp_linesearch_tol, json_params, "SQP-LinesearchTol", vp_);















  } catch (std::exception const &ex) {
    throw UnableToParseOptimizerParametersSectionException(
        "Unable to parse optimizer parameters: "
            + std::string(ex.what()));
  }

  return params;
}

Optimizer::Objective
Optimizer::parseObjective(QJsonObject &json_objective) {

  Objective obj;

  try {
    QString objective_type = json_objective["Type"].toString();

    if (QString::compare(objective_type, "WeightedSum") == 0) {

      obj.type = ObjectiveType::WeightedSum;
      obj.weighted_sum = QList<Objective::WeightedSumComponent>();

      QJsonArray json_components =
          json_objective["WeightedSumComponents"].toArray();

      for (int i = 0; i < json_components.size(); ++i) {

        Objective::WeightedSumComponent component;
        component.coefficient =
            json_components.at(i).toObject()["Coefficient"].toDouble();
        component.property =
            json_components.at(i).toObject()["Property"].toString();

        if (json_components.at(i).toObject()["IsWellProp"].toBool()) {
          component.is_well_prop = true;
          component.well =
              json_components.at(i).toObject()["Well"].toString();
        } else {
          component.is_well_prop = false;
        }

        component.time_step =
            json_components.at(i).toObject()["TimeStep"].toInt();
        obj.weighted_sum.append(component);
      }

    } else if (QString::compare(objective_type, "NPV") == 0) {
      obj.type = ObjectiveType::NPV;
      obj.NPV_sum = QList<Objective::NPVComponent>();

      QJsonArray json_components =
          json_objective["NPVComponents"].toArray();

      for (int i = 0; i < json_components.size(); ++i) {
        Objective::NPVComponent component;
        set_req_prop_double(component.coefficient, json_components[i].toObject(), "Coefficient", vp_);
        set_req_prop_string(component.property, json_components[i].toObject(), "Property", vp_);
        set_opt_prop_double(component.discount, json_components[i].toObject(), "DiscountFactor", vp_);
        set_req_prop_string(component.interval, json_components[i].toObject(), "Interval", vp_);
        set_opt_prop_bool(component.usediscountfactor, json_components[i].toObject(), "UseDiscountFactor", vp_);
        obj.NPV_sum.append(component);
      }

    } else if (QString::compare(objective_type, "Augmented") == 0) {
      stringstream so;

      obj.type = ObjectiveType::Augmented;
      obj.terms = vector<Objective::AugTerms>();

      QJsonArray terms = json_objective["Terms"].toArray();

      for (int ii = 0; ii < terms.size(); ++ii) {
        Objective::AugTerms term;
        auto json_term = terms[ii].toObject();
        set_req_prop_double(term.coefficient, json_term, "Coefficient", vp_);
        set_req_prop_bool(term.active, json_term, "Active", vp_);
        set_req_prop_string(term.prop_name, json_term, "PropName", vp_);

        if (!set_opt_prop_string(term.scaling, json_term, "Scaling", vp_)) {
          term.scaling = "None";
        }

        if (json_term.contains("Wells") && json_term.contains("Segments")) {
          QJsonArray term_wells = json_term["Wells"].toArray();
          for (auto wname : term_wells) {
            term.wells.push_back(wname.toString().toStdString());
          }

          QJsonArray term_segs = json_term["Segments"].toArray();
          for (int ii = 0; ii < term_segs.size(); ++ii) {
            vector<int> segs_vec;
            QJsonArray segs_array = term_segs.at(ii).toArray();
            for (int jj = 0; jj < segs_array.size(); ++jj) {
              segs_vec.push_back(segs_array.at(jj).toInt());
            }
            term.segments[term.wells[ii]] = segs_vec;
          }
        }
        so << term.showTerms();
        obj.terms.push_back(term);
      }
      if (vp_.vSET > 3) { ext_info(so.str(), md_, cl_, vp_.lnw); }

    } else {
      string em = "Objective type " + objective_type.toStdString() + " not recognized";
      throw UnableToParseOptimizerObjectiveSectionException(em);
    }

    if (json_objective.contains("UsePenaltyFunction")) {
      obj.use_penalty_function = json_objective["UsePenaltyFunction"].toBool();
    } else {
      obj.use_penalty_function = false;
    }

    if(json_objective.contains("SeparateHorizontalAndVertical")){
      obj.separatehorizontalandvertical =
          json_objective["SeparateHorizontalAndVertical"].toBool();
    } else {
      obj.separatehorizontalandvertical= false;
    }

    if (json_objective.contains("UseWellCost")) {
      obj.use_well_cost = json_objective["UseWellCost"].toBool();

      if (obj.separatehorizontalandvertical) {
        if (json_objective.contains("WellCostXY")) {
          obj.wellCostXY = json_objective["WellCostXY"].toDouble();
          if (json_objective.contains("WellCostZ")) {
            obj.wellCostZ = json_objective["WellCostZ"].toDouble();
          } else {
            throw UnableToParseOptimizerObjectiveSectionException(
                "Unable to parse optimizer objective a WellCostZ "
                "was not defined, while SeparateHorizontalAndVertical "
                "was invoked");
          }
        } else {
          throw UnableToParseOptimizerObjectiveSectionException(
              "Unable to parse optimizer objective a WellCostXY "
              "was not defined, while SeparateHorizontalAndVertical "
              "was invoked");
        }
      } else {
        if (json_objective.contains("WellCost")) {
          obj.wellCost = json_objective["WellCost"].toDouble();
          obj.wellCostXY = 0;
          obj.wellCostZ = 0;
        } else {
          throw UnableToParseOptimizerObjectiveSectionException(
              "Unable to parse optimizer objective a WellCost "
              "was not defined, while UseWellCost was invoked");
        }
      }

    } else {
      obj.use_well_cost = false;
      obj.wellCost = 0;
      obj.wellCostXY = 0;
      obj.wellCostZ = 0;
    }

  } catch (std::exception const &ex) {
    throw UnableToParseOptimizerObjectiveSectionException(
        "Unable to parse optimizer objective: " + std::string(ex.what()));
  }

  return obj;
}

Optimizer::OptimizerType Optimizer::parseType(QString &type) {
  OptimizerType opt_type;
  if (QString::compare(type, "Compass") == 0) {
    opt_type = OptimizerType::Compass;

  } else if (QString::compare(type, "APPS") == 0) {
    opt_type = OptimizerType::APPS;

  } else if (QString::compare(type, "GeneticAlgorithm") == 0) {
    opt_type = OptimizerType::GeneticAlgorithm;

  } else if (QString::compare(type, "EGO") == 0) {
    opt_type = OptimizerType::EGO;

  } else if (QString::compare(type, "ExhaustiveSearch2DVert") == 0) {
    opt_type = OptimizerType::ExhaustiveSearch2DVert;

  } else if (QString::compare(type, "PSO") == 0) {
    opt_type = OptimizerType::PSO;

  } else if (QString::compare(type, "CMA_ES") == 0) {
    opt_type = OptimizerType::CMA_ES;

  } else if (QString::compare(type, "VFSA") == 0) {
    opt_type = OptimizerType::VFSA;

  } else if (QString::compare(type, "SPSA") == 0) {
    opt_type = OptimizerType::SPSA;

  } else if (QString::compare(type, "Hybrid") == 0) {
    opt_type = OptimizerType::Hybrid;

  } else if (QString::compare(type, "TrustRegionOptimization") == 0) {
    opt_type = OptimizerType::TrustRegionOptimization;

  } else if (QString::compare(type, "DFTR") == 0) {
    opt_type = OptimizerType::DFTR;

  } else if (QString::compare(type, "SQP") == 0) {
    opt_type = OptimizerType::SQP;

  } else if (QString::compare(type, "APPS_DFTR") == 0) {
    opt_type = OptimizerType::APPS_DFTR;

  } else {
    em_ = "Optimizer type " + type.toStdString() + " not recognized.";
    throw OptimizerTypeNotRecognizedException(em_);
  }
  return opt_type;
}

QList<Optimizer::HybridComponent>
Optimizer::parseHybridComponents(QJsonObject &json_optimizer) {

  QList<HybridComponent> comps;
  for (auto json_comp : json_optimizer["HybridComponents"].toArray()) {
    HybridComponent comp;
    QString type = json_comp.toObject()["Type"].toString();
    QJsonObject json_params = json_comp.toObject()["Parameters"].toObject();
    comp.type = parseType(type);
    comp.parameters = parseParameters(json_params);
    comps.push_back(comp);
  }
  return comps;
}

Optimizer::Optimizer(Optimizer::HybridComponent hc) {
  type_ = hc.type;
  parameters_ = hc.parameters;
}

}