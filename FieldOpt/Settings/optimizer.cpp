/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb-pcg@gmail.com>

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

Optimizer::Optimizer(QJsonObject json_optimizer, VerbParams vp) {
  vp_=vp;

  // Get the root objects.
  QJsonObject json_parameters = json_optimizer["Parameters"].toObject();
  QJsonObject json_objective = json_optimizer["Objective"].toObject();
  QJsonArray json_constraints = json_optimizer["Constraints"].toArray();
  QString type = json_optimizer["Type"].toString();

  scaling_ = json_optimizer["ScaleVars"].toBool();

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

  optmzr_constraint.scaling_ = scaling_;

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
    string em = "A constraint must always specify either the Well or the Wells property.";
    throw runtime_error(em);
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
      set_opt_prop_double(optmzr_constraint.max, json_constraint, "Max");
    }
    if (json_constraint.contains("Min")) {
      set_opt_prop_double(optmzr_constraint.min, json_constraint, "Min");
    }

    // Packer- and ICV Constraints
  } else if (QString::compare(constraint_type, "ICVConstraint") == 0) {

    optmzr_constraint.type = ConstraintType::ICVConstraint;
    set_opt_prop_double(optmzr_constraint.max, json_constraint, "Max");

    if (optmzr_constraint.max >= 7.8540E-03) {
      string wm = "Maximum valve size is too big. Setting it to 7.8539-3.";
      ext_warn(wm, md_, cl_, vp_.lnw);
      optmzr_constraint.max = 7.8539E-3;
    }
    set_opt_prop_double(optmzr_constraint.min, json_constraint, "Min");

  } else if (CnstrCmp(constraint_type, "PackerConstraint")) {
    optmzr_constraint.type = ConstraintType::PackerConstraint;

  } else if (QString::compare(constraint_type, "Rate") == 0) {
    optmzr_constraint.type = ConstraintType::Rate;
    if (json_constraint.contains("Max"))
      set_opt_prop_double(optmzr_constraint.max, json_constraint, "Max");
    if (json_constraint.contains("Min"))
      set_opt_prop_double(optmzr_constraint.min, json_constraint, "Min");

  } else if (CnstrCmp(constraint_type, "Boundary2D")) {

    optmzr_constraint.type = ConstraintType::PseudoContBoundary2D;
    set_opt_prop_double(optmzr_constraint.box_imin, json_constraint, "Imin");
    set_opt_prop_double(optmzr_constraint.box_imax, json_constraint, "Imax");
    set_opt_prop_double(optmzr_constraint.box_jmin, json_constraint, "Jmin");
    set_opt_prop_double(optmzr_constraint.box_jmax, json_constraint, "Jmax");

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


Optimizer::Parameters Optimizer::parseParameters(QJsonObject &json_parameters) {
  Parameters params;
  string ind_val = "Invalid value for setting";
  string not_rec = "" + not_rec;

  try {

    // -----------------------------------------------------
    // GSS PARAMETERS
    // GSS parameters :: MaxEvaluations
    if (json_parameters.contains("MaxEvaluations"))
      params.max_evaluations = json_parameters["MaxEvaluations"].toInt();

    // GSS parameters :: AutoStepLengths
    if (json_parameters.contains("AutoStepLengths"))
      params.auto_step_lengths = json_parameters["AutoStepLengths"].toBool();

    // GSS parameters :: AutoStepInitScale
    if (json_parameters.contains("AutoStepInitScale"))
      params.auto_step_init_scale = json_parameters["AutoStepInitScale"].toDouble();

    // GSS parameters :: AutoStepConvScale
    if (json_parameters.contains("AutoStepConvScale"))
      params.auto_step_conv_scale = json_parameters["AutoStepConvScale"].toDouble();

    // GSS parameters :: InitialStepLength
    if (json_parameters.contains("InitialStepLength")) {
      set_req_prop_double(params.initial_step_length, json_parameters, "InitialStepLength");
    }

    // GSS parameters :: MinimumStepLength
    if (json_parameters.contains("MinimumStepLength")) {
      set_req_prop_double(params.minimum_step_length, json_parameters, "MinimumStepLength");
    }

    // GSS parameters :: ContractionFactor
    if (json_parameters.contains("ContractionFactor"))
      params.contraction_factor = json_parameters["ContractionFactor"].toDouble();
    else { params.contraction_factor = 0.5; }

    // GSS parameters :: ExpansionFactor
    if (json_parameters.contains("ExpansionFactor"))
      params.expansion_factor = json_parameters["ExpansionFactor"].toDouble();
    else { params.expansion_factor = 1.0; }

    // GSS parameters :: MaxQueueSize
    if (json_parameters.contains("MaxQueueSize"))
      params.max_queue_size = json_parameters["MaxQueueSize"].toDouble();
    else { params.max_queue_size = 2; }

    // GSS parameters :: Pattern
    if (json_parameters.contains("Pattern"))
      params.pattern = json_parameters["Pattern"].toString();
    else params.pattern = "Compass";

    // -----------------------------------------------------
    // GA PARAMETERS
    // GA parameters :: MaxGenerations
    if (json_parameters.contains("MaxGenerations"))
      params.max_generations = json_parameters["MaxGenerations"].toInt();
    else { params.max_generations = 50; }

    // GA parameters :: PopulationSize
    if (json_parameters.contains("PopulationSize"))
      params.population_size = json_parameters["PopulationSize"].toInt();
    else { params.population_size = -1; } // Will be properly set in optimizer.

    // GA parameters :: CrossoverProbability
    if (json_parameters.contains("CrossoverProbability"))
      params.p_crossover = json_parameters["CrossoverProbability"].toDouble();
    else { params.p_crossover = 0.1; }

    // GA parameters :: DiscardParameter
    if (json_parameters.contains("DiscardParameter"))
      params.discard_parameter = json_parameters["DiscardParameter"].toDouble();
    else { params.discard_parameter = -1; } // Will be properly set in optimizer

    // GA parameters :: DecayRate
    if (json_parameters.contains("DecayRate"))
      params.decay_rate = json_parameters["DecayRate"].toDouble();
    else { params.decay_rate = 4.0; }

    // GA parameters :: MutationStrength
    if (json_parameters.contains("MutationStrength"))
      params.mutation_strength = json_parameters["MutationStrength"].toDouble();
    else { params.mutation_strength = 0.25; }

    // GA parameters :: StagnationLimit
    if (json_parameters.contains("StagnationLimit"))
      params.stagnation_limit = json_parameters["StagnationLimit"].toDouble();
    else { params.stagnation_limit = 1e-10; }

    // GA parameters :: LowerBound
    if (json_parameters.contains("LowerBound"))
      params.lower_bound = json_parameters["LowerBound"].toDouble();
    else { params.lower_bound = -10; }

    // GA parameters :: UpperBound
    if (json_parameters.contains("UpperBound"))
      params.upper_bound = json_parameters["UpperBound"].toDouble();
    else { params.upper_bound = 10; }

    // -----------------------------------------------------
    // PSO PARAMETERS
    // PSO parameters :: PSO-LearningFactor1
    if(json_parameters.contains("PSO-LearningFactor1")){
      params.pso_learning_factor_1 = json_parameters["PSO-LearningFactor1"].toDouble();
    } else { params.pso_learning_factor_1 = 2; }

    // PSO parameters :: PSO-LearningFactor2
    if(json_parameters.contains("PSO-LearningFactor2")){
      params.pso_learning_factor_2 = json_parameters["PSO-LearningFactor2"].toDouble();
    } else { params.pso_learning_factor_2 = 2; }

    // PSO parameters :: PSO-SwarmSize
    if(json_parameters.contains("PSO-SwarmSize")){
      params.pso_swarm_size = json_parameters["PSO-SwarmSize"].toDouble();
    } else { params.pso_swarm_size = 50; }

    // PSO parameters :: PSO-VelocityScale
    if(json_parameters.contains("PSO-VelocityScale")){
      params.pso_velocity_scale = json_parameters["PSO-VelocityScale"].toDouble();
    } else { params.pso_velocity_scale = 1.0; }

    // -----------------------------------------------------
    // EGO PARAMEETERS
    // EGO Parameters :: EGO-InitGuesses
    if (json_parameters.contains("EGO-InitGuesses")) {
      params.ego_init_guesses = json_parameters["EGO-InitGuesses"].toInt();
    }

    // EGO Parameters :: EGO-InitSamplingMethod
    if (json_parameters.contains("EGO-InitSamplingMethod")) {
      QString method = json_parameters["EGO-InitSamplingMethod"].toString();
      if (QString::compare(method, "Random") == 0 || QString::compare(method, "Uniform") == 0)
        params.ego_init_sampling_method = json_parameters["EGO-InitSamplingMethod"].toString().toStdString();
      else {
        Printer::error("EGO-InitSamplingMethod " + method.toStdString() + " not recognized.");
        throw std::runtime_error("Failed reading EGO settings.");
      }
    }

    // EGO Parameters :: EGO-Kernel
    if (json_parameters.contains("EGO-Kernel")) {
      QStringList available_kernels = { "CovLinearard", "CovLinearone", "CovMatern3iso",
                                        "CovMatern5iso", "CovNoise", "CovRQiso", "CovSEard",
                                        "CovSEiso", "CovPeriodicMatern3iso", "CovPeriodic"};
      if (available_kernels.contains(json_parameters["EGO-Kernel"].toString())) {
        params.ego_kernel = json_parameters["EGO-Kernel"].toString().toStdString();
      }
      else {
        Printer::error("EGO-Kernel " + json_parameters["EGO-Kernel"].toString().toStdString() + " not recognized.");
        Printer::info("Available kernels: " + available_kernels.join(", ").toStdString());
        throw std::runtime_error("Failed reading EGO settings.");
      }
    }

    // EGO Parameters :: EGO-AF
    if (json_parameters.contains("EGO-AF")) {
      QStringList available_afs = { "ExpectedImprovement", "ProbabilityOfImprovement" };
      if (available_afs.contains(json_parameters["EGO-AF"].toString())) {
        params.ego_af = json_parameters["EGO-AF"].toString().toStdString();

      } else {
        Printer::error("EGO-AF " + json_parameters["EGO-AF"].toString().toStdString() + " not recognized.");
        Printer::info("Available acquisition functions: " + available_afs.join(", ").toStdString());
        throw std::runtime_error("Failed reading EGO settings.");
      }
    }

    // -----------------------------------------------------
    // CMA-ES PARAMETERS
    if (json_parameters.contains("ImproveBaseCase")) {
      params.improve_base_case = json_parameters["ImproveBaseCase"].toBool();
    }

    // -----------------------------------------------------
    // VFSA PARAMETERS
    if (json_parameters.contains("VFSA-EvalsPrIteration")) {
      params.vfsa_evals_pr_iteration = json_parameters["VFSA-EvalsPrIteration"].toInt();
    }
    if (json_parameters.contains("VFSA-MaxIterations")) {
      params.vfsa_max_iterations = json_parameters["VFSA-MaxIterations"].toInt();
    }
    if (json_parameters.contains("VFSA-Parallel")) {
      params.vfsa_parallel = json_parameters["VFSA-Parallel"].toBool();
    }
    if (json_parameters.contains("VFSA-InitTemp")) {
      params.vfsa_init_temp = json_parameters["VFSA-InitTemp"].toDouble();
    }
    if (json_parameters.contains("VFSA-TempScale")) {
      params.vfsa_temp_scale = json_parameters["VFSA-TempScale"].toDouble();
    }

    // -----------------------------------------------------
    // SPSA PARAMETERS
    if (json_parameters.contains("SPSA-MaxIterations")) {
      params.spsa_max_iterations = json_parameters["SPSA-MaxIterations"].toInt();
    }
    if (json_parameters.contains("SPSA-gamma")) {
      params.spsa_gamma = json_parameters["SPSA-gamma"].toDouble();
    }
    if (json_parameters.contains("SPSA-alpha")) {
      params.spsa_alpha = json_parameters["SPSA-alpha"].toDouble();
    }
    if (json_parameters.contains("SPSA-c")) {
      params.spsa_c = json_parameters["SPSA-c"].toDouble();
    }
    if (json_parameters.contains("SPSA-A")) {
      params.spsa_A = json_parameters["SPSA-A"].toDouble();
    } else params.spsa_A = 0.1 * params.spsa_max_iterations;
    if (json_parameters.contains("SPSA_a")) {
      params.spsa_a = json_parameters["SPSA_a"].toDouble();
    }
    if (json_parameters.contains("SPSA-InitStepMagnitude")) {
      params.spsa_init_step_magnitude = json_parameters["SPSA-InitStepMagnitude"].toDouble();
    }

    // -----------------------------------------------------
    //  HYBRID PARAMETERS
    // Hybrid parameters :: HybridSwitchMode
    if (json_parameters.contains("HybridSwitchMode")) {
      if (json_parameters["HybridSwitchMode"].toString() == "OnConvergence") {
        params.hybrid_switch_mode = "OnConvergence";
      } else {
        throw std::runtime_error("HybridSwitchMode " + not_rec);
      }
    }

    // Hybrid parameters :: HybridTerminationCondition
    if (json_parameters.contains("HybridTerminationCondition")) {
      if (json_parameters["HybridTerminationCondition"].toString() == "NoImprovement") {
        params.hybrid_termination_condition = "NoImprovement";
      } else {
        throw std::runtime_error("HybridTerminationCondition " + not_rec);
      }
    }

    // Hybrid parameters :: HybridMaxIterations
    if (json_parameters.contains("HybridMaxIterations")) {
      if (json_parameters["HybridMaxIterations"].toInt() >= 1) {
        params.hybrid_max_iterations =
            json_parameters["HybridMaxIterations"].toInt();
      } else {
        throw std::runtime_error(ind_val + "HybridMaxIterations");
      }
    }

    // -----------------------------------------------------
    //  TRUST REGION PARAMETERS
    // Trust Region parameters :: Initial radius
    if (json_parameters.contains("InitialTrustRegionRadius")) {
      if (json_parameters["InitialTrustRegionRadius"].toDouble() >= 0.0) {
        params.tr_initial_radius =
            json_parameters["InitialTrustRegionRadius"].toDouble();
      } else {
        throw std::runtime_error(ind_val + "InitialTrustRegionRadius");
      }
    }

    // Trust Region parameters:: Trust region radius tolerance
    if (json_parameters.contains("TrustRegionRadiusTolerance")) {
      if (json_parameters["TrustRegionRadiusTolerance"].toDouble() >= 0.0) {
        params.tr_tol_radius = json_parameters["TrustRegionRadiusTolerance"].toDouble();
      } else {
        throw std::runtime_error(ind_val + "TrustRegionRadiusTolerance");
      }
    }

    // Trust Region parameters :: Max trust region radius
    if (json_parameters.contains("MaxTrustRegionRadius")) {
      if (json_parameters["MaxTrustRegionRadius"].toDouble() >= 0.0) {
        params.tr_radius_max =
            json_parameters["MaxTrustRegionRadius"].toDouble();
      } else {
        throw std::runtime_error(ind_val + "MaxTrustRegionRadius");
      }
    }

    // Trust Region parameters :: Lower bound
    if (json_parameters.contains("TrustRegionLowerBound")) {
      if (json_parameters.contains("TrustRegionUpperBound")) {
        if (json_parameters["TrustRegionLowerBound"].toDouble() <
            json_parameters["TrustRegionUpperBound"].toDouble()) {
          params.tr_lower_bound = json_parameters["TrustRegionLowerBound"].toDouble();
        } else {
          throw std::runtime_error(ind_val + "TrustRegionLowerBound");
        }
      } else {
        params.tr_lower_bound =
            json_parameters["TrustRegionLowerBound"].toDouble();
      }
    }

    // Trust Region parameters :: Upper bound
    if (json_parameters.contains("TrustRegionUpperBound")) {
      if (json_parameters.contains("TrustRegionLowerBound")) {
        if (json_parameters["TrustRegionLowerBound"].toDouble() <
            json_parameters["TrustRegionUpperBound"].toDouble()) {
          params.tr_upper_bound =
              json_parameters["TrustRegionUpperBound"].toDouble();
        } else {
          throw std::runtime_error(ind_val + "TrustRegionUpperBound");
        }
      } else {
        params.tr_upper_bound =
            json_parameters["TrustRegionUpperBound"].toDouble();
      }
    }

    // Trust Region parameters :: RNG seed
    if (json_parameters.contains("RNGSeed")) {
      params.rng_seed = json_parameters["RNGSeed"].toInt();
    } else {
      params.rng_seed = std::time(0);
    }

    // Trust Region parameters :: Max iter
    if (json_parameters.contains("TRMaxIter")) {
      params.tr_iter_max = json_parameters["TRMaxIter"].toInt();
    } else {
      params.tr_iter_max = 10000;
    }

    // Trust Region parameters :: Criticality Mu
    if (json_parameters.contains("CriticalityMu")) {
      if (json_parameters["CriticalityMu"].toDouble() >= 0.0) {
        params.tr_criticality_mu =
            json_parameters["CriticalityMu"].toDouble();
      } else {
        throw std::runtime_error(ind_val + "CriticalityMu");
      }
    }

    // Trust Region parameters :: Criticality Omega
    if (json_parameters.contains("CriticalityOmega")) {
      if (json_parameters["CriticalityOmega"].toDouble() >= 0.0) {
        params.tr_criticality_omega =
            json_parameters["CriticalityOmega"].toDouble();
      } else {
        throw std::runtime_error(ind_val + "CriticalityOmega");
      }
    }

    // Trust Region parameters :: Criticality Beta
    if (json_parameters.contains("CriticalityBeta")) {
      if (json_parameters["CriticalityBeta"].toDouble() >= 0.0) {
        params.tr_criticality_beta =
            json_parameters["CriticalityBeta"].toDouble();
      } else {
        throw std::runtime_error(ind_val + "CriticalityBeta");
      }
    }

    // Trust Region parameters :: Criticality problem name
    if (json_parameters.contains("ProblemName")) {
      if (json_parameters["ProblemName"].toDouble() >= 0.0) {
        params.tr_prob_name =
            json_parameters["ProblemName"].toString().toStdString();
      } else {
        throw std::runtime_error(ind_val + "ProblemName");
      }
    }

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
        set_req_prop_double(component.coefficient, json_components[i].toObject(), "Coefficient");
        set_req_prop_string(component.property, json_components[i].toObject(), "Property");
        set_opt_prop_double(component.discount, json_components[i].toObject(), "DiscountFactor");
        set_req_prop_string(component.interval, json_components[i].toObject(), "Interval");
        set_opt_prop_bool(component.usediscountfactor, json_components[i].toObject(), "UseDiscountFactor");
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
        set_req_prop_double(term.coefficient, json_term, "Coefficient");
        set_req_prop_bool(term.active, json_term, "Active");
        set_req_prop_string(term.prop_name, json_term, "PropName");

        if (!set_opt_prop_string(term.scaling, json_term, "Scaling")) {
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
      ext_info(so.str(), md_, cl_);

    } else {
      throw UnableToParseOptimizerObjectiveSectionException(
          "Objective type " + objective_type.toStdString()
              + " not recognized");
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
  } else {
    throw OptimizerTypeNotRecognizedException(
        "The optimizer type " + type.toStdString()
            + " was not recognized.");
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