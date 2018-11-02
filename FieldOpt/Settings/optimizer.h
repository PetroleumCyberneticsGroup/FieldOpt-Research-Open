/******************************************************************************
   Copyright (C) 2015-2017 Einar J.M. Baumann <einar.baumann@gmail.com>

   This file is part of the FieldOpt project.

   FieldOpt is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   FieldOpt is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FieldOpt.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef SETTINGS_OPTIMIZER_H
#define SETTINGS_OPTIMIZER_H

#include "settings.h"

#include <QList>
#include <QString>
#include <QStringList>

namespace Settings {

/*!
 * \brief The Optimizer class contains optimizer-specific settings. Optimizer settings objects
 * may _only_ be created by the Settings class. They are created when reading a
 * JSON-formatted "driver file".
 */
class Optimizer
{
  friend class Settings;

 public:
  Optimizer(){}
  Optimizer(QJsonObject json_optimizer);
  enum OptimizerType { Compass, APPS, ExhaustiveSearch2DVert, GeneticAlgorithm, EGO, Hybrid };
  enum OptimizerMode { Maximize, Minimize };
  enum ConstraintType { BHP, Rate, SplinePoints,
    WellSplineLength, WellSplineInterwellDistance, WellSplineDomain,
    CombinedWellSplineLengthInterwellDistance,
    CombinedWellSplineLengthInterwellDistanceReservoirBoundary,
    ReservoirBoundary, PseudoContBoundary2D,
    PackerConstraint, ICVConstraint
  };
  enum ConstraintWellSplinePointsType { MaxMin, Function};
  enum ObjectiveType { WeightedSum, NPV};

  struct Parameters {
    // Common parameters
    int max_evaluations; //!< Maximum number of evaluations allowed before terminating the optimization run.
    int rng_seed;        //!< Seed to be used for random number renerators in relevant algorithms.

    // GSS parameters
    double initial_step_length; //!< The initial step length in the algorithm when applicable.
    double minimum_step_length; //!< The minimum step length in the algorithm when applicable.
    double contraction_factor;  //!< The contraction factor for GSS algorithms.
    double expansion_factor;    //!< The expansion factor for GSS algorithms.
    double max_queue_size;      //!< Maximum size of evaluation queue.
    bool auto_step_lengths = false;     //!< Automatically determine appropriate step lengths from bound constraints.
    double auto_step_init_scale = 0.25; //!< Scaling factor for auto-determined initial step lengths (e.g. 0.25*(upper-lower)
    double auto_step_conv_scale = 0.01; //!< Scaling factor for auto-determined convergence step lengths (e.g. 0.01*(upper-lower)
    QString pattern;                     //!< The pattern to be used for GSS algorithms.

    // GA parameters
    int max_generations;      //!< Max iterations. Default: 50
    int population_size;      //!< Optional. Can be determined automatically. Default: min(10*nvars, 100).
    double discard_parameter; //!< Fraction to be discarded during selection. Defaults: 1/population.
    double p_crossover;       //!< Crossover probability. Default: 0.1.
    double decay_rate;        //!< Decay rate. Default: 4.0.
    double mutation_strength; //!< Mutation strength. Default: 0.25.
    double stagnation_limit;  //!< Stagnation limit. Default: 1e-10.
    double lower_bound;       //!< Simple lower bound. This is applied to _all_ variables. Default: -10.0.
    double upper_bound;       //!< Simple upper bound. This is applied to _all_ variables. Default: +10.0.

    // EGO Parameters
    int ego_init_guesses = -1; //!< Number of initial guesses to be made (default is two times the number of variables).
    std::string ego_init_sampling_method = "Random"; //!< Sampling method to be used for initial guesses (Random or Uniform)
    std::string ego_kernel = "CovMatern5iso";        //!< Which kernel function to use for the gaussian process model.
    std::string ego_af = "ExpectedImprovement";      //!< Which acquisiton function to use.
  };

  struct Objective {
    ObjectiveType type; //!< The objective definition type (e.g. WeightedSum, NPV)
    bool use_penalty_function; //!< Whether or not to use penalty function (default: false).
    bool use_well_cost; //!<Whether or not to use costs associated to wells in calculation of the objective.
    bool separatehorizontalandvertical; //!<Whether or not to use different values in the horizontal or vertical direction
    double wellCostXY; //!<Cost associated with drilling in the horizontal plane [$/m]
    double wellCostZ; //!<Cost associated with drilling in the vertical plane [$/m]
    double wellCost; //!<Cost associated with drilling the well, independent of direction [$/m]
    struct WeightedSumComponent {
      double coefficient; QString property; int time_step;
      bool is_well_prop; QString well; }; //!< A component of a weighted sum formulatied objective function
    struct NPVComponent{
      double coefficient; QString property; QString interval;
      bool usediscountfactor; QString well; double discount; };
    QList<WeightedSumComponent> weighted_sum; //!< The expression for the Objective function formulated as a weighted sum
    QList<NPVComponent> NPV_sum;

  };

  struct Constraint {
    struct RealCoordinate { double x; double y; double z; }; //!< Used to express (x,y,z) coordinates.
    struct RealMaxMinLimit { RealCoordinate max; RealCoordinate min; }; //!< Used to define a box-shaped 3D area. Max and min each define a corner.
    ConstraintType type; //!< The constraint type (e.g. BHP or SplinePoints positions).
    QString well; //!< The name of the well this Constraint applies to.
    QStringList wells; //!< List of well names if the constraint applies to more than one.
    double max; //!< Max limit when using constraints like BHP.
    double min; //!< Min limit when using constraints like BHP.
    double box_imin, box_imax, box_jmin, box_jmax, box_kmin, box_kmax; //!< Min max limits for geometrix box constraints.
    double min_length, max_length;
    double min_distance;
    double min_md, max_md;
    long double penalty_weight; //!< The weight to be used when considering the constraint in a penalty function. (default: 0.0)
    int max_iterations;
    ConstraintWellSplinePointsType spline_points_type; //!< How the SplinePoints constraint is given when SplinePoints constraint type is selected.
    QList<RealMaxMinLimit> spline_points_limits; //!< Box limits a spline point needs to be within to be valid when SplinePoints constraint type is selected.
  };

  struct HybridComponent {
    OptimizerType type;
    Parameters parameters;
  };
  Optimizer(HybridComponent hc); //!< Create a basic Optimizer Settings object from a HybridComponent object.

  OptimizerType type() const { return type_; } //!< Get the Optimizer type (e.g. Compass).
  OptimizerMode mode() const { return mode_; } //!< Get the optimizer mode (maximize/minimize).
  Parameters parameters() const { return parameters_; } //!< Get the optimizer parameters.
  Objective objective() const { return objective_; } //!< Get the optimizer objective function.
  QList<Constraint> constraints() const { return constraints_; } //!< Get the optimizer constraints.
  QList<HybridComponent> HybridComponents() { return hybrid_components_; } // Get the list of hybrid-optimizer components when using the HYBRID type.


 private:
  QList<Constraint> constraints_;
  OptimizerType type_;
  Parameters parameters_;
  Objective objective_;
  OptimizerMode mode_;
  QList<HybridComponent> hybrid_components_;

  OptimizerType parseType(QString &type);
  Constraint parseSingleConstraint(QJsonObject json_constraint);
  OptimizerMode parseMode(QJsonObject &json_optimizer);
  Parameters parseParameters(QJsonObject &json_parameters);
  Objective parseObjective(QJsonObject &json_objective);
  QList<HybridComponent> parseHybridComponents(QJsonObject &json_optimizer);
};

}


#endif // SETTINGS_OPTIMIZER_H
