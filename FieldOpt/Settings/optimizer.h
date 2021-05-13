/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
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

#ifndef SETTINGS_OPTIMIZER_H
#define SETTINGS_OPTIMIZER_H

#include "settings.h"

#include <QList>
#include <QString>
#include <QStringList>

namespace Settings {

using std::vector;
using std::string;
using std::map;
using Printer::num2str;
using Printer::pad_text;
using Printer::E;

/*!
 * \brief The Optimizer class contains optimizer-specific
 * settings. Optimizer settings objects may _only_ be
 * created by the Settings class. They are created when
 * reading a JSON-formatted "driver file".
 */
class Optimizer
{
  friend class Settings;

 public:
  Optimizer(){}
  Optimizer(QJsonObject json_optimizer, VerbParams vp);

  Paths* paths() { return paths_; };
  void setPaths(Paths *p) { paths_ = p; };

  enum OptimizerType {
    Compass, APPS, ExhaustiveSearch2DVert, GeneticAlgorithm,
    EGO, PSO, VFSA, SPSA, CMA_ES, Hybrid, TrustRegionOptimization,
    DFTR, SQP
  };

  enum OptimizerMode { Maximize, Minimize };

  enum ConstraintType { BHP, Rate, SplinePoints,
    WSplineLength, WSplineInterwDist, WellSplineDomain,
    MxWSplineLengthInterwDist, MxWSplineLengthInterwDistResBound,
    ReservoirBoundary, PseudoContBoundary2D, PolarXYZBoundary,
    WellXYZBox,
    ReservoirXYZBoundary, ReservoirBoundaryToe,
    PackerConstraint, ICVConstraint, PolarWellLength,
    PolarAzimuth, PolarElevation, PolarSplineBoundary
  };

  enum ConstraintWellSplinePointsType { MaxMin, Function};
  enum ObjectiveType { WeightedSum, NPV, Augmented};

  struct Parameters {
    // Common parameters
    //!< Max # of evals allowed b/f terminating optimization run
    int max_evaluations;

    //!< Seed to be used for random number generators
    int rng_seed;

    // GSS parameters --------------------------------------
    //!< The initial step length in the algorithm when applicable.
    double initial_step_length;

    //!< The minimum step length in the algorithm when applicable.
    double minimum_step_length;

    //!< The contraction factor for GSS algorithms.
    double contraction_factor;

    //!< The expansion factor for GSS algorithms.
    double expansion_factor;

    //!< Maximum size of evaluation queue.
    int max_queue_size;

    //!< Automatically determine appropriate step
    //!< lengths from bound constraints.
    bool auto_step_lengths = false;

    //!< Scaling factor for auto-determined initial
    //!< step lengths (e.g. 0.25*(upper-lower)
    double auto_step_init_scale = 0.25;

    //!< Scaling factor for auto-determined convergence
    //!< step lengths (e.g. 0.01*(upper-lower)
    double auto_step_conv_scale = 0.01;

    //!< Pattern to be used for GSS algorithms.
    QString pattern;

    // GA parameters ---------------------------------------
    //!< Max iterations. Default: 50
    int max_generations;

    //!< Optional. Can be determined automatically.
    //!< Default: min(10*nvars, 100).
    int population_size;

    //!< Fraction to be discarded during selection.
    //!< Defaults: 1/population.
    double discard_parameter;

    //!< Crossover probability. Default: 0.1.
    double p_crossover;

    //!< Decay rate. Default: 4.0.
    double decay_rate;

    //!< Mutation strength. Default: 0.25.
    double mutation_strength;

    //!< Stagnation limit. Default: 1e-10.
    double stagnation_limit;

    //!< Lower bound. Applied to _all_ variables. Default: -10.0.
    double lower_bound;

    //!< Upper bound. Applied to _all_ variables. Default: +10.0.
    double upper_bound;

    // PSO parameters --------------------------------------
    //!< Learning factor (c1), from the swarms
    //!< best known perturbation. Default: 2
    double pso_learning_factor_1;

    //!< Learning factor (c2), from the individual
    //!< particle's best known perturbation. Default: 2
    double pso_learning_factor_2;

    //!< # of particles in the swarm. Default: 50
    double pso_swarm_size;

    //!< Scaling factor for particle velocities. Default: 1.0
    double pso_velocity_scale;

    // EGO Parameters --------------------------------------
    //!< # of init guesses to be made (default
    //!< is two times the number of variables).
    int ego_init_guesses = -1;

    //!< Sampling method to be used for initial guesses
    //!< (Random or Uniform)
    std::string ego_init_sampling_method = "Random";

    //!< Which kernel function to use for gaussian process model.
    std::string ego_kernel = "CovMatern5iso";

    //!< Which acquisiton function to use.
    std::string ego_af = "ExpectedImprovement";

    // Trust Region parameters -----------------------------
    double tr_init_rad = 1; // Initial TR radius
    // Tols
    double tr_tol_f = 1e-6;
    double tr_eps_c = 1e-5;
    double tr_eta_0 = 0;
    double tr_eta_1 = 0.05;
    // Thesholds
    double tr_piv_thld = 0.0625;
    double tr_add_thld = 100;
    double tr_xch_thld = 1000;
    // Radii
    double tr_rad_max = 1e3;
    double tr_rad_fac = 6;
    double tr_rad_tol = 1e-5;
    // Gamma factors
    double tr_gamma_inc = 2;
    double tr_gamma_dec = 0.5;
    // Criticality
    double tr_crit_mu = 100;
    double tr_crit_omega = 0.5;
    double tr_crit_beta = 10;
    // Bounds
    double tr_lower_bnd = -std::numeric_limits<double>::infinity();
    double tr_upper_bnd = std::numeric_limits<double>::infinity();

    //!< # initial guesses provided to build the TR (def.=1)
    int tr_num_init_x = -1;
    int tr_iter_max = 10000;

    //!< Sampling method for initial guesses (Random or Uniform)
    std::string tr_init_smpln = "Random";
    std::string tr_basis = "diagonalHessian";

    std::string tr_prob_name = "prob0"; // dbg

    // SQP [SNOPT] parameters ------------------------------
    double sqp_ftol = 1e-6;
    double sqp_upper_bnd = 1.0;
    double sqp_lower_bnd = -1.0;

    double sqp_linesearch_tol = 0.9;
    double sqp_major_feasi_tol = 1e-12;
    double sqp_major_iter_lim = 100;
    double sqp_optimality_tol = 1e-12;













    // VFSA Parameters -------------------------------------
    //!< Number of evaluations to be performed pr. iteration (temperature). Default: 1.
    int vfsa_evals_pr_iteration = 1;

    //!< Maximum number of iterations to be performed. Default: 50.
    int vfsa_max_iterations = 50;

    //!< Run generate evals_pr_iteration cases immedeately in each generation? Default: false.
    bool vfsa_parallel = false;

    //!< Initial temperature (same used for all dimensions). Default: 1.0.
    double vfsa_init_temp = 1.0;

    //!< Constant used in scaling temperature. Default: 1.0.
    double vfsa_temp_scale = 1.0;

    // CMA-ES Parameters -----------------------------------
    bool improve_base_case = false;

    // SPSA Parameters -------------------------------------
    //!< # number of iterations to be performed. Default: 50.
    int spsa_max_iterations = 50;

    //!< Affects step length. Recommended 0.602 <= alpha <= 1.0.
    //!< Default: 0.602.
    double spsa_alpha = 0.602;

    //!< Affects step length. Recommended 0.101 <= gamma <= 1/6.
    //!< Default: 0.101.
    double spsa_gamma = 0.101;

    //!< Used to account for noise. Default: 1.0.
    double spsa_c = 1.0;

    //!< Affects step length in early iterations.
    //!< Default: 10% of max iterations.
    double spsa_A = 5;

    //!< Affects step lengths. Default and recommended:
    //!< automatically compute from spsa_init_step_magnitude.
    double spsa_a = 0.0;

    //!< Smallest desired step magnitude in early iterations.
    double spsa_init_step_magnitude = 0.0;

    // Hybrid parameters -----------------------------------
    /*!
     * @brief How switching b/e component optimizers is handled.
     *
     * Default: OnFinished -- switch between components when
     * IsFinished() == true
     *
     * Example: "Optimizer": { "Type": "Hybrid", "Parameters": { "HybridSwitchMode": "OnConvergence" } }
     */
    std::string hybrid_switch_mode = "OnConvergence";

    /*!
     * @brief Termination condition for hybrid optimizer.
     *
     * Default: NoImprovement -- terminate when a component has not managed to
     *                           improve upon the result of the previous component.
     *
     * Example: "Optimizer": { "Type": "Hybrid", "Parameters": { "HybridTerminationCondition": "NoImprovement" } }
     */
    std::string hybrid_termination_condition = "NoImprovement";

    /*!
     * @brief Max iterations for the hybrid optimizer.
     *
     * One iteration implies running each component to completion once.
     * If
     *
     * Default: 2 -- Each component will be executed twice (unleass the HybridTerminationCondition
     *               is met first).
     *
     * Example: "Optimizer": { "Type": "Hybrid", "Parameters": { "HybridMaxIterations": 2 } }
     */
    int hybrid_max_iterations = 2;
  };

  struct Objective {

    //!< Objective function type (e.g. WeightedSum, NPV)
    ObjectiveType type;

    //!< Whether or not to use penalty function (default: false).
    bool use_penalty_function;

    //!<Whether or not to use costs associated to
    //!< wells in calculation of the objective.
    bool use_well_cost;

    //!<Whether or not to use different values
    //!< in the horizontal or vertical direction
    bool separatehorizontalandvertical;

    //!<Cost associated with drilling in the horizontal plane [$/m]
    double wellCostXY;

    //!<Cost associated with drilling in the vertical plane [$/m]
    double wellCostZ;

    //!<Cost associated with drilling the well,
    //!< independent of direction [$/m]
    double wellCost;

    //!< Weighted sum component
    struct WeightedSumComponent {
      double coefficient;
      QString property;
      int time_step;
      bool is_well_prop;
      QString well;
    };

    //!< NPV component
    struct NPVComponent{
      double coefficient;
      string property;
      string interval = "";
      bool usediscountfactor = false;
      string well;
      double discount = 0.0;
    };

    //!< Weighted sum formulation
    QList<WeightedSumComponent> weighted_sum;

    //!< NPV formulation
    QList<NPVComponent> NPV_sum;

    //!< Augmented function term
    struct AugTerms {
      string prop_name; // -> defines prop_type
      double coefficient;
      string prop_spec = "";
      vector<string> wells;
      map<string, vector<int>> segments;
      bool active;
      string scaling;

      string showTerms() {
        stringstream ss;
        string pn, pc, ps;
        int pad = 30;

        pn = "prop_name: " + prop_name;
        pad_text(pn, pad);
        ss << pn + " ";

        pc = "coefficient:" + num2str(coefficient, 4, 0, 10);
        pad_text(pc, pad-4);
        ss << pc + " ";

        ss << "active: " << to_string(active) + "    ";

        ps = "prop_spec: " + prop_spec;

        // if (prop_spec != "") {
        for(string w : wells) {
          ps += segments[w].size() + ", well: " + w + " w/ segs: [ ";
          for  (int ii=0; ii < segments[w].size(); ++ii) {
            ps += num2str(segments[w][ii], 0) + " ";
          }
          ps += "]";
        }
        // } else {
        // ss += string(10, ' ');
        // }

        pad_text(ps, pad+30);
        ss << ps <<  "|";

        return ss.str();
      }
    };

    //!< Augmented formulation
    vector<AugTerms> terms;

  };

  struct Constraint {

    //!< Used to express (x,y,z) coordinates.
    struct RealCoordinate { double x; double y; double z; };

    //!< Used to define a box-shaped 3D area. Max, min each define a corner.
    struct RealMaxMinLimit { RealCoordinate max; RealCoordinate min; };

    //!< The constraint type (e.g. BHP or SplinePoints positions).
    ConstraintType type;

    //!< The name of the well this Constraint applies to.
    QString well;

    //!< List of well names if the constraint applies to more than one.
    QStringList wells;

    //!< Max limit when using constraints like BHP, Rate, ICD
    double max;

    //!< Min limit when using constraints like BHP, Rate, ICD.
    double min;

    //!< Min max limits for geometric ijk box constraints.
    double box_imin, box_imax, box_jmin;
    double box_jmax, box_kmin, box_kmax;

    //!< Min max limits for geometric xyz box constraints.
    double box_xyz_xmin, box_xyz_ymin, box_xyz_zmin;
    double box_xyz_xmax, box_xyz_ymax, box_xyz_zmax;

    double min_length, max_length;
    double min_distance;
    double min_md, max_md;

    //!< The weight to be used when considering the
    //!< constraint in a penalty function. (default: 0.0)
    long double penalty_weight;
    int max_iterations;

    //!< How the SplinePoints constraint is given
    //!< when SplinePoints constraint type is selected.
    ConstraintWellSplinePointsType spline_points_type;

    //!< Box limits a spline point needs to be within to be
    //!< valid when SplinePoints constraint type is selected.
    QList<RealMaxMinLimit> spline_points_limits;

    bool scaling_ = true;

  };

  struct HybridComponent {
    OptimizerType type;
    Parameters parameters;
  };

  //!< Create a basic Optimizer Settings object from a HybridComponent object.
  Optimizer(HybridComponent hc);

  //!< Get the Optimizer type (e.g. Compass).
  OptimizerType type() const { return type_; }

  //!< Get the optimizer mode (maximize/minimize).
  OptimizerMode mode() const { return mode_; }

  //!< Set the optimizer mode (used by HybridOptimizer)
  void set_mode(const OptimizerMode mode) { mode_ = mode; }

  //!< Get the optimizer parameters.
  Parameters parameters() const { return parameters_; }

  //!< Get the optimizer objective function.
  Objective objective() const { return objective_; }

  //!< Get the optimizer constraints.
  QList<Constraint> constraints() const { return constraints_; }

  // Get the list of hybrid-optimizer components when using the HYBRID type.
  QList<HybridComponent> HybridComponents() { return hybrid_components_; }

  //!< Change the RNG seed (used by HybridOptimizer).
  void SetRngSeed(const int seed) { parameters_.rng_seed = seed; }

  void setTRProbName(std::string pn) { parameters_.tr_prob_name = pn; }
  VerbParams verbParams() { return vp_; };

  int scaleVars() { return scale_vars_; }
  double ScaleObjf() { return scale_objf_; }

 private:
  QList<Constraint> constraints_;
  OptimizerType type_;
  Parameters parameters_;
  Objective objective_;
  Paths *paths_;

  int scale_vars_ = 1;
  double scale_objf_ = 1.0;

  string im_, wm_, em_;
  string md_ = "Settings";
  string cl_ = "Optimizer";
  VerbParams vp_;

  //!< Optimization mode (maximize or minimize). Def: Maximize
  OptimizerMode mode_ = OptimizerMode::Maximize;
  QList<HybridComponent> hybrid_components_;

  OptimizerType parseType(QString &type);
  Constraint parseSingleConstraint(QJsonObject json_constraint);
  OptimizerMode parseMode(QJsonObject &json_optimizer);
  Parameters parseParameters(QJsonObject &json_params);
  Objective parseObjective(QJsonObject &json_objective);
  QList<HybridComponent> parseHybridComponents(QJsonObject &json_optimizer);

  bool CnstrCmp(QString ctype, QString con) {
    return (QString::compare(ctype, con) == 0);
  }


};

}


#endif // SETTINGS_OPTIMIZER_H
