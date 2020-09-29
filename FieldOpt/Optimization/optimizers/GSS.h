/***********************************************************
Created by einar on 11/3/16.
Copyright (C) 2015-2019
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 02/19/19 Einar J. M. Baumann
<einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
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

#ifndef FIELDOPT_GSS_H
#define FIELDOPT_GSS_H

#include "Optimization/optimizer.h"
#include <Eigen/Core>
#include <vector>

using namespace Eigen;
using namespace std;

namespace Optimization {
namespace Optimizers {

/*!
 * @brief The _abstract_ GSS class implements a generalized form of the
 * Generating Set Search algorithm presented by Kolda, Torczon and Lewis
 * in the 2003 paper Optimization By Direct Search: New perspectives on
 * some classical and modern methods (DOI: 10.1137/S003614450242889)
 *
 * This implementation also borrows some from the Kolda's 2005 paper on
 * the Asynchronous Parallel Pattern Search (APPS) (DOI: 0.1137/040603589)
 * in that it allows for multiple current step lengths (e.g. one per search
 * direction). This is to allow for more flexibility in the implementation
 * of specific algorithms (e.g. APPS).
 *
 * @note that this is an abstract class. It must be extended by some other
 * class (e.g. like the CompassSearch class) to provide specific patterns
 * and contraction/expansion parameter values.
 */
class GSS : public Optimizer {
 public:
  /*!
   * @brief General constructor for GSS algorithms. Sets the step_tol_ property and calls the primary
   * Optimizer constructor.
   *
   * The following properties must be set in the constructor by classes extending this class:
   *
   *      contr_fac_  : The contraction factor.
   *      expan_fac_  : The expansion factor.
   *      directions_ : The set of search directions to be used.
   *      step_lengths_ : The set of step lengts to be used (one per step direction).
   */
  GSS(Settings::Optimizer *settings,
      Case *base_case,
      Model::Properties::VarPropContainer *variables,
      Reservoir::Grid::Grid *grid,
      Logger *logger,
      CaseHandler *case_handler=0,
      Constraints::ConstraintHandler *constraint_handler=0
  );

  /*!
   * \brief IsFinished Check if the optimization is finished.
   *
   * The algorithm has two termination conditions: max #
   * of objective function evaluations and min step length.
   * \return True if the algorithm has finished, otherwise false.
   */
  TerminationCondition IsFinished();

 protected:
  //!< # of variables in the problem. Used in initialization.
  int num_vars_;

  //!< Step length convergence tolerance.
  VectorXd step_tol_;

  //!< Step length contraction factor.
  double contr_fac_;

  //!< Step length expansion factor.
  double expan_fac_;

  //!< Vector of step lengths.
  VectorXd step_lengths_;

  //!< Vector of search directions.
  vector<VectorXi> directions_;

  /*!
   * @brief Contract the search pattern: step_lengths_ * contr_fac_
   *
   * @param dirs (optional) The direction indices to expand.
   * If not provided, the expansion will be applied to all
   * directions.
   */
  void contract(vector<int> dirs = vector<int>{-1});

  /*!
   * @brief Expand the search pattern: step_lengths_ * expan_fac_
   *
   * @param dirs (optional) The direction indices to expand.
   * If not provided, the expansion will be applied to all
   * directions.
   */
  void expand(vector<int> dirs = vector<int>{-1});

  /*!
   * @brief Set all step lengths _proportionally_ equal
   * to the length provided.
   *
   * This will only work properly if bound constraints are
   * provided for all variables. If not, all directions will
   * be set to the same value (as described in the original
   * paper): max(step_min, len).
   *
   * The modification to the operation described in the paper
   * is intended to account for variables with very different
   * bounds (e.g. when optimizing using the PolarSpline well
   * formulation, or well placement in general in thin
   * reservoirs).
   *
   * Given that the variable at dir_idx has bounds [min_a, max_a],
   * a number S is computed as S = (len - min_a)/(max_a - min_a).
   * The step length at dir_idx is set to len, and the remaining
   * step lengths len_i for each direction i are computed as
   * len_i = max(step_tol_i, min_i + S * (max_i - min_i)).
   *
   * @param dir_idx The direction index that produced an improvement.
   * @param len The step length that produced an improvement.
   */
  void set_step_lengths(int dir_idx, double len);

  /*!
   * @brief Generate a set of trial points.
   * @param dirs (optional) The direction indices in which perturbations
   * should be created.
   *
   * @return A list of new trial points.
   */
  QList<Case *> generate_trial_points(vector<int> dirs = vector<int>{-1});

  /*!
   * @brief Check if the algorithm has converged, i.e. if all current step lengths
   * are below the step length convergence tolerance.
   * @return
   */
  bool is_converged();

  /*!
   * @brief Remove the case that has the worst origin from the evaluation queue.
   * @return Return a pointer to the case that is removed.
   */
  Case *dequeue_case_with_worst_origin();

  string im_ = "", wm_ = "", em_ = "";
  string cl_ = "GSS";
  string md_ = "Optimization::Optimizers";
  Settings::VerbParams vp_;

  VectorXd vd_ = VectorXd::Zero(0);
  string fl_ = "dbg_aug_obj.py";

 private:

  /*!
   * @brief Create a perturbation from a point in the specified direction index.
   * @tparam T An Eigen::VectorX object.
   * @param base The point to perturb from.
   * @param dir The direction index in which to perturb.
   * @return A perturbation (trial point).
   */
  template <typename T>
  Matrix<T, Dynamic, 1> perturb(Matrix<T, Dynamic, 1> base, int dir);

};

}
}
#endif //FIELDOPT_GSS_H
