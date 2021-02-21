/***********************************************************
Created by einar on 11/21/16.
Copyright (C) 2019
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 02/19/19
Einar J. M. Baumann <einar.baumann@gmail.com>

Modified 2017-2021 Mathias Bellout
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

#ifndef FIELDOPT_APPS_H
#define FIELDOPT_APPS_H

#include <set>
#include "GSS.h"

namespace Optimization {
namespace Optimizers {

class APPS : public GSS {
 public:
  APPS(Settings::Optimizer *settings, Case *base_case,
       Model::Properties::VarPropContainer *variables,
       Reservoir::Grid::Grid *grid,
       Logger *logger,
       CaseHandler *case_handler=nullptr,
       Constraints::ConstraintHandler *constraint_handler=nullptr
  );

 protected:
  void handleEvaluatedCase(Case *c) override;

  void iterate() override;

 private:
  //!< Maximum length of queue.
  int max_queue_length_;

  //!< Set containing the indices of all active search dirs.
  set<int> active_;

  //!< Mark the direction indices in the vector as active.
  void set_active(vector<int> dirs);

  //!< Mark the direction indices in the vector as inactive.
  void set_inactive(vector<int> dirs);

  //!< Reset the list of active search directions.
  void reset_active();

  //!< Get vector containing all _inactive_ search
  //!< directions with step length greater than step_tol_.
  vector<int> inactive();

  /*!
   * @brief Handle a successful iteration.
   *
   * Will be called by handleEvaluatedCase() if the most
   * recently evaluated case is an improvement on the
   * tentative_best_case_.
   * @param c Most recently evaluated case.
   */
  void successful_iteration(Case *c);

  /*!
   * @brief Handle an unsuccessful iteration.
   *
   * Will be called by handleEvaluatedCase() if the most
   * recently evaluated case is _not_ and improvement on
   * the tentative_best_case_.
   * @param c Most recently evaluated case.
   */
  void unsuccessful_iteration(Case *c);

  /*!
   * @brief Prune the evaluation queue to max_queue_length_.
   */
  void prune_queue();

  /*!
   * @brief Print the state of the optimizer.
   * Detail level depends on the verbosity setting.
   */
  void print_state(string header);
};

}
}


#endif //FIELDOPT_APPS_H
