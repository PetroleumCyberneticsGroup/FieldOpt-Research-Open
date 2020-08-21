/***********************************************************
Created by bellout on 8/19/20.

Copyright (C) 2020
Mathias Bellout <chakibbb-pcg@gmail.com>

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

#ifndef FIELDOPT_OPTIMIZATION_OBJECTIVE_AUGMENTED_H_
#define FIELDOPT_OPTIMIZATION_OBJECTIVE_AUGMENTED_H_

#include "objective.h"
#include "Settings/model.h"
#include "Simulation/results/results.h"

namespace Optimization {
namespace Objective {

using std::string;
using std::vector;
using namespace Simulation::Results;

class Augmented : public Objective {
 public:
/*!
 * \brief Augmented
 * \param
 */

  Augmented(Settings::Optimizer *settings,
            Simulation::Results::Results *results,
            Model::Model *model);

  struct Term {
    string prop_name_str;
    Results::Property prop_name;
    Results::Property prop_type;
    Results::Property prop_spec;
    double coeff;
  };

  void setUpAugTerms();

  double value() const override;

 private:
  Simulation::Results::Results *results_;
  Settings::Optimizer *settings_;

  vector<Augmented::Term*> terms_;

  string cl_ = "Augmented";

};

}
}

#endif //FIELDOPT_OPTIMIZATION_OBJECTIVE_AUGMENTED_H_
