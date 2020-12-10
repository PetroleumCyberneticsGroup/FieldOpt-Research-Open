/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
<mathias.bellout@gmail.com>

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

#ifndef FIELDOPT_COMPARTMENT_H
#define FIELDOPT_COMPARTMENT_H

#include "wellbore/completions/icd.h"
#include "wellbore/completions/packer.h"

#include "Model/wells/wellbore/wellblock.h"

namespace Model {
namespace Wells {

using namespace Wellbore::Completions;

/*!
 * The Compartment struct models a compartment made up of
 * two packers and an ICD. The ICD is placed at the beginning
 * of the compartment (i.e. closest to the heel of the well).
 * This class is used by segmented wells to build segments.
 */
class Compartment {
 public:
  Compartment();
  Compartment(double start_md, double end_md,
              double start_tvd, double end_tvd,
              double well_length,
              std::vector<Wellbore::WellBlock *> block_range,
              std::vector<int> block_idx_range,
              const Settings::Model::Well &well_settings,
              Properties::VarPropContainer *variable_container,
              std::vector<Compartment> &compartments_);

  double GetLength(const double &well_length) const;
  std::vector<Wellbore::WellBlock *> GetBlockRange() const;
  std::vector<int> GetBlockIdxange() const;
  double GetTVDDifference() const;

  /*!
   * Update the ICD MD to fit new start packer MDs
   */
  void Update();

  Packer *start_packer;
  Packer *end_packer;
  ICD *icd;

  std::vector<Wellbore::WellBlock *> block_range_;
  std::vector<int> block_idx_range_;

  string im_ = "", wm_ = "", em_ = "";
  string md_ = "Model";
  string cl_ = "Well";
  Settings::VerbParams vp_;
};

}
}
#endif //FIELDOPT_COMPARTMENT_H
