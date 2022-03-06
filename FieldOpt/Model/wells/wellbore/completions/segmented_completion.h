/***********************************************************
Copyright (C) 2015-2018
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

#ifndef FIELDOPT_SEGMENTED_COMPLETION_H
#define FIELDOPT_SEGMENTED_COMPLETION_H

#include "Model/properties/var_prop_container.h"
#include "completion.h"

namespace Model {
namespace Wells {
namespace Wellbore {
namespace Completions {

/*!
 * Parent class for completions typically used with the segmented well model.
 */
class SegmentedCompletion : public Completion {

 protected:
  SegmentedCompletion(Settings::Model::Well::Completion completion_settings,
                      Properties::VarPropContainer *variable_container);


 public:
  // Getters
  double placement() const;
  double md(const double &well_length) const;
  double tvd() const;
  double roughness() const;
  double diam() const;

  // Setters
  void setPlacement(const double &placement);
  void setTvd(const double &tvd);

 protected:
  /*!
   * Placement of a packer of ICD as a fraction of well length.
   */
  Properties::ContinuousProperty *placement_;

  /*!
   * True vertical depth. For ICDs and packers, this is the vertical depth
   * at the discrete location of the component along the well trajectory.
   * For tubing and annulus, this is the depth at the end-point for the
   * segment.
   */
  double tvd_; //!< True vertical depth

  /*!
   * Roughness of the component. For ICDs, this property will be set for
   * the generated parent segment, which will, by ECLs default, also be
   * used by the component.
   */
  double roughness_;

  /*!
   * Diameter of the component. For ICDs, this property will be set for
   * the generated parent segment, which will, by ECLs default, also be
   * used by the component.
   */
  double diameter_;
};

}
}
}
}

#endif //FIELDOPT_SEGMENTED_COMPLETION_H
