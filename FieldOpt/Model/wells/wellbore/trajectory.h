/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

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

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "wellblock.h"
#include "wellspline.h"
#include "completions/completion.h"
#include "completions/perforation.h"
#include "Reservoir/grid/eclgrid.h"
#include "Settings/model.h"
#include "Model/properties/var_prop_container.h"
#include "Model/properties/property.h"
#include "pseudo_cont_vert.h"

#include <QList>
#include <WellIndexCalculation/wicalc_rixx.h>

namespace Model {
namespace Wells {
namespace Wellbore {

using Printer::ext_info;
using Printer::ext_warn;
using Printer::info;

using Printer::num2str;
using std::runtime_error;

using WType=Settings::Model::WellDefinitionType;

class WellSpline;

/*!
 * \brief The Trajectory class describes the trajectory of
 * the wellbore as a set of well blocks in the reservoir grid.
 *
 * A trajectory may be described either as a set of well
 * blocks directly, or as a spline which a set of well
 * blocks is calculated from. Either way it is always
 * perceived as a set of blocks from the outside.
 *
 * \todo Initialize ICDs.
 */
class Trajectory
{
 public:
  Trajectory(::Settings::Model::Well well_settings,
             Properties::VarPropContainer *variable_container,
             Reservoir::Grid::Grid *grid,
             Reservoir::WellIndexCalculation::wicalc_rixx *wic);

  //!< Get the well block at index (i,j,k).
  WellBlock *GetWellBlock(int i, int j, int k);

  //!< Get a list containing all well blocks.
  QList<WellBlock *> *GetWellBlocks();

  //!< Update the well blocks, in particular the ones defined by a spline.
  void UpdateWellBlocks();
  int GetTimeSpentInWic() const;
  Settings::Model::WellDefinitionType GetDefinitionType();

  //!< Get length of the wellbore (measured depth from the heel to the toe)
  double GetLength() const;

  //!< Get the wellblock surrounding the given MD.
  WellBlock * GetWellBlockByMd(double md, string dmsg="");
  std::vector<WellBlock *> GetWellBlocksByMdRange(double start_md,
                                                  double end_md,
                                                  string dmsg="") const;

  WellBlock* dummy_block = new WellBlock(0,0,0);

  //!< Get the measured depth for the entry point to the block.
  double GetEntryMd(const WellBlock *wb) const;

  //!< Get the measured depth for the exit point from the block.
  double GetExitMd(const WellBlock *wb) const;

  //!< Get length of trajectory (summed distance b/e defining points).
  double GetSplineLength() const;
  WellSpline *GetWellSpline() const { return well_spline_; }


 private:
  Settings::Model::WellDefinitionType definition_type_;
  QList<WellBlock *> *well_blocks_;

  //!< Used to defined trajectories with a spline.
  //!< When used, this generates the well blocks.
  WellSpline *well_spline_;

  //!< A pseudo-continuous vertical well.
  PseudoContVert *pseudo_cont_vert_;

  void initializeWellBlocks(Settings::Model::Well well,
                            Properties::VarPropContainer *variable_container);

  // Calculate direction of penetration for all well blocks
  void calculateDirectionOfPenetration();

  /*!
   * Convert the list of well blocks in the well settings
   * object to a well spline, with the number of spline
   * points specified in the settings object (defaulting to 2).
   *
   * The spline points will be selected naively. E.g. for
   * three spline points: the heel point will be the center
   * of the first well block in the list; the center point
   * will be the center of the middle well block in the list
   * ( list[size/2] ) and the toe point the center of the
   * last well block;
   */
  void convertWellBlocksToWellSpline(Settings::Model::Well &well_settings,
                                     Reservoir::Grid::Grid *grid);

  //!< Indicates if the well should only be able to vary
  //!< in the x-y plane (z variables will not be created).
  bool is_2d_;

  void printWellBlocks();

  void E(string m) const {
    m = "[mod: " + md_ + "] [cls: " + cl_ + "] " + m;
    throw runtime_error(m);
  };

  string im_ = "", wm_ = "", em_ = "";
  string md_ = "Model/wells/wellbore";
  string cl_ = "Trajectory";
  Settings::VerbParams vp_;

};

}
}
}

#endif // TRAJECTORY_H
