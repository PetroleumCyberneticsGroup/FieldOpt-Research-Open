/***********************************************************
Copyright (C) 2015-2018
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

#ifndef COMPDAT_H
#define COMPDAT_H

#include "ecldriverpart.h"
#include "Model/wells/well.h"

namespace Simulation {
namespace ECLDriverParts {

/*!
 * \brief The Compdat class Generates the string for
 * the WELSPECS section in the ECLIPSE simulator.
 *
 * The information is taken from both the Model::Wells::Well
 * and Model::Wells::Wellbore::WellBlock classes.
 */
class Compdat : public ECLDriverPart
{
 public:
  /*!
   * Generate the compdat table for all wells.
   *
   * This will result in a non-viable deck unless all wells
   * are opened (i.e. have defined controls) at the first
   * time step.
   */
  explicit Compdat(QList<Model::Wells::Well *> *wells);

  Compdat()= default;

  /*!
   * Generate the compdat table for all wells opened at a
   * certain time step.
   * @param timestep Time step at which to generate table
   * for newly opened wells.
   */
  Compdat(QList<Model::Wells::Well *> *wells, double timestep);

  QString GetPartString() const override;

 private:
  QList<QStringList> createWellEntries(Model::Wells::Well *well);
  QStringList createBlockEntry(QString well_name,
                               double wellbore_radius,
                               Model::Wells::Wellbore::WellBlock *well_block);

  // void E(string m) const {
  //   m = "[mod: " + md_ + "] [cls: " + cl_ + "] " + m;
  //   throw runtime_error(m);
  // };

  string im_ = "", wm_ = "", em_ = "";
  string md_ = "Simulation/sim_interfaces/driver_file_writers/ecl_driver_parts";
  string cl_ = "Compdat";
  Settings::VerbParams vp_;

};

}
}
#endif // COMPDAT_H
