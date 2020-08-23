/***********************************************************
Created: 12.11.2015 2015 by einar

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

#ifndef ECLDRIVERFILEWRITER_H
#define ECLDRIVERFILEWRITER_H

#include "Settings/settings.h"
#include "Settings/simulator.h"
#include "Model/model.h"
#include "driver_parts/ecl_driver_parts/schedule_insets.h"

namespace Simulation {
class ECLSimulator;
}

namespace Simulation {

/*!
 * \brief The EclDriverFileWriter class writes driver files
 * that can be executed by the ECL100 reservoir simulator.
 * This class should _only_ be used by the ECLSimulator class.
 */
class EclDriverFileWriter
{
 private:
  friend class ::Simulation::ECLSimulator;
  EclDriverFileWriter(::Settings::Settings *settings,
                      Model::Model *model);
  void WriteDriverFile(QString schedule_file_path);
  std::string buildActionStrings();

  Model::Model *model_;
  ::Settings::Settings *settings_;
  ECLDriverParts::ScheduleInsets insets_;
  bool use_actionx_;

  string im_ = "", wm_ = "", em_ = "";
  ::Settings::VerbParams vp_;
  string md_ = "Simulation::sim_interfaces::driver_file_writers";
  string cl_ = "EclDriverFileWriter";
};

}


#endif // ECLDRIVERFILEWRITER_H
