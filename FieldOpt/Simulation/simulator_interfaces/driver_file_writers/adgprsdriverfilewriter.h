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

#ifndef ADGPRSDRIVERFILEWRITER_H
#define ADGPRSDRIVERFILEWRITER_H

#include "Settings/settings.h"
#include "Settings/simulator.h"
#include "Model/model.h"

namespace Simulation {
class AdgprsSimulator;
}

namespace Simulation {

/*!
 * \brief The AdgprsDriverFileWriter class is responsible for writing AD-GPRS driver files representing
 * the model. It uses many of the parts created for the ECL simulator, as many of the keywords are
 * more or less identical.
 */
class AdgprsDriverFileWriter
{
 private:
  friend class ::Simulation::AdgprsSimulator;
  AdgprsDriverFileWriter(::Settings::Settings *settings, Model::Model *model);
  void WriteDriverFile(QString output_dir);

  Model::Model *model_;
  Settings::Settings *settings_;
  ::Settings::VerbParams vp_;
};

}

#endif // ADGPRSDRIVERFILEWRITER_H
