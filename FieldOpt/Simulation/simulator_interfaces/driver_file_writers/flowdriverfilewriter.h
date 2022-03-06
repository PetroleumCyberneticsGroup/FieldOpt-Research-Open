/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2021 Mathias Bellout
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

#ifndef FIELDOPT_FLOWDRIVERFILEWRITER_H
#define FIELDOPT_FLOWDRIVERFILEWRITER_H

#include "Settings/settings.h"
#include "Settings/simulator.h"
#include "Model/model.h"

namespace Simulation {
class FlowSimulator;
}

namespace Simulation {
class FlowDriverFileWriter {
  friend class FlowSimulator;
 private:
  friend class ::Simulation::FlowSimulator;
  FlowDriverFileWriter(::Settings::Settings *settings, Model::Model *model);
  void WriteDriverFile(QString output_dir);

  Model::Model *model_;
  ::Settings::Settings *settings_;

  string md_ = "Simulation::simulator_interfaces::driver_file_writer";
  string cl_ = "FlowDriverFileWriter";
  ::Settings::VerbParams vp_;

  QString output_driver_file_name_; //!< Path to the driver file to be written.
  QString GetCompdatString();
};
}

#endif // FIELDOPT_FLOWDRIVERFILEWRITER_H
