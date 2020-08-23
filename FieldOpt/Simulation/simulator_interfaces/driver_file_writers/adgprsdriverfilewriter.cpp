/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

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

#include "adgprsdriverfilewriter.h"
#include "Simulation/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/welspecs.h"
#include "Simulation/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/compdat.h"
#include "Simulation/simulator_interfaces/driver_file_writers/driver_parts/adgprs_driver_parts/wellstre.h"
#include "Simulation/simulator_interfaces/driver_file_writers/driver_parts/adgprs_driver_parts/adgprs_wellcontrols.h"
#include <iostream>
#include "Utilities/filehandling.hpp"

namespace Simulation {

using namespace Utilities::FileHandling;
using namespace Utilities::FileHandling;
using std::runtime_error;

AdgprsDriverFileWriter::AdgprsDriverFileWriter(Settings::Settings *settings,
                                               Model::Model *model) {
  model_ = model;
  settings_ = settings;
  vp_ = settings_->global()->verbParams();
}

void AdgprsDriverFileWriter::WriteDriverFile(QString output_dir) {
  auto welspecs = ECLDriverParts::Welspecs(model_->wells());
  auto compdat = ECLDriverParts::Compdat(model_->wells());
  model_->SetCompdatString(compdat.GetPartString());

  auto wellstre = AdgprsDriverParts::Wellstre(model_->wells(), settings_->simulator()->fluid_model());
  auto wellcontrols = AdgprsDriverParts::WellControls(model_->wells(), settings_->model()->control_times());

  if (!FileExists(output_dir+"/include/wells.in", vp_))
    throw runtime_error("Unable to find include/wells.in file to write to.");
  else WriteStringToFile(welspecs.GetPartString(), output_dir+"/include/welspecs.in");

  if (!FileExists(output_dir+"/include/compdat.in",vp_))
    throw runtime_error("Unable to find include/compdat.in file to write to.");
  else WriteStringToFile(compdat.GetPartString(), output_dir+"/include/compdat.in");

  if (!FileExists(output_dir+"/include/controls.in", vp_)) {
    throw runtime_error("Unable to find include/controls.in file to write to.");
  } else {
    QString str = wellstre.GetPartString() + wellcontrols.GetPartString();
    WriteStringToFile(str, output_dir+"/include/controls.in");
  }

}

}
