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

#include <simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/welspecs.h>
#include <simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/compdat.h>
#include <simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/wellcontrols.h>
#include <Utilities/filehandling.hpp>
#include "flowdriverfilewriter.h"

namespace Simulation {
FlowDriverFileWriter::FlowDriverFileWriter(::Settings::Settings *settings,
                                           Model::Model *model) {
  model_ = model;
  settings_ = settings;
  vp_ = settings_->global()->verbParams();
}

void FlowDriverFileWriter::WriteDriverFile(QString output_dir) {
  auto welspecs = ECLDriverParts::Welspecs(model_->wells());
  auto compdat = ECLDriverParts::Compdat(model_->wells());
  auto wellcontrols = ECLDriverParts::WellControls(model_->wells(), settings_->model()->control_times());
  model_->SetCompdatString(compdat.GetPartString());

  if (!Utilities::FileHandling::FileExists(output_dir+"/include/wells.in", vp_, md_, cl_)
    || !Utilities::FileHandling::FileExists(output_dir+"/include/welspecs.in", vp_, md_, cl_))
    throw std::runtime_error("Unable to find include/wells.in or include/welspecs.in file to write to.");
  else Utilities::FileHandling::WriteStringToFile(welspecs.GetPartString(), output_dir+"/include/welspecs.in");

  if (!Utilities::FileHandling::FileExists(output_dir+"/include/compdat.in", vp_, md_, cl_))
    throw std::runtime_error("Unable to find include/compdat.in file to write to.");
  else Utilities::FileHandling::WriteStringToFile(compdat.GetPartString(), output_dir+"/include/compdat.in");

  if (!Utilities::FileHandling::FileExists(output_dir+"/include/controls.in", vp_, md_, cl_))
    throw std::runtime_error("Unable to find include/controls.in file to write to.");
  else Utilities::FileHandling::WriteStringToFile(wellcontrols.GetPartString(), output_dir+"/include/controls.in");
}

QString FlowDriverFileWriter::GetCompdatString() {
  return ECLDriverParts::Compdat(model_->wells()).GetPartString();
}

}
