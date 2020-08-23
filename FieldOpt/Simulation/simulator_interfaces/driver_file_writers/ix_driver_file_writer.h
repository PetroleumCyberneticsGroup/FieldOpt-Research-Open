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

#ifndef FIELDOPT_IX_DRIVER_FILE_WRITER_H
#define FIELDOPT_IX_DRIVER_FILE_WRITER_H

#include "Model/model.h"
#include "Settings/settings.h"

namespace Simulation {

/*!
 * \brief This class implements driver file writing for Schlumberger's
 * INTERSECT reservoir simulator. It writes to the [DECK_NAME]_fm_edits.ixf file.
 */
class IXDriverFileWriter {
 public:

  IXDriverFileWriter(Settings::Settings *settings,
                     Model::Model *model);
  void WriteDriverFile(std::string fm_edits_path);

 private:
  Model::Model *model_;
  ::Settings::Settings *settings_;

  string im_ = "", wm_ = "", em_ = "";
  ::Settings::VerbParams vp_;
  string md_ = "Simulation::sim_interfaces::driver_file_writers";
  string cl_ = "IXDriverFileWriter";

};

}

#endif //FIELDOPT_IX_DRIVER_FILE_WRITER_H
