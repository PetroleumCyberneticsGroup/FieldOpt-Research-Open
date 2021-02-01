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

#include "ix_driver_file_writer.h"
#include "driver_parts/ix_driver_parts/report_tuning.hpp"
#include "driver_parts/ix_driver_parts/ix_control.hpp"
#include "driver_parts/ix_driver_parts/flow_control_device.hpp"

namespace Simulation {

using Utilities::FileHandling::WriteStringToFile;
using Printer::ext_info;

IXDriverFileWriter::IXDriverFileWriter(Settings::Settings *settings,
                                       Model::Model *model) {
  model_ = model;
  settings_ = settings;
  vp_ = settings_->global()->verbParams();
  // settings_->global()->showVerbParams("[ @IXDriverFileWriter ]");
}

void IXDriverFileWriter::WriteDriverFile(std::string fm_edits_path) {

  std::string fm_edits = "MODEL_DEFINITION\n\n";
  fm_edits += IXParts::FieldManagementStandardReport();
  fm_edits += IXParts::EclReports();

  bool use_segments = false;
  auto icd_segments = std::vector<int>();
  for (auto well : *model_->wells()) {
    if (well->IsSegmented()) {
      use_segments = true;
      auto icd_seg_idxs = well->GetICDSegmentIndices();
      for (auto seg : icd_seg_idxs) {
        // Check if segment is already in list
        if(std::find(icd_segments.begin(), icd_segments.end(), seg) == icd_segments.end()) {
          icd_segments.push_back(seg); // Add if segment is not already there
        }
      }
    }
  }

  fm_edits += IXParts::XYPlotSummaryReport(use_segments, icd_segments);
  for (auto well : *model_->wells()) {
    fm_edits += IXParts::CreateControlEntries(well);
    if (well->HasSimpleICVs()) {
      fm_edits += IXParts::AllWellDevices(well);
    }
  }

  im_ = "fm_edits: " + fm_edits;
  ext_info(im_, md_, cl_, vp_.lnw);

  im_ = "  fm_edits_path: " + fm_edits_path;
  ext_info(im_, md_, cl_, vp_.lnw);
  // if(vp_.vSIM >= 3) {
  // }

  WriteStringToFile(QString::fromStdString(fm_edits),
                    QString::fromStdString(fm_edits_path));

}
}
