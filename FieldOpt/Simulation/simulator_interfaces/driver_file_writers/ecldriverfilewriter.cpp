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

#include "ecldriverfilewriter.h"
#include "driver_parts/ecl_driver_parts/schedule_section.h"
#include "driver_parts/ecl_driver_parts/actionx.hpp"
#include "Simulation/simulator_interfaces/simulator_exceptions.h"

#include "Utilities/filehandling.hpp"
#include <Utilities/printer.hpp>
#include "Utilities/verbosity.h"

namespace Simulation {

using namespace ECLDriverParts;
using namespace Utilities::FileHandling;

EclDriverFileWriter::EclDriverFileWriter(Settings::Settings *settings, Model::Model *model)
{
  model_ = model;
  settings_ = settings;
  use_actionx_ = settings->simulator()->use_actionx();

  if (settings->paths().IsSet(Paths::SIM_SCH_INSET_FILE)) {
    insets_ = ECLDriverParts::ScheduleInsets(settings->paths().GetPath(Paths::SIM_SCH_INSET_FILE));
  }
}

void EclDriverFileWriter::WriteDriverFile(QString schedule_file_path)
{
  if (VERB_SIM >= 2) {
    auto fp = schedule_file_path.toStdString();
    Printer::ext_info("Writing driver file to " + fp + ".", "Simulation", "EclDriverFileWriter");
  }
  assert(FileExists(schedule_file_path));

  if (use_actionx_ == false) {
    Schedule schedule = ECLDriverParts::Schedule(model_->wells(), settings_->model()->control_times(), insets_);
    model_->SetCompdatString(schedule.GetPartString());
    Utilities::FileHandling::WriteStringToFile(schedule.GetPartString(), schedule_file_path);
  }
  else {
    Utilities::FileHandling::WriteStringToFile(QString::fromStdString(buildActionStrings()), schedule_file_path);
  }
  if (VERB_SIM >= 2) {
    Printer::ext_info("Wrote driver string to" + schedule_file_path.toStdString(), "Simulation", "EclDriverFileWriter");
  }
}
std::string EclDriverFileWriter::buildActionStrings() {
    if (VERB_SIM >= 2) {
        Printer::ext_info("Generating action strings", "Simulation", "EclDriverFileWriter");
    }
    std::string actions = "";

    std::string icv_actions = "";
    for (auto well : *model_->wells()) {
        if (well->HasSimpleICVs()) {
            auto wsegv = ECLDriverParts::Wsegvalv(well);
            icv_actions += wsegv.GetPartString().toStdString() + "\n";
        }
    }

    std::string ctrl_actions = "";
    if (model_->wells()->first()->controls()->size() > 0) {
        Schedule schedule = ECLDriverParts::Schedule(model_->wells(), settings_->model()->control_times(), insets_);
        auto schedule_time_entries = schedule.GetScheduleTimeEntries();
        for (auto entry : schedule_time_entries) {
            int ctrl_time = 0;
            entry.control_time == 0 ? ctrl_time = 1 : ctrl_time = entry.control_time;
            for (QString entry_part : entry.well_controls.GetWellEntryList()) {
                ctrl_actions += ActionX::ACTIONX("CTR_" + Printer::num2str(entry.control_time),
                        ECLDriverParts::ActionX::ACTX_LHQuantity::Day,
                        ECLDriverParts::ActionX::ACTX_Operator::GE, 
                        ctrl_time,
                        entry_part.toStdString()
                );
                ctrl_actions += "\n";
            }
        }
    }
    else if (VERB_SIM >= 2) {
        Printer::ext_warn("First well did not have controls; assuming none do.", "Simulation", "ECLDriverFileWriter");
    }

    actions += ActionX::ACTIONX("ICVS_T0", ECLDriverParts::ActionX::ACTX_LHQuantity::Day,
        ECLDriverParts::ActionX::ACTX_Operator::GE, 0, icv_actions);
    actions += ctrl_actions;
    return actions;
}

}
