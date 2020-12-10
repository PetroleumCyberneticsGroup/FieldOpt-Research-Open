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

#include "adgprs_wellcontrols.h"

namespace Simulation {
namespace AdgprsDriverParts {

WellControls::WellControls(QList<Model::Wells::Well *> *wells, QList<double> control_times)
    : Simulation::ECLDriverParts::WellControls(wells, control_times)
{}

QString WellControls::GetPartString()
{
    return Simulation::ECLDriverParts::WellControls::GetPartString();
}

QString WellControls::createTimeEntry(int time, int prev_time)
{
    if (time == 0) return "";

    // Find prev time step
    int prev_step = 0;
    for (int i = 0; i < time_entries_.keys().size(); ++i) {
        if (time_entries_.keys()[i] == time) {
            prev_step = time_entries_.keys()[i-1];
            break;
        }
    }
    int dt = time-prev_step;
    return QString("TSTEP\n\t%1/\n").arg(dt);
}


}}
