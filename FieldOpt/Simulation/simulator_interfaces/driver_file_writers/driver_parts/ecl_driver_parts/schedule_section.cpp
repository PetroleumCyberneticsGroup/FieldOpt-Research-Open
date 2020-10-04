/***********************************************************
Created: 18.11.2015 2015 by einar

Copyright (C) 2015
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

#include "schedule_section.h"
#include <ctime>
#include "Utilities/execution.hpp"

namespace Simulation {
namespace ECLDriverParts {

using Printer::num2str;

Schedule::Schedule(QList<Model::Wells::Well *> *wells,
                   const QList<int>& control_times,
                   const QList<int>& start_date,
                   ScheduleInsets &insets,
                   Settings::Settings *settings) : ECLDriverPart(settings) {

  schedule_time_entries_ = QList<ScheduleTimeEntry>();
  for (int ts : control_times) {
    time_DMY ts_DMY;
    GetControlDate(ts_DMY, start_date, ts);

    ScheduleTimeEntry time_entry =
      ScheduleTimeEntry(ts,
                        ts_DMY,
                        Welspecs(wells, ts),
                        Compdat(wells, ts),
                        WellControls(wells, control_times, ts),
                        Welsegs(wells, ts),
                        Compsegs(wells, ts),
                        Wsegvalv(wells, ts));
    schedule_time_entries_.append(time_entry);
  }

  if (insets.HasInset(-1)) {
    schedule_.append(QString::fromStdString(insets.GetInset(-1)));
  }

  stringstream ss;
  for (auto time_entry : schedule_time_entries_) {

    schedule_.append(time_entry.welspecs.GetPartString());

    if (insets.HasInset(time_entry.control_time)) {
      schedule_.append(QString::fromStdString(insets.GetInset(time_entry.control_time)));
    }

    schedule_.append(time_entry.compdat.GetPartString());
    schedule_.append(time_entry.welsegs.GetPartString());
    schedule_.append(time_entry.compsegs.GetPartString());
    schedule_.append(time_entry.wsegvalv.GetPartString());
    schedule_.append(time_entry.well_controls.GetPartString());

    if (vp_.vRUN >= 3) {
      ss << "TIME: " << time_entry.control_time << "|";
      ss << "WELLSPEC |" << time_entry.welspecs.GetPartString().toStdString();
      ss << "COMPDAT |" << time_entry.compdat.GetPartString().toStdString();
      ss << "WELSEGS |" << time_entry.welsegs.GetPartString().toStdString();
      ss << "COMPSEGS |" << time_entry.compsegs.GetPartString().toStdString();
      ss << "WSEGVALV |" << time_entry.wsegvalv.GetPartString().toStdString();
      ss << "WCNTRLS |" << time_entry.well_controls.GetPartString().toStdString();
    }

  }
  schedule_.append("\n\n");
}

QString Schedule::GetPartString() const {
  return schedule_;
}

Schedule::ScheduleTimeEntry::ScheduleTimeEntry(int control_time,
                                               time_DMY control_time_DMY,
                                               Welspecs welspecs,
                                               Compdat compdat,
                                               WellControls well_controls,
                                               Welsegs welsegs,
                                               Compsegs compsegs,
                                               Wsegvalv wsegvalv) {
  this->control_time = control_time;
  this->control_time_DMY = control_time_DMY;
  this->welspecs = welspecs;
  this->compdat = compdat;
  this->well_controls = well_controls;
  this->welsegs = welsegs;
  this->compsegs = compsegs;
  this->wsegvalv = wsegvalv;
}

void Schedule::GetControlDate(time_DMY &ts_DMY, QList<int> sd, int control_time) const {

  const time_t sec_per_day = 24*60*60;

  tm start_date;
  time_t start_date_secs;
  time_t new_date_secs;

  start_date.tm_sec = 0; // seconds of minutes from 0 to 61
  start_date.tm_hour = 0; // hours of day from 0 to 24
  start_date.tm_min = 0; // minutes of hour from 0 to 59

  start_date.tm_mday = sd.at(0); // day of month from 1 to 31
  start_date.tm_mon = sd.at(1) - 1; // month of year from 0 to 11
  start_date.tm_year = sd.at(2) - 1900;	// year, begining 1900

  start_date.tm_wday = 0;	// days since sunday
  start_date.tm_yday = 0; // days since January 1st
  start_date.tm_isdst = 0; // hours of daylight savings time

  start_date_secs = mktime(&start_date);
  new_date_secs = start_date_secs + control_time * sec_per_day;

  ts_DMY.day = localtime(&new_date_secs)->tm_mday;
  ts_DMY.month = localtime(&new_date_secs)->tm_mon + 1;
  ts_DMY.year = localtime(&new_date_secs)->tm_year + 1900;

  ts_DMY.day_str = to_string(ts_DMY.day);
  ts_DMY.year_str = to_string(ts_DMY.year);
  if (ts_DMY.month == 1) { ts_DMY.month_str = "JAN"; }
  else if (ts_DMY.month ==  2) { ts_DMY.month_str = "FEB"; }
  else if (ts_DMY.month ==  3) { ts_DMY.month_str = "MAR"; }
  else if (ts_DMY.month ==  4) { ts_DMY.month_str = "APR"; }
  else if (ts_DMY.month ==  5) { ts_DMY.month_str = "MAY"; }
  else if (ts_DMY.month ==  6) { ts_DMY.month_str = "JUN"; }
  else if (ts_DMY.month ==  7) { ts_DMY.month_str = "JUL"; }
  else if (ts_DMY.month ==  8) { ts_DMY.month_str = "AUG"; }
  else if (ts_DMY.month ==  9) { ts_DMY.month_str = "SEP"; }
  else if (ts_DMY.month == 10) { ts_DMY.month_str = "OCT"; }
  else if (ts_DMY.month == 11) { ts_DMY.month_str = "NOV"; }
  else if (ts_DMY.month == 12) { ts_DMY.month_str = "DEC"; }

  if (vp_.vRUN >= 3) {
    stringstream im;
    im << "TIME: " << num2str(control_time);
    im << " -- D: " << ts_DMY.day_str << " " << ts_DMY.month_str << " " << ts_DMY.year_str;
    string linux_date = "date '+%d %b %Y' -d \"1 Oct 2019+" + num2str(control_time) + " days\" | tr -d '\\n'";
    im << " -- " << Utilities::Unix::ExecStr(linux_date.c_str());

    // im << " -- asctime date: " << asctime(localtime(&new_date_secs));
    info(im.str(), vp_.lnw);
  }
}

}
}
