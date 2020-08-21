/***********************************************************
Copyright (C) 2016
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020- Mathias Bellout
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

#ifndef ADGPRSRESULTS_H
#define ADGPRSRESULTS_H

#include "results.h"
#include <QHash>
#include "Hdf5SummaryReader/hdf5_summary_reader.h"

namespace Simulation {
namespace Results {

/*!
 * \brief AdgprsResults provides access to the results
 * of an ADGPRS simulation.
 *
 * TODO: Support well properties. This is problematic
 * because in the summary the wells only have a number
 * that I can't reliably map to a name.
 */
class AdgprsResults : public Results
{
 public:
  explicit AdgprsResults(Settings::Settings *settings);

  double GetValue(int well_nr, Property prop);
  double GetValue(int well_nr, Property prop, int time_index);

  // Results interface
 public:
  void ReadResults(QString file_path);
  void DumpResults() override;
  double GetValue(Property prop) override;
  double GetValue(Property prop, QString well) override;
  double GetValue(Property prop, int time_index) override;
  double GetValue(Property prop, QString well, int time_index) override;

  std::vector<double> GetValueVector(Property prop) override;
  VectorXd GetValueVectorXd(Property prop) override;

 private:
  QString file_path_;
  Hdf5SummaryReader *smry_reader_;
};

}}

#endif // ADGPRSRESULTS_H
