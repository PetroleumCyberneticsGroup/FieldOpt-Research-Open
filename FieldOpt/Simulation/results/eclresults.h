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

#ifndef ECLRESULTS_H
#define ECLRESULTS_H

#include "results.h"
#include "ERTWrapper/eclsummaryreader.h"
#include <QHash>
#include <set>

namespace Simulation {
namespace Results {

/*!
 * \brief The ECLResults class uses the ECLSummaryReader class in the ERTWrapper library
 * to read properties from ECLIPSE summary files.
 */
class ECLResults : public Results
{
 public:
  explicit ECLResults(Settings::Settings *settings);

  void ReadResults(QString file_path);
  static void ReadResults(const QString& file_path,
                          QString build_dir);
  void DumpResults() override;

  double GetValue(Property prop) override;
  double GetValue(Property prop, int time_index) override;
  double GetValue(Property prop, QString well) override;
  double GetValue(Property prop, QString well, int time_index) override;

  vector<double> GetValueVector(Property prop) override;
  vector<double> GetValueVector(Property prop, string wn) override;
  vector<vector<double>> GetValVectorSeg(Property prop, string wn);

  VectorXd GetValueVectorXd(Property prop) override;
  VectorXd GetValueVectorXd(Property prop, string wn);
  vector<VectorXd> GetValVectorSegXd(Property prop, string wn);

 private:
  QString file_path_;
  ERTWrapper::ECLSummary::ECLSummaryReader *smry_reader_;

  Settings::Settings* settings_;

  string im_ = "", wm_ = "", em_ = "";
  string cl_ = "eclresults";
};

}
}

#endif // ECLRESULTS_H
