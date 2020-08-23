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

#ifndef FIELDOPT_IX_SIMULATOR_H
#define FIELDOPT_IX_SIMULATOR_H

#include "Simulation/simulator_interfaces/simulator.h"

namespace Simulation {

class IXSimulator : public Simulator {

 public:
  IXSimulator(Settings::Settings *settings, Model::Model *model);
  void Evaluate() override;
  bool Evaluate(int timeout, int threads) override;
  bool Evaluate(const Settings::Ensemble::Realization &realization, int timeout, int threads) override;
  void WriteDriverFilesOnly() override;
  void CleanUp() override;
 protected:
  void UpdateFilePaths() override;
 private:
  void copyDriverFiles();
  void setResultPath();

  QString deck_name_;
  QString result_path_;

  string cl_ = "IXSimulator";
};

}

#endif //FIELDOPT_IX_SIMULATOR_H
