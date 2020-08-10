/***********************************************************
Copyright (C) 2020-2021 Mathias Bellout
<chakibbb-pcg@gmail.com>

Created by bellout on 8/10/20.

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

#include "global.h"
#include "settings_exceptions.h"

using namespace Utilities::FileHandling;

namespace Settings {

Global::Global(QJsonObject json_global) {

  name_ = json_global["Name"].toString().toStdString();
  bookkeeper_tol_ = json_global["BookkeeperTol"].toDouble();
  if (bookkeeper_tol_ < 0.0) {
    throw UnableToParseGlobalSectionException(
      "Bookkeeper tolerance must larger than zero.");
  }

  QJsonObject json_verb = json_global["VerbConf"].toObject();
  parseVerbParams(json_verb);
  // showVerbParams();
}

VerbParams Global::parseVerbParams(QJsonObject json_verb) {
  verb_params_.lnw=json_verb["iBoxLineWidth"].toInt();

  verb_params_.vMOD=json_verb["Model"].toInt();
  verb_params_.vOPT=json_verb["Optim"].toInt();
  verb_params_.vWIC=json_verb["WICal"].toInt();
  verb_params_.vSIM=json_verb["Simul"].toInt();
  verb_params_.vRUN=json_verb["Runnr"].toInt();
  verb_params_.vRES=json_verb["Resvr"].toInt();
  verb_params_.vSET=json_verb["Setng"].toInt();
  verb_params_.vUTI=json_verb["Utils"].toInt();
}

void Global::showVerbParams() {

  cout << "lnw:" << verb_params_.lnw << endl;

  cout << "vMOD:" << verb_params_.vMOD << endl;
  cout << "vOPT:" << verb_params_.vOPT << endl;
  cout << "vWIC:" << verb_params_.vWIC << endl;
  cout << "vSIM:" << verb_params_.vSIM << endl;
  cout << "vRUN:" << verb_params_.vRUN << endl;
  cout << "vRES:" << verb_params_.vRES << endl;
  cout << "vSET:" << verb_params_.vSET << endl;
  cout << "vUTI:" << verb_params_.vUTI << endl;

}

}