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

#ifndef FIELDOPT_SETTINGS_GLOBAL_H_
#define FIELDOPT_SETTINGS_GLOBAL_H_

#include "settings.h"
#include "Settings/paths.h"

#include <QStringList>

namespace Settings {

struct VerbParams {
  int lnw=140;
  int vMOD=0;
  int vOPT=0;
  int vWIC=0;
  int vSIM=0;
  int vRUN=0;
  int vRES=0;
  int vSET=0;
  int vUTI=0;
};

class Global {
  friend class Settings;

 public:
  Global(QJsonObject json_global);
  VerbParams verbParams() { return verb_params_; };
  void showVerbParams();

 private:
  string name_;
  double bookkeeper_tol_;

  VerbParams verb_params_;
  VerbParams parseVerbParams(QJsonObject json_verb);
};

}

#endif //FIELDOPT_SETTINGS_GLOBAL_H_


