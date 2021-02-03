/***********************************************************
Copyright (C) 2020-2021 Mathias Bellout
<chakibbb.pcg@gmail.com>

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

class Global {
  friend class Settings;
  // Since Settings is friend of Global,
  // Settings can access private members of Global

 public:
  Global(QJsonObject json_global);

  //!< Name to be used for the run. Output file
  //!< and folder names are derived from this.
  string name() const { return name_; }

  //!< Get the value for the bookkeeper tolerance.
  //!< Used by the Bookkeeper in the Runner library.
  double bookkeeper_tol() const { return bookkeeper_tol_; }

  VerbParams verbParams() { return verb_params_; };
  void showVerbParams(string sin="");
  void showVerbParams(VerbParams &vp, string sin="");

 protected:
  string name_;
  double bookkeeper_tol_;

 private:
  string im_ = "", wm_ = "", em_ = "";
  string md_ = "Settings";
  string cl_ = "Global";
  VerbParams verb_params_;
  VerbParams parseVerbParams(QJsonObject json_verb);
};

}

#endif //FIELDOPT_SETTINGS_GLOBAL_H_


