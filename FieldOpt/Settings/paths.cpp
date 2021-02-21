/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb.pcg@gmail.com>

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

#include "paths.h"
#include "Utilities/printer.hpp"

using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::stringstream;
using std::runtime_error;

using Printer::ext_warn;
using Printer::error;

Paths::Paths() {}

void Paths::SetPath(Paths::Path path, const string& path_string,
                    bool skip_check, Settings::VerbParams vp) {
  vp_=vp;

  if (path >= 0 && !FileExists(path_string, vp_, md_, cl_) && ! skip_check) {
    stringstream ss;
    ss << "Cannot set " << GetPathDesc(path);
    ss << " path to non-existing file (" << path_string << ")";
    error(ss.str());

    string em = Paths::GetPathDesc(path) + " not found at " + path_string;
    throw runtime_error(em);

  } else if (path < 0 && !DirExists(path_string, vp_, md_, cl_) && ! skip_check) {
    cerr << "Cannot set " << GetPathDesc(path);
    cerr << " path to non-existing directory (";
    cerr << path_string << ")" << endl;

    string em = Paths::GetPathDesc(path) + " not found at " + path_string;
    throw runtime_error(em);

  } else {
    paths_[path] = path_string;
  }
}

void Paths::CopyPath(Paths::Path path0, Paths::Path path1) {
  paths_[path0] = paths_[path1];
}


bool Paths::IsSet(Paths::Path path) {
  return paths_.count(path) > 0;
}

string Paths::GetPath(Paths::Path path) {
  if (!IsSet(path)) {
    string wm = "Getting unset variable (" + GetPathDesc(path) + ")";
    ext_warn(wm, md_, cl_, vp_.lnw);
    return "UNSET";

  } else {
    return paths_[path];
  }
}

QString Paths::GetPathQstr(Paths::Path path) {
  return QString::fromStdString(GetPath(path));
}

void Paths::ShowPaths() {
  for (auto & path : paths_) {
    string pth = Paths::GetPathDesc(path.first);
    Printer::pad_text(pth, 30);
    cout << FYELLOW <<  pth << ": " << AEND
         << path.second << endl;
  }
}

const string &Paths::GetPathDesc(Paths::Path path) const {
  return path_descriptions.at(path);
}
