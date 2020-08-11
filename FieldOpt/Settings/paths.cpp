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

#include "paths.h"
#include "Utilities/printer.hpp"

using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::stringstream;
using std::runtime_error;

Paths::Paths() {}

// =========================================================
void Paths::SetPath(Paths::Path path, string path_string,
    bool skip_check) {
  
  if (path >= 0 && !FileExists(path_string) && ! skip_check) {
    
    stringstream ss;
    ss << "Cannot set " << GetPathDesc(path)
       << " path to non-existing file (" << path_string << ")";
    Printer::error(ss.str());
    throw runtime_error(
        Paths::GetPathDesc(path)
        + " not found at " + path_string);
    
  } else if (path < 0 && !DirExists(path_string) && ! skip_check) {

    cerr << "Cannot set " << GetPathDesc(path)
    << " path to non-existing directory (" << path_string
    << ")" << endl;
    throw runtime_error(
        Paths::GetPathDesc(path)
        + " not found at " + path_string);
    
  } else {
    paths_[path] = path_string;
  }
}

bool Paths::IsSet(Paths::Path path) {
  return paths_.count(path) > 0;
}

string Paths::GetPath(Paths::Path path) {
  if (!IsSet(path)) {
   Printer::ext_warn("Getting unset variable ("
       + GetPathDesc(path) + ")", "Settings",
                     "Paths");
    return "UNSET";
    
  } else {
    return paths_[path];
  }
}

QString Paths::GetPathQstr(Paths::Path path) {
 return QString::fromStdString(GetPath(path));
}

void Paths::ShowPaths() {
  for (auto it = paths_.begin(); it != paths_.end(); it++) {
    string pth = Paths::GetPathDesc(it->first);
    Printer::pad_text(pth, 30);
    cout << FYELLOW <<  pth << ": " << AEND
    << it->second << endl;
  }
}

const string &Paths::GetPathDesc(Paths::Path path) const {
  return path_descriptions.at(path);
}
