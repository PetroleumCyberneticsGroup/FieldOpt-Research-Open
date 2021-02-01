/***********************************************************
Copyright (C) 2015-2016
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#ifndef ECLGRIDREADER_EXCEPTIONS_H
#define ECLGRIDREADER_EXCEPTIONS_H

#include <stdexcept>
#include <string>

using std::string;
using std::runtime_error;

namespace ERTWrapper {

class GridNotReadException : public runtime_error {
 public:
  explicit GridNotReadException(const string& message)
    : runtime_error(message) {}
};

class InvalidIndexException : public runtime_error {
 public:
  explicit InvalidIndexException(const string& message)
    : runtime_error(message) {}
};

class SmryFileNotFoundAtPathExc : public runtime_error {
 public:
  explicit SmryFileNotFoundAtPathExc(const string &path)
    : runtime_error("No valid simulation case found at path  " + path) {}
};

class SmryVarDoesNotExistExc : public runtime_error {
 public:
  explicit SmryVarDoesNotExistExc(const string& message)
    : runtime_error(message) {}
};

class SmryTimeStepDoesNotExistExc : public runtime_error {
 public:
  explicit SmryTimeStepDoesNotExistExc(const string& message)
    : runtime_error(message) {}
};

}

#endif // ECLGRIDREADER_EXCEPTIONS_H
