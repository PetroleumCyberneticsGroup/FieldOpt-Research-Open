/***********************************************************
Created: 12.11.2015 2015 by einar

Copyright (C) 2015-2017
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

#ifndef SIMULATOR_EXCEPTIONS
#define SIMULATOR_EXCEPTIONS

#include <stdexcept>
#include <string>
#include <QString>

using std::string;

namespace Simulation {

class OutputDirectoryDoesNotExistException : public std::runtime_error {
 public:
  OutputDirectoryDoesNotExistException(const QString path)
    : std::runtime_error("Specified output dir does not exist: " + path.toStdString()) {}
};

class DriverFileDoesNotExistException : public std::runtime_error {
 public:
  DriverFileDoesNotExistException(const QString path)
    : std::runtime_error("Specified driver file does not exist: " + path.toStdString()) {}
};

class DriverFileInvalidException : public std::runtime_error {
 public:
  DriverFileInvalidException(const QString message)
    : std::runtime_error(message.toStdString()) {}
};

class UnableToWriteDriverFileException : public std::runtime_error {
 public:
  UnableToWriteDriverFileException(const QString message)
    : std::runtime_error(message.toStdString()) {}
};

class UnableToFindKeywordException : public std::runtime_error {
 public:
  UnableToFindKeywordException(const QString keyword)
    : std::runtime_error("Unable to find keyword in driver file: " + keyword.toStdString()) {}
};

}

#endif // SIMULATOR_EXCEPTIONS
