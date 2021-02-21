/***********************************************************
Copyright (C) 2015
Einar J.M. Baumann <einar.baumann@gmail.com>
Created: 30.11.2015 2015 by einar

Modified 2021 Mathias Bellout
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

#ifndef OPTIMIZATION_EXCEPTIONS
#define OPTIMIZATION_EXCEPTIONS

#include <stdexcept>
#include <string>
#include <QString>

namespace Optimization {

class ObjectiveFunctionException : public std::runtime_error {
 public:
  explicit ObjectiveFunctionException(const string &message)
    : std::runtime_error(message) {}
};

class VariableException : public std::runtime_error {
 public:
  explicit VariableException(const string &message)
    : std::runtime_error(message) {}
};

class CaseHandlerException : public std::runtime_error {
 public:
  explicit CaseHandlerException(const string &message)
    : std::runtime_error(message) {}
};

class OptimizerInitializationException : public std::runtime_error {
 public:
  explicit OptimizerInitializationException(const string &message)
    : std::runtime_error("Unable to initialize the optimizer: " + message) {}
};

}

#endif // OPTIMIZATION_EXCEPTIONS

