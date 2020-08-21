/***********************************************************
Created: 12.10.2015 2015 by einar

Copyright (C) 2015
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

#ifndef RESULTS_EXCEPTIONS
#define RESULTS_EXCEPTIONS

#include <stdexcept>
#include <string>
#include <boost/lexical_cast.hpp>

using std::string;
using std::runtime_error;

namespace Simulation { namespace Results {

class RsltFileNotFoundExc : public std::runtime_error {
 public:
  RsltFileNotFoundExc(const string &path)
    : runtime_error("No valid result file found at path" + path) {}
};

class RsltsNotAvailExc : public std::runtime_error {
 public:
  RsltsNotAvailExc()
    : runtime_error("Results are not currently available.") {}
};

class RsltPropKeyDoesNotExistExc : public std::runtime_error {
 public:
  explicit RsltPropKeyDoesNotExistExc(const string &simulator,
                                      const string &prop)
    : runtime_error("Requested property key (" + prop
                      + ") cannot be retrieved from " + simulator
                      + "summaries.") {}
};

class ResultTimeIndexInvalidException : public std::runtime_error {
 public:
  ResultTimeIndexInvalidException(const int index)
    : runtime_error("Time index "
    + boost::lexical_cast<std::string>(index)
      + " is not valid for the current results.") {}
};

}}

#endif // RESULTS_EXCEPTIONS
