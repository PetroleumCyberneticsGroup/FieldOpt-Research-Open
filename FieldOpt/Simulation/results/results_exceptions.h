/***********************************************************
Created: 12.10.2015 2015 by einar

Copyright (C) 2015
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#ifndef RESULTS_EXCEPTIONS
#define RESULTS_EXCEPTIONS

#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/lexical_cast.hpp>

using std::cout;
using std::endl;
using std::string;
using std::runtime_error;

namespace Simulation {
namespace Results {

class RsltFileNotFoundExc : public std::runtime_error {
 public:
  string em_;
  explicit RsltFileNotFoundExc(const string& path)
    : runtime_error("[ECLResults] No valid result file found") {
    em_ = "Path: " + path + "\n";
    em_ += "Enter dir and run simulation to create results files, ";
    em_ += "e.g., run @e300 -file OLPZ_BCXX_R37_F37_W01 -local";
    cout << endl << em_ << endl;
  }
};

class RsltsNotAvailExc : public std::runtime_error {
 public:
  string em_;
  explicit RsltsNotAvailExc(const string& err_src)
    : runtime_error("[ECLResults] Results not available.") {
    em_ = "Check results in " + err_src + "." + "\n";
    cout << endl << em_ << endl;
  }
};

class RsltPropKeyDoesNotExistExc : public std::runtime_error {
 public:
  string em_;
  explicit RsltPropKeyDoesNotExistExc(const string &simulator,
                                      const string &prop)
    : runtime_error("[ECLResults] Requested property key cannot be retrieved.") {
    em_ = "Key: " + prop + " cannot be retrieved from " + simulator + " summaries.";
    cout << endl << em_ << endl;
  }
};

class RsltTimeIndexInvalidExc : public std::runtime_error {
 public:
  string em_;
  explicit RsltTimeIndexInvalidExc(const int index)
    : runtime_error("[ECLResults] Time index invalid.") {
    em_ = "Time index " + std::to_string(index) + " is not";
    em_ += " valid for the current results.";
    cout << endl << em_ << endl;
  }
};

}}

#endif // RESULTS_EXCEPTIONS
