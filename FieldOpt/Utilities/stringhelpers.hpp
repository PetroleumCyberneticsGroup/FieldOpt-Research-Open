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

#ifndef STRINGHELPERS_FUNCTIONS_H
#define STRINGHELPERS_FUNCTIONS_H

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <Eigen/Core>
#include <iomanip>

using namespace std;

template <typename T>
inline string vec_to_str(vector<T> vec) {
  stringstream str;
  str << std::fixed << scientific << std::setprecision(3);
  if (vec.size() == 0)
    return "";
  str << vec[0];
  if (vec.size() > 1) {
    for (int i = 1; i < vec.size(); ++i) {
      str << ", " << vec[i];
    }
  }
  return str.str();
}

inline string eigenvec_to_str(Eigen::VectorXd vec) {
  string str = boost::lexical_cast<string>(vec(0));
  for (int i = 1; i < vec.size(); ++i) {
    str = str + ", " + boost::lexical_cast<string>(vec(i));
  }
  return str;
}

#endif // STRINGHELPERS_FUNCTIONS_H
