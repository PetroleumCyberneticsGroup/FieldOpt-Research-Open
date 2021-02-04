/***********************************************************
Created: 22.09.2015 2015 by einar

Copyright (C) 2015-2015
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

#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <QPair>
#include <QList>
#include "Settings/model.h"
#include "Model/model.h"


namespace Optimization {
namespace Objective {

/*!
 * \brief The Objective class defines an interface for
 * objective functions. It cannot be instantiated on its
 * own.
 */
class Objective
{
 public:
  /*!
   * \brief value Get evaluated value for objective function.
   */
  virtual double value() const = 0;

  virtual double value(bool base_case) = 0;

 protected:
  Objective(Settings::Optimizer *settings);

  Settings::VerbParams vp_;
  string md_ = "Optimization";
  string cl_ = "Objective";
};

}
}

#endif // OBJECTIVE_H
