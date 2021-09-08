/***********************************************************
Copyright (C) 2016
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020 Mathias Bellout
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

#ifndef BOOKKEEPER_H
#define BOOKKEEPER_H

#include "Settings/settings.h"
#include "Optimization/case_handler.h"

namespace Runner {

/*!
 * \brief The Bookkeeper class is used to check whether a
 * case has already been evaluated, or is currently being
 * evaluated.
 *
 * If a case has been evaluated, the Bookkeeper can set it's
 * objective function value to the already known value.
 *
 * The Bookkeeper uses the case_handler from the optimizer
 * to keep track of which cases have been evaluated.
 *
 * \todo Handle the case where a case is currently being
 * evaluated; i.e. there exists a case in the "under
 * evaluation" list which is equal to the case being
 * checked, but has a different UUID.
 */
class Bookkeeper
{
 public:
  Bookkeeper(Settings::Settings *settings,
             Optimization::CaseHandler *case_handler);

  /*!
   * \brief IsEvaluated Check if a case has already been evaluated. If the set_obj parameter
   * is set to true, the objective value of the case will be set to that of the existing
   * case.
   * \param c The case to check.
   * \param set_obj Automatically set the objective value if it is known.
   * \return True if objective value of these variable values has already been calculated; otherwise false.
   */
  bool IsEvaluated(Optimization::Case *c, bool set_obj=false);

 private:
  double tol_;
  Optimization::CaseHandler *case_handler_;
};

}

#endif // BOOKKEEPER_H
