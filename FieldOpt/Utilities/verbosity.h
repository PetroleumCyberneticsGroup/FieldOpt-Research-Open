/***********************************************************
Created by einar on 8/24/18.
Copyright (C) 2018
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

#ifndef FIELDOPT_VERBOSITY_H
#define FIELDOPT_VERBOSITY_H

/*!
 * @brief Verbosity level settings for all FieldOpt modules.
 * FieldOpt needs to be recompiled if settings are changed.
 *
 * Level descriptions (they are additive):
 *  0 - Silent (only errors)
 *  1 - Only warnings and rudimentary progression
 *      (e.g, once pr. iteration).
 *  2 - More progression printing.
 *  3 - Debug-level printing. Not recommended to use
 *      on more than one module at a time.
 */

//                    Module
#define VERB_MOD 1 // Model
#define VERB_OPT 1 // Optimization
#define VERB_WIC 1 // WellIndexCalculation
#define VERB_SIM 1 // Simulation
#define VERB_RUN 1 // Runner
#define VERB_RES 1 // Reservoir
#define VERB_SET 1 // Settings

#define LINEWDTH 165

namespace Settings {

struct VerbParams {
  int lnw = 141;
  int vMOD = 0;
  int vOPT = 3;
  int vWIC = 0;
  int vSIM = 0;
  int vRUN = 0;
  int vRES = 0;
  int vSET = 0;
  int vUTI = 0;
};

}

#endif //FIELDOPT_VERBOSITY_H
