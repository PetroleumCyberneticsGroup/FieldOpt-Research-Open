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

#ifndef ECLSIMULATOR_H
#define ECLSIMULATOR_H

#include "simulator.h"
#include "driver_file_writers/ecldriverfilewriter.h"
#include "Model/model.h"
#include <QStringList>

namespace Simulation {

// =========================================================
/*!
 * \brief The ECLSimulator class implements simulation
 * of models using the ECLIPSE reservoir simulator.
 *
 * This class should not be used directly except for
 * instantiation. All other actions should be called
 * through the Simulator class. The intended use is
 * as follows:
 *
 * \code
 *  Simulator sim = new ECLSimulator();
 *  sim.SetOutputDirectory("some/path");
 *  sim.Evaluate();
 *  sim.CleanUp();
 * \endcode
 *
 * \todo Support custom execution commands.
 */
class ECLSimulator : public Simulator {
 public:
  ECLSimulator(Settings::Settings *settings, Model::Model *model);
  
  /*!
   * \brief Evaluate Executes the simulation of
   * the current model. The evaluation is blocking.
   */
  void Evaluate() override;
  bool Evaluate(int timeout, int threads = 1) override;
  
  bool Evaluate(const Settings::Ensemble::Realization &realization,
      int timeout, int threads = 1) override;
  
  void WriteDriverFilesOnly() override;
  /*!
   * \brief CleanUp Deletes files created during
   * the simulation. All files except the .DATA,
   * .UNSMRY, .SMSPEC  and .LOG are deleted.
   */
  void CleanUp() override;
 
 private:

  //!< Driver file name without the final .DATA
  QString deck_name_;
  Settings::Settings *settings_;
  void copyDriverFiles();
  
  // Simulator interface
 protected:
  void UpdateFilePaths() override;
};

}
#endif // ECLSIMULATOR_H
