/***********************************************************
Copyright (C) 2021
Mathias Bellout <chakibbb-pcg@gmail.com>

bellout - Tue Aug 31 2021 15:22:17 week 35 CET+0200

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

#ifndef FIELDOPT_RUNNER_RUNNERS_BILVL_SYNCHR_MPI_RNR_H_
#define FIELDOPT_RUNNER_RUNNERS_BILVL_SYNCHR_MPI_RNR_H_

#include "mpi_runner.h"
#include "overseer.h"
#include "worker.h"

namespace Runner {
namespace MPI {

class BilevelSynchrMPIRunner : public MPIRunner, public Loggable {
 public:
  LogTarget GetLogTarget() override;
  explicit BilevelSynchrMPIRunner(RuntimeSettings *rts);

  map<string, string> GetState() override;
  QUuid GetId() override;
  void Execute() override;

 private:
  MPI::Overseer *overseer_;
  MPI::Worker *worker_;

  bool model_update_done_ = false;
  bool simulation_done_ = false;

  double bl_ps_fdiff_ = 0.0;

  void initialDistribution();
};

}
}

#endif //FIELDOPT_RUNNER_RUNNERS_BILVL_SYNCHR_MPI_RNR_H_
