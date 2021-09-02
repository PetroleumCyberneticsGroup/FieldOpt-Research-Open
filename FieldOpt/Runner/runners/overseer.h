/***********************************************************
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

#ifndef FIELDOPT_OVERSEER_H
#define FIELDOPT_OVERSEER_H

#include "mpi_runner.h"
#include "Utilities/time.hpp"
#include <chrono>

namespace Runner {
namespace MPI {
/*!
 * @brief The Overseer class takes care of distributing cases
 * between workers. The runner taken as an is primarily used
 * for the common MPI helpers.
 */
class Overseer {
 public:
  explicit Overseer(MPIRunner *runner);

  /*!
   * @brief Assign a Case to a Worker. If no workers
   * are free, the scheduler will wait untill one is.
   * @param c The case to be assigned to a Worker for evaluation.
   * @param preferred_worker Prefer to assign the case to
   * the worker with this rank (default: no preference)
   */
  void AssignCase(Optimization::Case *c, int preferred_worker=-1);

  /*!
   * @brief Wait to receive an evaluated case.
   * @return An evaluated case object.
   */
  Optimization::Case *RecvEvaluatedCase();

  /*!
   * @brief Wait for a message with the TERMINATE tag from
   * each of the workers to confirm termination before moving
   * on to finalization.
   */
  void EnsureWorkerTermination();

  /*!
   * @brief Terminate all workers by sending a message with
   * the TERMINATE tag.
   */
  void TerminateWorkers();

  /*!
   * @brief Get the number of workers that are currently not
   * performing any work.
   */
  int NumberOfFreeWorkers();

  /*!
   * @brief Get the number of workers that are currently
   * executing simulations.
   */
  int NumberOfBusyWorkers();

  /*!
   * @brief Get the ranks of all free workers.
   */
  std::vector<int> GetFreeWorkerRanks() const;

  /*!
   * @brief The WorkerStatus struct holds information about
   * the current status of all workers in the network.
   */
  struct WorkerStatus {
    WorkerStatus() { rank = -1; }

    explicit WorkerStatus(int r) { rank = r;}

    //!< Rank of the process the worker is running on.
    int rank;

    //!< Indicates if the worker is currently performing simulations.
    bool working = false;

    //!< The last time a job was sent to the worker.
    QDateTime working_since;

    //!< Number of seconds since last work was sent to the process.
    int working_seconds() const {
      return time_since_seconds(working_since);
    }

    /*!
     * @brief Start the worker. Should be called whenever
     * work is sent to the worker. This marks is as working.
     */
    void start() {
      working = true;
      working_since = QDateTime::currentDateTime();
    }
    /*!
     * @brief Stop the worker. This should be called whenever
     * results are received from the worker. This marks the
     * worker as not working.
     */
    void stop() {
      working = false;
    }
  };

  /*!
       * @brief Get the status for the longest running worker.
       * @return A copy of the status object for the longest running worker.
       */
  WorkerStatus * GetLongestRunningWorker();

  //!< The message tag for the last received case.
  MPIRunner::MsgTag last_case_tag;

  void addToBestF(double f) { bestF_.push_back(f); }

  double getBestF() {
    return *std::max_element(bestF_.begin(), bestF_.end());
  }

 private:
  MPIRunner *runner_;

  std::vector<double> bestF_{0.0};

  //!< A map of the workers. The key is the rank of the process.
  QHash<int, WorkerStatus*> workers_;

  //!< Get a worker not marked as working.
  WorkerStatus * getFreeWorker();

  /*!
   * @brief Get a string summarizing the status for all workers.
   */
  std::string workerStatusSummary();

  //!< Time stamp for the start of the previous simulation.
  std::chrono::system_clock::time_point last_sim_start_;
};
}
}


#endif //FIELDOPT_OVERSEER_H
