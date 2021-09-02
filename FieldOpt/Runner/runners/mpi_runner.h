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

#ifndef FIELDOPT_MPIRUNNER_H
#define FIELDOPT_MPIRUNNER_H

#include "abstract_runner.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

namespace mpi = boost::mpi;

namespace Runner {
namespace MPI {
class Worker;
class Overseer;
}
}

namespace Runner {
namespace MPI {
/*!
 * @brief The MPIRunner class is the base class
 * for runners using MPI for parallelization.
 */
class MPIRunner : public AbstractRunner {

  friend class Worker;
  friend class Overseer;

 public:
  mpi::environment &env() { return env_; }

  mpi::communicator &world() { return world_; }

  int rank() const { return rank_; }

  /*!
   * @brief Tags used when sending and receiving.
   *
   * CASE_UNEVAL: Used when sending unevaluated cases.
   *
   * CASE_EVAL_SUCCESS:
   * Used when sending successfully evaluated cases.
   *
   * CASE_EVAL_INVALID:
   * Used when sending cases that were deemed invalid.
   *
   * CASE_EVAL_TIMEOUT:
   * Used when sending cases whose simulation was terminated
   * by a timeout condition.
   *
   * MODEL_SYNC:
   * Used when sending model synchronization objects.
   *
   * ANY_TAG: This will match any tag.
   *
   * TERMINATE: Sent by overseer to terminate a worker.
   */
  enum MsgTag : int {
    CASE_UNEVAL = 1, CASE_EVAL_SUCCESS = 2, CASE_EVAL_INVALID = 3,
    CASE_EVAL_TIMEOUT = 4, MODEL_SYNC = 10, TERMINATE = 100,
    ANY_TAG = MPI_ANY_TAG
  };

  std::map<int, std::string> tag_to_string = {
    {-1, "any tag"},
    {1, "unevaluated case"},
    {2, "successfully evaluated case"},
    {3, "invalid case"},
    {4, "timed out case"},
    {10, "model synchronization object"},
    {100, "termination signal"}
  };

  /*!
   * @brief The Message struct should be used when sending
   * and receiving any message.
   */
  struct Message {
    Message() {
      c = nullptr;
      this->tag = MPI_ANY_TAG;
      this->source = MPI_ANY_SOURCE;
      this->destination = MPI_ANY_SOURCE;
    }

    void set_status(mpi::status stat) {
      this->status = stat;
      this->source = stat.source();
      this->tag = stat.tag();
    }

    MsgTag get_tag() {
      switch (tag) {
        case 1: return CASE_UNEVAL;
        case 2: return CASE_EVAL_SUCCESS;
        case 3: return CASE_EVAL_INVALID;
        case 4: return CASE_EVAL_TIMEOUT;
        case 10: return MODEL_SYNC;
        case 100: return TERMINATE;
      }
      return ANY_TAG;
    }

    //!< Case associated with the message (if any).
    Optimization::Case *c;

    int tag; //!< Tag for the message.

    int source; //!< Rank of process sending the message.

    int destination; //!< Rank of process receiving message.

    mpi::status status; //!< Status object for message.
  };

  /*!
   * @brief Send a message potentially containing a case.
   * @param message The message to be sent.
   */
  void SendMessage(Message &message);


  /*!
   * @brief Receive a message potentially containing a Case.
   *
   * If source or tag is specified in the message parameter,
   * only messages with this tag and/or source will be received.
   *
   * If not, a message will be received from any source and/or
   * with any tag, and the values will be entered in the Message
   * object.
   *
   * If a case is received, the c field in the parameter
   * message object will be set to it.
   * @param message
   * @return
   */
  void RecvMessage(Message &message);

  /*!
   * @brief Create a ModelSynchronizationObject
   * and send it to all other processes.
   *
   * This should be called by the process with rank 0, to
   * make the variable UUIDs match across all processes.
   */
  void BroadcastModel();

  /*!
   * @brief Receive the ModelSynchronizationObject broadcast
   * by the root process. This is applied to the model object.
   */
  void RecvModelSynchronizationObject();

  int SimulatorDelay() const;

 protected:
  explicit MPIRunner(RuntimeSettings *rts);

  mpi::environment env_;
  mpi::communicator world_;

  int rank_;
  int scheduler_rank_ = 0;
  int simulator_delay_;

  /*!
   * @brief Print a message to the console.
   * @param message The message to be printed.
   * @param min_verb Minimum verbosity level requred
   * for the message to be printed.
   */
  void printMessage(std::string message, int min_verb=1);

  string im_ = "", wm_ = "", em_ = "";
  string md_ = "Runner";
  string cl_ = "MPIRunner";

};
}
}


#endif //FIELDOPT_MPIRUNNER_H
