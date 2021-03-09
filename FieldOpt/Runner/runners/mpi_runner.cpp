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

#include "mpi_runner.h"
#include "Optimization/case_transfer_object.h"
#include "Model/model_synchronization_object.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/mpi/status.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/serialization/list.hpp>

#include <iostream>
#include "Utilities/verbosity.h"
#include "Utilities/printer.hpp"

BOOST_IS_MPI_DATATYPE(boost::uuids::uuid)

namespace Runner {
namespace MPI {

MPIRunner::MPIRunner(RuntimeSettings *rts) : AbstractRunner(rts) {
  if (rts->verbParams().vRUN > 3) {
    im_ = "world_.size(): " + num2str(world_.size()) + " | ";
    im_ += "world_.rank(): " + num2str(world_.rank());
    ext_info(im_, md_, cl_, rts->verbParams().lnw, 2, 2);
  }
  rank_ = world_.rank();
  simulator_delay_ = rts->simulation_delay();
}

void MPIRunner::SendMessage(Message &message) {
  std::string s;
  if (message.c != nullptr) {
    auto cto = Optimization::CaseTransferObject(message.c);
    std::ostringstream oss;
    boost::archive::text_oarchive oa(oss);
    oa << cto;
    s = oss.str();
  } else {
    s = "";
  }
  world_.send(message.destination, message.tag, s);
  im_ = "Sent a message to " + num2str(message.destination) + " with tag ";
  im_ += num2str(message.tag) + " (" + tag_to_string[message.tag] + ")";
  printMessage(im_, 2);
}

void MPIRunner::RecvMessage(Message &message) {
  Optimization::CaseTransferObject cto;

  im_ = "Waiting to receive a message with tag " + num2str(message.tag) + " (";
  im_ = tag_to_string[message.tag] + ") " + " from source " + num2str(message.source);
  printMessage(im_, 2);

  std::string s;
  mpi::status status = world_.recv(message.source, ANY_TAG, s);
  message.set_status(status);
  message.tag = status.tag();

  auto handle_received_case = [&]() mutable {
    std::istringstream iss(s);
    boost::archive::text_iarchive ia(iss);
    ia >> cto;
    message.c = cto.CreateCase();
  };

  if (message.tag == TERMINATE) {
    message.c = nullptr;
    printMessage("Received termination signal.", 2);
    return;

  } else if (message.tag == MODEL_SYNC) {
    im_ = "Received message with MODEL_SYNC tag. RecvMessage ";
    im_ += "method cannot handle this. Throwing exception.";
    printMessage(im_);

    em_ = "RecvMessage is unable to handle model synchronization objects. ";
    em_ += "This should be handled by the RecvModelSynchronizationObject method.";
    E(em_, md_, cl_);

  } else if (message.tag == CASE_UNEVAL) {
    printMessage("Received an unevaluated case.", 2);
    handle_received_case();

  } else if (message.tag == CASE_EVAL_SUCCESS) {
    printMessage("Received a successfully evaluated case.", 2);
    handle_received_case();

  } else if (message.tag == CASE_EVAL_INVALID) {
    printMessage("Received an invalid case.", 2);
    handle_received_case();

  } else if (message.tag == CASE_EVAL_TIMEOUT) {
    printMessage("Received a case that was terminated due to timeout.", 2);

  } else {
    printMessage("Received message with an unrecognized tag. Throwing exception.");
    E("RecvMessage received a message with an unrecognized tag.", md_, cl_);
  }
}

void MPIRunner::BroadcastModel() {
  if (rank() != 0) {
    E("BroadcastModel should only be called on the root process.", md_, cl_);
  }
  auto mso = Model::ModelSynchronizationObject(model_);
  std::ostringstream oss;
  boost::archive::text_oarchive oa(oss);
  oa << mso;
  std::string s = oss.str();
  for (int r = 1; r < world_.size(); ++r) {
    world_.send(r, MODEL_SYNC, s);
  }
}

void MPIRunner::RecvModelSynchronizationObject() {
  if (rank() == 0) {
    E("RecvModelSynchronizationObject should not be called on the root process.", md_, cl_);
  }
  Model::ModelSynchronizationObject mso;
  std::string s;
  world_.recv(0, MODEL_SYNC, s);
  std::istringstream iss(s);
  boost::archive::text_iarchive ia(iss);
  ia >> mso;
  mso.UpdateVariablePropertyIds(model_);
}

int MPIRunner::SimulatorDelay() const {
  return simulator_delay_;
}

void MPIRunner::printMessage(std::string message, int min_verb) {
  if (vp_.vRUN >= min_verb) {
    std::string time_stamp = QDateTime::currentDateTime().toString("hh:mm").toStdString();
    std::stringstream ss;
    ss << "MPI Rank " << world_.rank() << "@" << time_stamp << ": " << message;
    Printer::info(ss.str());
  }
}
}
}
