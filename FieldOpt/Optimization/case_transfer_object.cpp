/***********************************************************
Copyright (C) 2015-2017
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

#include <boost/lexical_cast.hpp>
#include "case_transfer_object.h"
#include <QString>

namespace Optimization {
using namespace boost::uuids;
using namespace std;


CaseTransferObject::CaseTransferObject(Optimization::Case *c) {
  id_ = qUuidToBoostUuid(c->id_);
  objective_function_value_ = c->objective_function_value_;

  // binary_variables_ = qHashToStdMap(c->binary_variables_);
  // integer_variables_ = qHashToStdMap(c->integer_variables_);
  // real_variables_ = qHashToStdMap(c->real_variables_);

  // binary_variables_ = qMapToStdMap(c->binary_variables_);
  // integer_variables_ = qMapToStdMap(c->integer_variables_);
  // real_variables_ = qMapToStdMap(c->real_variables_);

  binary_variables_ = qListToStdList(c->binary_variables_);
  integer_variables_ = qListToStdList(c->integer_variables_);
  real_variables_ = qListToStdList(c->real_variables_);

  wic_time_secs_ = c->GetWICTime();
  sim_time_secs_ = c->GetSimTime();
  ensemble_realization_ = c->GetEnsembleRealization().toStdString();

  status_eval_ = c->state.eval;
  status_cons_ = c->state.cons;
  status_queue_ = c->state.queue;
  status_err_msg_ = c->state.err_msg;
}

Case *CaseTransferObject::CreateCase() {
  auto c = new Case();

  // c->binary_variables_ = stdMapToQhash(binary_variables_);
  // c->integer_variables_ = stdMapToQhash(integer_variables_);
  // c->real_variables_ = stdMapToQhash(real_variables_);

  // c->binary_variables_ = stdMapToQmap(binary_variables_);
  // c->integer_variables_ = stdMapToQmap(integer_variables_);
  // c->real_variables_ = stdMapToQmap(real_variables_);

  c->binary_variables_ = stdListToQlist(binary_variables_);
  c->integer_variables_ = stdListToQlist(integer_variables_);
  c->real_variables_ = stdListToQlist(real_variables_);

  c->id_ = boostUuidToQuuid(id_);
  c->objective_function_value_ = objective_function_value_;
  c->SetWICTime(wic_time_secs_);
  c->SetSimTime(sim_time_secs_);
  c->SetEnsembleRealization(QString::fromStdString(ensemble_realization_));
  c->state.eval = static_cast<Case::CaseState::EvalStatus>(status_eval_);
  c->state.cons = static_cast<Case::CaseState::ConsStatus>(status_cons_);
  c->state.queue = static_cast<Case::CaseState::QueueStatus>(status_queue_);
  c->state.err_msg = static_cast<Case::CaseState::ErrorMessage >(status_err_msg_);
  return c;
}

QUuid CaseTransferObject::boostUuidToQuuid(const uuid buuid) const {
  return QUuid(boostUuidToQstring(buuid));
}

uuid CaseTransferObject::qUuidToBoostUuid(const QUuid quuid) const {
  uuid uuid(string_generator()(quuid.toString().toStdString()));
  return uuid;
}

QString CaseTransferObject::boostUuidToQstring(const uuid buuid) const {
  auto uuid_string = boost::lexical_cast<string>(buuid);
  auto uuid_qstring = QString(uuid_string.c_str());
  uuid_qstring.append("}");
  uuid_qstring.prepend("{");
  return uuid_qstring;
}

// template<typename T>
// std::map<uuid, T> CaseTransferObject::qHashToStdMap(const QHash<QUuid, T> &qhash) const {
//   std::map<uuid, T> stdmap = std::map<uuid, T>();
//   for (auto key : qhash.keys()) {
//     stdmap[qUuidToBoostUuid(key)] = qhash[key];
//   }
//   return stdmap;
// }

// template<typename T>
// std::map<uuid, T> CaseTransferObject::qMapToStdMap(const QMap<QUuid, T> &qhash) const {
//   std::map<uuid, T> stdmap = std::map<uuid, T>();
//   for (auto key : qhash.keys()) {
//     stdmap[qUuidToBoostUuid(key)] = qhash[key];
//   }
//   return stdmap;
// }

template<typename T>
std::list<pair<uuid, T>> CaseTransferObject::qListToStdList(const QList<QPair<QUuid, T>> &qhash) const {
  std::list<pair<uuid, T>> stdlist = std::list<pair<uuid, T>>();
  for (int ii=0; ii < qhash.size(); ii++) {
    pair<uuid, T> qpair = pair<uuid, T>(qUuidToBoostUuid(qhash.at(ii).first), qhash.at(ii).second);
    stdlist.push_back(qpair);
  }
}

// template<typename T>
// QHash<QUuid, T> CaseTransferObject::stdMapToQhash(const std::map<uuid, T> &stmap) const {
//   QHash<QUuid, T> qhash = QHash<QUuid, T>();
//   for (auto const &ent : stmap) {
//     qhash[boostUuidToQuuid(ent.first)] = ent.second;
//   }
//   return qhash;
// }

// template<typename T>
// QMap<QUuid, T> CaseTransferObject::stdMapToQmap(const std::map<uuid, T> &stmap) const {
//   QMap<QUuid, T> qhash = QMap<QUuid, T>();
//   for (auto const &ent : stmap) {
//     qhash[boostUuidToQuuid(ent.first)] = ent.second;
//   }
//   return qhash;
// }

template<typename T>
QList<QPair<QUuid, T>> CaseTransferObject::stdListToQlist(const std::list<pair<uuid, T>> &stmap) const {
  QList<QPair<QUuid, T>> qhash = QList<QPair<QUuid, T>>();
  for (auto const &ent : stmap) {
    QPair<QUuid, T> qpair = QPair<QUuid, T>(boostUuidToQuuid(ent.first), ent.second);
    qhash.append(qpair);
  }
  return qhash;
}

}
