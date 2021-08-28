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

#include <boost/lexical_cast.hpp>
#include <iostream>
#include "model_synchronization_object.h"

namespace Model {

using namespace std;
using namespace boost::uuids;

ModelSynchronizationObject::ModelSynchronizationObject(Model *model) {
  discrete_variable_ids_ = createNameToIdMapping(model->variables()->GetDiscreteVariables());
  continous_variable_ids_ = createNameToIdMapping(model->variables()->GetContinuousVariables());
  binary_variable_ids_ = createNameToIdMapping(model->variables()->GetBinaryVariables());
}

template<typename T>
list<pair<string, uuid>>
ModelSynchronizationObject::createNameToIdMapping(const QList<QPair<QUuid, T>> *qhash) const {
  auto name_id_map = list<pair<string, uuid>>();
  for (auto prop : *qhash) {
    pair<string, uuid> pval = pair<string, uuid>(prop.second->name().toStdString(), qUuidToBoostUuid(prop.first));
    name_id_map.push_back(pval);
  }
  return name_id_map;
}

// map<string, uuid>
// ModelSynchronizationObject::createNameToIdMapping(const QMap<QUuid, Properties::DiscreteProperty *> *qhash) const {
//   auto name_id_map = map<string, uuid>();
//   for (auto prop : qhash->values())
//     name_id_map[prop->name().toStdString()] = qUuidToBoostUuid(prop->id());
//   return name_id_map;
// }

list<pair<string, uuid>>
ModelSynchronizationObject::createNameToIdMapping(const QList<QPair<QUuid, Properties::DiscreteProperty *>> *qhash) const {
  auto name_id_map = list<pair<string, uuid>>();
  for (auto prop : *qhash) {
    pair<string, uuid> pval = pair<string, uuid>(prop.second->name().toStdString(), qUuidToBoostUuid(prop.first));
    name_id_map.push_back(pval);
  }
  return name_id_map;
}

// map<string, uuid>
// ModelSynchronizationObject::createNameToIdMapping(const QMap<QUuid, Properties::ContinuousProperty *> *qhash) const {
//   auto name_id_map = map<string, uuid>();
//   for (auto prop : qhash->values())
//     name_id_map[prop->name().toStdString()] = qUuidToBoostUuid(prop->id());
//   return name_id_map;
// }

list<pair<string, uuid>>
ModelSynchronizationObject::createNameToIdMapping(const QList<QPair<QUuid, Properties::ContinuousProperty *>> *qhash) const {
  auto name_id_map = list<pair<string, uuid>>();
  for (auto prop : *qhash) {
    pair<string, uuid> pval = pair<string, uuid>(prop.second->name().toStdString(), qUuidToBoostUuid(prop.first));
    name_id_map.push_back(pval);
  }
  return name_id_map;
}

// map<string, uuid>
// ModelSynchronizationObject::createNameToIdMapping(const QMap<QUuid, Properties::BinaryProperty *> *qhash) const {
//   auto name_id_map = map<string, uuid>();
//   for (auto prop : qhash->values())
//     name_id_map[prop->name().toStdString()] = qUuidToBoostUuid(prop->id());
//   return name_id_map;
// }

list<pair<string, uuid>>
ModelSynchronizationObject::createNameToIdMapping(const QList<QPair<QUuid, Properties::BinaryProperty *>> *qhash) const {
  auto name_id_map = list<pair<string, uuid>>();
  for (auto prop : *qhash) {
    pair<string, uuid> pval = pair<string, uuid>(prop.second->name().toStdString(), qUuidToBoostUuid(prop.first));
    name_id_map.push_back(pval);
  }
  return name_id_map;
}

// QMap<QString, QUuid> ModelSynchronizationObject::convertToQtMapping(const std::map<string, uuid> map) {
//   QMap<QString, QUuid> qmap = QMap<QString, QUuid>();
//   for (auto const &ent : map) {
//     qmap[QString::fromStdString((ent.first))] = boostUuidToQuuid(ent.second);
//   }
//   return qmap;
// }

QList<QPair<QString, QUuid>> ModelSynchronizationObject::convertToQtMapping(const std::list<pair<string, uuid>> map) {
  QList<QPair<QString, QUuid>> qmap = QList<QPair<QString, QUuid>>();
  for (auto const &ent : map) {
    QPair<QString, QUuid> pval = QPair<QString, QUuid>(QString::fromStdString((ent.first)), boostUuidToQuuid(ent.second));
    qmap.append(pval);
  }
  return qmap;
}

// QMap<QString, QUuid> ModelSynchronizationObject::GetDiscreteVariableMap() {
//   return convertToQtMapping(discrete_variable_ids_);
// }

QList<QPair<QString, QUuid>> ModelSynchronizationObject::GetDiscreteVariableMap() {
  return convertToQtMapping(discrete_variable_ids_);
}

// QMap<QString, QUuid> ModelSynchronizationObject::GetContinousVariableMap() {
//   return convertToQtMapping(continous_variable_ids_);
// }

QList<QPair<QString, QUuid>> ModelSynchronizationObject::GetContinousVariableMap() {
  return convertToQtMapping(continous_variable_ids_);
}

// QMap<QString, QUuid> ModelSynchronizationObject::GetBinaryVariableMap() {
//   return convertToQtMapping(binary_variable_ids_);
// }

QList<QPair<QString, QUuid>> ModelSynchronizationObject::GetBinaryVariableMap() {
  return convertToQtMapping(binary_variable_ids_);
}

void ModelSynchronizationObject::UpdateVariablePropertyIds(Model *model) {
  auto vpc = model->var_container_;

  // Binary variables
  if (vpc->BinaryVariableSize() > 0) {
    // for (auto name : GetBinaryVariableMap().keys()) {
    for (int ii=0; ii < GetBinaryVariableMap().size(); ii++) {
      auto var = GetBinaryVariableMap().at(ii);
      auto name = var.first;
      auto prop = vpc->GetBinaryVariable(name);
      // QUuid new_uuid = GetBinaryVariableMap()[name];
      QUuid new_uuid = var.second;
      QUuid old_uuid = prop->id();
      prop->UpdateId(new_uuid);
      // vpc->binary_variables_->insert(new_uuid, prop);
      // vpc->binary_variables_->remove(old_uuid);
      vpc->binary_variables_->append(qMakePair(new_uuid, prop));
      vpc->binary_variables_->removeOne(qMakePair(old_uuid, prop));
    }
  }

  // Continuous variables
  if (vpc->ContinuousVariableSize() > 0) {
    for (int ii=0; ii < GetContinousVariableMap().size(); ii++) {
      auto var = GetContinousVariableMap().at(ii);
      auto name = var.first;
      auto prop = vpc->GetContinousVariable(name);
      // QUuid new_uuid = GetContinousVariableMap()[name];
      QUuid new_uuid = var.second;
      QUuid old_uuid = prop->id();
      prop->UpdateId(new_uuid);
      // vpc->continuous_variables_->insert(new_uuid, prop);
      // vpc->continuous_variables_->remove(old_uuid);
      vpc->continuous_variables_->append(qMakePair(new_uuid, prop));
      vpc->continuous_variables_->removeOne(qMakePair(old_uuid, prop));
    }
  }

  // Discrete variables
  if (vpc->DiscreteVariableSize() > 0) {
    for (int ii=0; ii < GetDiscreteVariableMap().size(); ii++) {
      auto var = GetDiscreteVariableMap().at(ii);
      auto name = var.first;
      auto prop = vpc->GetDiscreteVariable(name);
      // QUuid new_uuid = GetDiscreteVariableMap()[name];
      QUuid new_uuid = var.second;
      QUuid old_uuid = prop->id();
      prop->UpdateId(new_uuid);
      // vpc->discrete_variables_->insert(new_uuid, prop);
      // vpc->discrete_variables_->remove(old_uuid);
      vpc->discrete_variables_->append(qMakePair(new_uuid, prop));
      vpc->discrete_variables_->removeOne(qMakePair(old_uuid, prop));
    }
  }
}

QUuid ModelSynchronizationObject::boostUuidToQuuid(const uuid buuid) const {
  return QUuid(boostUuidToQstring(buuid));
}

uuid ModelSynchronizationObject::qUuidToBoostUuid(const QUuid quuid) const {
  uuid buuid(string_generator()(quuid.toString().toStdString()));
  return buuid;
}

QString ModelSynchronizationObject::boostUuidToQstring(const uuid buuid) const {
  const std::string uuid_string = boost::lexical_cast<std::string>(buuid);
  QString uuid_qstring = QString(uuid_string.c_str());
  uuid_qstring.append("}");
  uuid_qstring.prepend("{");
  return uuid_qstring;
}


}
