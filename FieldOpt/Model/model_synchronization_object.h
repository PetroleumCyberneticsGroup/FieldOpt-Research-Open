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

#ifndef FIELDOPT_MODEL_SYNCHRONIZATION_OBJECT_H
#define FIELDOPT_MODEL_SYNCHRONIZATION_OBJECT_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_serialize.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <map>
#include <list>
#include <string>

#include "model.h"

using namespace boost::uuids;
using namespace std;

namespace Model {
class Model;

/*!
 * @brief The ModelSynchronization class is used by MPI-based
 * runners to synchronize the variable IDs across process
 * instances. This is necessary to make the transferred Case
 * objects work.
 */
class ModelSynchronizationObject {
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & binary_variable_ids_;
    ar & discrete_variable_ids_;
    ar & continous_variable_ids_;
  }

 public:
  ModelSynchronizationObject() {}
  ModelSynchronizationObject(Model *model);

  // QMap<QString, QUuid> GetDiscreteVariableMap();
  // QMap<QString, QUuid> GetContinousVariableMap();
  // QMap<QString, QUuid> GetBinaryVariableMap();
  QList<QPair<QString, QUuid>> GetDiscreteVariableMap();
  QList<QPair<QString, QUuid>> GetContinousVariableMap();
  QList<QPair<QString, QUuid>> GetBinaryVariableMap();

  void UpdateVariablePropertyIds(Model *model);

 private:
  //!< Mapping from variable name to variable UUID, for discrete variables.
  std::list<std::pair<string, uuid>> discrete_variable_ids_;

  //!< Mapping from variable name to variable UUID, for continuous variables.
  std::list<std::pair<string, uuid>> continous_variable_ids_;

  //!< Mapping from variable name to variable UUID, for binary variables.
  // std::map<string, uuid> binary_variable_ids_;
  std::list<std::pair<string, uuid>> binary_variable_ids_;

  template<typename T>
  std::list<std::pair<string, uuid>> createNameToIdMapping(const QList<QPair<QUuid, T>> *qhash) const;

  //!< Create a standard library hash map from a QHash
  // std::map<string, uuid> createNameToIdMapping(const QMap<QUuid, Properties::DiscreteProperty *> *qhash) const;
  std::list<std::pair<string, uuid>> createNameToIdMapping(const QList<QPair<QUuid, Properties::DiscreteProperty *>> *qhash) const;

  //!< Create a standard library hash map from a QHash
  // std::map<string, uuid> createNameToIdMapping(const QMap<QUuid, Properties::ContinuousProperty *> *qhash) const;
  std::list<std::pair<string, uuid>> createNameToIdMapping(const QList<QPair<QUuid, Properties::ContinuousProperty *>> *qhash) const;

  //!< Create a standard library hash map from a QHash
  // std::map<string, uuid> createNameToIdMapping(const QMap<QUuid, Properties::BinaryProperty *> *qhash) const;
  std::list<std::pair<string, uuid>> createNameToIdMapping(const QList<QPair<QUuid, Properties::BinaryProperty *>> *qhash) const;

  //!< Convert a std/boost based mapping to a Qt based mapping to be used by the rest of the model.
  // QMap<QString, QUuid> convertToQtMapping(const std::map<string, uuid> map);
  QList<QPair<QString, QUuid>> convertToQtMapping(const std::list<std::pair<string, uuid>> map);

  //!< Create a QUuid from a boost uuid
  QUuid boostUuidToQuuid(const uuid buuid) const;

  //!< Create a boost uuid from a QUuid
  uuid qUuidToBoostUuid(const QUuid quuid) const;

  //!< Get the string representation of a boost uuid
  QString boostUuidToQstring(const uuid buuid) const;
};

}


#endif //FIELDOPT_MODEL_SYNCHRONIZATION_OBJECT_H
