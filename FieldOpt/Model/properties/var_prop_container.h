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

#ifndef VARIABLE_PROPERTY_CONTAINER_H
#define VARIABLE_PROPERTY_CONTAINER_H

#include <QHash>
#include <QMap>
#include <QPair>
#include <QString>

#include "property.h"
#include "binary_property.h"
#include "discrete_property.h"
#include "continous_property.h"

namespace Model {
class ModelSynchronizationObject;
}

namespace Model {
namespace Properties {

/*!
 * \brief VarPropContainer class facilitates the handling
 * of variable properties.
 *
 * This class has members that hold all continuous, discrete
 * and binary properties that are considered _variable_ in
 * the model, as well as and methods to maintain the lists.
 *
 * This class also maintains unique IDs in the form of
 * hashtable keys for every variable property of each type.
 */
class VarPropContainer
{
  friend class ::Model::ModelSynchronizationObject;
 public:
  VarPropContainer();

  //!< Add a property to the container and mark it as variable
  void AddVariable(BinaryProperty *var);

  //!< Add a property to the container and mark it as variable
  void AddVariable(DiscreteProperty *var);

  //!< Add a property to the container and mark it as variable
  void AddVariable(ContinuousProperty *var);

  //!< Get the number of binary variables.
  int BinaryVariableSize() const { return binary_variables_->size(); }

  //!< Get the number of discrete variables.
  int DiscreteVariableSize() const { return discrete_variables_->size(); }

  //!< Get the number of continuous variables.
  int ContinuousVariableSize() const { return continuous_variables_->size(); }

  //!< Get the binary variable with index id.
  BinaryProperty *GetBinaryVariable(QUuid id) const;

  //!< Get the discrete variable with index id.
  DiscreteProperty *GetDiscreteVariable(QUuid id) const;

  //!< Get the continous variable with index id.
  ContinuousProperty *GetContinuousVariable(QUuid id) const;

  //!< Get the binary variable with the specified name.
  BinaryProperty *GetBinaryVariable(QString name) const;

  //!< Get the discrete variable with the specified name.
  DiscreteProperty *GetDiscreteVariable(QString name) const;

  //!< Get the continous variable with the specified name.
  ContinuousProperty *GetContinousVariable(QString name) const;

  //!< Set the value of a binary variable.
  void SetBinaryVariableValue(QUuid id, bool val);

  //!< Set the value of a binary variable.
  void SetDiscreteVariableValue(QUuid id, int val);

  //!< Set the value of a binary variable.
  void SetContinuousVariableValue(QUuid id, double val);

  //!< Get all binary variables
  // QHash<QUuid, BinaryProperty *> *GetBinaryVariables() const { return binary_variables_; }
  // QMap<QUuid, BinaryProperty *> *GetBinaryVariables() const { return binary_variables_; }
  QList<QPair<QUuid, BinaryProperty *>> *GetBinaryVariables() const { return binary_variables_; }

  //!< Get all discrete variables
  // QHash<QUuid, DiscreteProperty *> *GetDiscreteVariables() const { return discrete_variables_; }
  // QMap<QUuid, DiscreteProperty *> *GetDiscreteVariables() const { return discrete_variables_; }
  QList<QPair<QUuid, DiscreteProperty *>> *GetDiscreteVariables() const { return discrete_variables_; }

  //!< Get all continous variables
  // QHash<QUuid, ContinuousProperty *> *GetContinuousVariables() const { return continuous_variables_; }
  // QMap<QUuid, ContinuousProperty *> *GetContinuousVariables() const { return continuous_variables_; }
  QList<QPair<QUuid, ContinuousProperty *>> *GetContinuousVariables() const { return continuous_variables_; }

  //!< Get a hashmap containing all binary variable values.
  //!< The key represents each variable's ID.
  // QHash<QUuid, bool> GetBinaryVariableValues() const;
  // QMap<QUuid, bool> GetBinaryVariableValues() const;
  QList<QPair<QUuid, bool>> GetBinaryVariableValues() const;

  //!< Get a hashmap containing all discrete variable values.
  //!< The key represents each variable's ID.
  // QHash<QUuid, int> GetDiscreteVariableValues() const;
  // QMap<QUuid, int> GetDiscreteVariableValues() const;
  QList<QPair<QUuid, int>> GetDiscreteVariableValues() const;

  //!< Get a hashmap containing all continuous variable values.
  //!< The key represents each variable's ID.
  // QHash<QUuid, double> GetContinuousVariableValues() const;
  // QMap<QUuid, double> GetContinuousVariableValues() const;
  QList<QPair<QUuid, double>> GetContinuousVariableValues() const;

  //!< Get all control (rate/bhp) variables.
  QList<ContinuousProperty *> GetWellControlVariables() const;

  //!< Get all BHP variables.
  QList<ContinuousProperty *> GetWellBHPVariables() const;

  //!< Get all BHP variables.
  QList<ContinuousProperty *> GetWellRateVariables() const;

  //!< Get all control variables for a specific well
  QList<ContinuousProperty *> GetWellControlVariables(QString well_name) const;

  //!< Get all BHP variables for a specific well.
  QList<ContinuousProperty *> GetWellBHPVariables(QString well_name) const;

  //!< Get all BHP variables for a specific well.
  QList<ContinuousProperty *> GetWellRateVariables(QString well_name) const;

  //!< Get all variables for the spline defining a well.
  QList<ContinuousProperty *> GetWellSplineVariables(QString well_name) const;

  //!< Get all variables defining a polar well spline.
  QList<ContinuousProperty *> GetPolarSplineVariables(QString well_name) const;

  //!< Get well block position variables.
  QList<DiscreteProperty *> GetWellBlockVariables() const;

  //!< Get well block position variables for a well.
  QList<DiscreteProperty *> GetWellBlockVariables(QString well_name) const;

  //!< Get x and y pseudo-continuous vertical well variables for a well. First element is x; second is y.
  QList<ContinuousProperty *> GetPseudoContVertVariables(QString well_name) const;

  //!< Get all transmissibility variables.
  QList<ContinuousProperty *> GetTransmissibilityVariables() const;

  //!< Get all transmissibility variables for a well.
  QList<ContinuousProperty *> GetTransmissibilityVariables(QString well_name) const;

  /*!
   * @brief Get a list of properties in the same order as
   * they occur in the input ID vector.
   * @param ids IDs for properties to get.
   * @return A vector of properties.
   */
  QList<ContinuousProperty *> GetContinuousProperties(QList<QUuid> ids) const;

  //!< Check that all variable names are unique. If they are not, throw an error.
  void CheckVariableNameUniqueness();

 private:
//  QHash<QUuid, BinaryProperty *> *binary_variables_;
//  QHash<QUuid, DiscreteProperty *> *discrete_variables_;
//  QHash<QUuid, ContinuousProperty *> *continuous_variables_;

//  QMap<QUuid, BinaryProperty *> *binary_variables_;
//  QMap<QUuid, DiscreteProperty *> *discrete_variables_;
//  QMap<QUuid, ContinuousProperty *> *continuous_variables_;

  QList<QPair<QUuid, BinaryProperty *>> *binary_variables_;
  QList<QPair<QUuid, DiscreteProperty *>> *discrete_variables_;
  QList<QPair<QUuid, ContinuousProperty *>> *continuous_variables_;

};

}
}

#endif // VARIABLE_PROPERTY_CONTAINER_H
