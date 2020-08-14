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

#include "var_prop_container.h"
#include <QStringList>
#include <iostream>

namespace Model {
namespace Properties {

//VarPropContainer::VarPropContainer() {
//  binary_variables_ = new QHash<QUuid, BinaryProperty *>();
//  discrete_variables_ = new QHash<QUuid, DiscreteProperty *>();
//  continuous_variables_ = new QHash<QUuid, ContinuousProperty *>();
//}

VarPropContainer::VarPropContainer() {
  binary_variables_ = new QMap<QUuid, BinaryProperty *>();
  discrete_variables_ = new QMap<QUuid, DiscreteProperty *>();
  continuous_variables_ = new QMap<QUuid, ContinuousProperty *>();
}

void VarPropContainer::AddVariable(BinaryProperty *var) {
  var->SetVariable();
  binary_variables_->insert(var->id(), var);
}

void VarPropContainer::AddVariable(DiscreteProperty *var) {
  var->SetVariable();
  discrete_variables_->insert(var->id(), var);
}

void VarPropContainer::AddVariable(ContinuousProperty *var) {
  var->SetVariable();
  continuous_variables_->insert(var->id(), var);
}

BinaryProperty *VarPropContainer::GetBinaryVariable(QUuid id) const {
  string tm = "Binary variable not found.";
  if (!binary_variables_->contains(id))
    throw VariableIdDoesNotExistException(tm);
  return binary_variables_->value(id);
}

DiscreteProperty *VarPropContainer::GetDiscreteVariable(QUuid id) const {

  string tm = "Discrete variable not found.";
  if (!discrete_variables_->contains(id))
    throw VariableIdDoesNotExistException(tm);
  return discrete_variables_->value(id);
}

ContinuousProperty *VarPropContainer::GetContinousVariable(QUuid id) const {

  string tm = "Continuous variable not found.";
  if (!continuous_variables_->contains(id))
    throw VariableIdDoesNotExistException(tm);
  return continuous_variables_->value(id);
}

void VarPropContainer::SetBinaryVariableValue(QUuid id, bool val) {
  string tm = "Binary variable not found.";
  if (!binary_variables_->contains(id))
    throw VariableIdDoesNotExistException(tm);
  else binary_variables_->value(id)->setValue(val);
}

void VarPropContainer::SetDiscreteVariableValue(QUuid id, int val) {
  string tm = "Integer variable not found.";
  if (!discrete_variables_->contains(id))
    throw VariableIdDoesNotExistException(tm);
  else discrete_variables_->value(id)->setValue(val);
}

void VarPropContainer::SetContinousVariableValue(QUuid id, double val) {
  string tm = "Continuous variable not found.";
  if (!continuous_variables_->contains(id))
    throw VariableIdDoesNotExistException(tm);
  else continuous_variables_->value(id)->setValue(val);
}

//QHash<QUuid, bool> VarPropContainer::GetBinaryVariableValues() const {
//  QHash<QUuid, bool> binary_values = QHash<QUuid, bool>();
//  for (QUuid key : binary_variables_->keys())
//    binary_values[key] = binary_variables_->value(key)->value();
//  return binary_values;
//}

QMap<QUuid, bool> VarPropContainer::GetBinaryVariableValues() const {
  QMap<QUuid, bool> binary_values = QMap<QUuid, bool>();
  for (QUuid key : binary_variables_->keys())
    binary_values[key] = binary_variables_->value(key)->value();
  return binary_values;
}

//QHash<QUuid, int> VarPropContainer::GetDiscreteVariableValues() const {
//  QHash<QUuid, int> discrete_values = QHash<QUuid, int>();
//  for (QUuid key : discrete_variables_->keys()) {
//    discrete_values[key] = discrete_variables_->value(key)->value();
//  }
//  return discrete_values;
//}

QMap<QUuid, int> VarPropContainer::GetDiscreteVariableValues() const {
  QMap<QUuid, int> discrete_values = QMap<QUuid, int>();
  for (QUuid key : discrete_variables_->keys()) {
    discrete_values[key] = discrete_variables_->value(key)->value();
  }
  return discrete_values;
}

//QHash<QUuid, double> VarPropContainer::GetContinuousVariableValues() const {
//  QHash<QUuid, double> continous_values = QHash<QUuid, double>();
//  for (QUuid key : continuous_variables_->keys()) {
//    continous_values[key] = continuous_variables_->value(key)->value();
//  }
//  return continous_values;
//}

QMap<QUuid, double> VarPropContainer::GetContinuousVariableValues() const {
  QMap<QUuid, double> continous_values = QMap<QUuid, double>();
  for (QUuid key : continuous_variables_->keys()) {
    continous_values[key] = continuous_variables_->value(key)->value();
  }
  return continous_values;
}

void VarPropContainer::CheckVariableNameUniqueness() {
  QList<QString> names = QList<QString>();
  string tm = "Encountered non-unique variable name: ";

  for (auto var : discrete_variables_->values()) {
    if (var->name().size() > 0) {
      if (names.contains(var->name()))
        throw std::runtime_error(tm + var->name().toStdString());
      else
        names.append(var->name());
    }
  }

  for (auto var : continuous_variables_->values()) {
    if (var->name().size() > 0) {
      if (names.contains(var->name()))
        throw std::runtime_error(tm + var->name().toStdString());
      else
        names.append(var->name());
    }
  }

  for (auto var : binary_variables_->values()) {
    if (var->name().size() > 0) {
      if (names.contains(var->name()))
        throw std::runtime_error(tm + var->name().toStdString());
      else
        names.append(var->name());
    }
  }
}

QList<ContinuousProperty *> VarPropContainer::GetWellSplineVariables(const QString well_name) const {
  QList<ContinuousProperty *> spline_vars;
  for (auto var : continuous_variables_->values()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::SplinePoint &&
      QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      spline_vars.append(var);
    }
  }
  return spline_vars;
}

QList<ContinuousProperty *> VarPropContainer::GetPolarSplineVariables(const QString well_name) const {
  QList<ContinuousProperty *> polar_vars;
  for (auto var : continuous_variables_->values()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::PolarSpline &&
      well_name == var->propertyInfo().parent_well_name) {
      polar_vars.append(var);
    }
  }
  return polar_vars;
}

QList<ContinuousProperty *> VarPropContainer::GetPseudoContVertVariables(const QString well_name) const {
  QList<ContinuousProperty *> pcv_vars;

  for (auto var : continuous_variables_->values()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::PseudoContVert &&
      QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      pcv_vars.append(var);
    }
  }

  if (pcv_vars.size() != 2) {
    std::cerr << "Error: Incorrect number (" << pcv_vars.size()
              << ") found for well " << well_name.toStdString() << std::endl;
    throw std::runtime_error("Incorrect number of pseudocontinuous variables found.");
  }

  if (pcv_vars[0]->propertyInfo().coord != Property::Coordinate::x)
    pcv_vars.swap(0, 1);
  return pcv_vars;
}


QList<ContinuousProperty *> VarPropContainer::GetWellControlVariables() const {
  QList<ContinuousProperty *> well_controls;
  for (auto var : continuous_variables_->values()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::BHP||
      var->propertyInfo().prop_type == Property::PropertyType::Rate) {
      well_controls.append(var);
    }
  }
  std::sort(well_controls.begin(), well_controls.end(), [](ContinuousProperty *p1, ContinuousProperty *p2) {
    int time_1 = p1->propertyInfo().index;
    int time_2 = p2->propertyInfo().index;
    return time_1 < time_2;
  });
  return well_controls;
}

QList<ContinuousProperty *> VarPropContainer::GetWellControlVariables(const QString well_name) const {
  QList<ContinuousProperty *> well_controls;
  for (auto var : GetWellControlVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      well_controls.append(var);
    }
  }
  return well_controls;
}

QList<ContinuousProperty *> VarPropContainer::GetWellBHPVariables() const {
  QList<ContinuousProperty *> bhp_variables;
  for (auto var : GetWellControlVariables()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::BHP)
      bhp_variables.append(var);
  }
  return bhp_variables;
}

QList<ContinuousProperty *> VarPropContainer::GetWellRateVariables() const {
  QList<ContinuousProperty *> rate_variables;
  for (auto var : GetWellControlVariables()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::Rate)
      rate_variables.append(var);
  }
  return rate_variables;
}

QList<ContinuousProperty *> VarPropContainer::GetWellBHPVariables(QString well_name) const {
  QList<ContinuousProperty *> well_bhp_variables;
  for (auto var : GetWellBHPVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      well_bhp_variables.append(var);
    }
  }
  return well_bhp_variables;
}

QList<ContinuousProperty *> VarPropContainer::GetWellRateVariables(QString well_name) const {
  QList<ContinuousProperty *> well_rate_variables;
  for (auto var : GetWellRateVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      well_rate_variables.append(var);
    }
  }
  return well_rate_variables;
}

QList<DiscreteProperty *> VarPropContainer::GetWellBlockVariables() const {
  QList<DiscreteProperty *> wb_vars;
  for (auto var : discrete_variables_->values()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::WellBlock)
      wb_vars.append(var);
  }
  return wb_vars;
}

QList<DiscreteProperty *> VarPropContainer::GetWellBlockVariables(const QString well_name) const {
  QList<DiscreteProperty *> wb_vars;
  for (auto var : GetWellBlockVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0)
      wb_vars.append(var);
  }
  return wb_vars;
}

QList<ContinuousProperty *> VarPropContainer::GetTransmissibilityVariables() const {
  QList<ContinuousProperty *> trans_vars;
  for (auto var : continuous_variables_->values()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::Transmissibility)
      trans_vars.append(var);
  }
  return trans_vars;
}

QList<ContinuousProperty *> VarPropContainer::GetTransmissibilityVariables(const QString well_name) const {
  QList<ContinuousProperty *> trans_vars;
  for (auto var : GetTransmissibilityVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0)
      trans_vars.append(var);
  }
  return trans_vars;
}

BinaryProperty *VarPropContainer::GetBinaryVariable(QString name) const {
  for (auto var : binary_variables_->values()) {
    if (QString::compare(var->name(), name) == 0)
      return var;
  }
  throw std::runtime_error("Unable to find binary variable with name " + name.toStdString());
}

DiscreteProperty *VarPropContainer::GetDiscreteVariable(QString name) const {
  string tm = "Unable to find discrete variable with name ";
  for (auto var : discrete_variables_->values()) {
    if (QString::compare(var->name(), name) == 0)
      return var;
  }
  throw std::runtime_error(tm + name.toStdString());
}

ContinuousProperty *VarPropContainer::GetContinousVariable(QString name) const {
  string tm = "Unable to find continous variable with name ";
  for (auto var : continuous_variables_->values()) {
    if (QString::compare(var->name(), name) == 0)
      return var;
  }
  throw std::runtime_error(tm + name.toStdString());
}

QList<ContinuousProperty *> VarPropContainer::GetContinuousProperties(QList<QUuid> ids) const {
  QList<ContinuousProperty *> props;
  for (int i = 0; i < ids.size(); ++i) {
    props.append(continuous_variables_->value(ids[i]));
  }
  return props;
}

}
}
