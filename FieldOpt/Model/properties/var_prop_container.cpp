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

#include "var_prop_container.h"
#include <QStringList>
#include <iostream>

namespace Model {
namespace Properties {

// MAKE VARIABLE CONTAINER ---------------------------------
//VarPropContainer::VarPropContainer() {
//  binary_variables_ = new QHash<QUuid, BinaryProperty *>();
//  discrete_variables_ = new QHash<QUuid, DiscreteProperty *>();
//  continuous_variables_ = new QHash<QUuid, ContinuousProperty *>();
//}

//VarPropContainer::VarPropContainer() {
//  binary_variables_ = new QMap<QUuid, BinaryProperty *>();
//  discrete_variables_ = new QMap<QUuid, DiscreteProperty *>();
//  continuous_variables_ = new QMap<QUuid, ContinuousProperty *>();
//}

VarPropContainer::VarPropContainer(Settings::VerbParams vp) {
  binary_variables_ = new QList<QPair<QUuid, BinaryProperty *>>();
  discrete_variables_ = new QList<QPair<QUuid, DiscreteProperty *>>();
  continuous_variables_ = new QList<QPair<QUuid, ContinuousProperty *>>();
  vp_ = vp;
}

// ADD VARIABLE --------------------------------------------
void VarPropContainer::AddVariable(BinaryProperty *var) {
  var->SetVariable();
  // binary_variables_->insert(var->id(), var);
  binary_variables_->append(qMakePair(var->id(), var));
}

void VarPropContainer::AddVariable(DiscreteProperty *var) {
  var->SetVariable();
  // discrete_variables_->insert(var->id(), var);
  discrete_variables_->append(qMakePair(var->id(), var));
}

void VarPropContainer::AddVariable(ContinuousProperty *var) {
  var->SetVariable();
  // continuous_variables_->insert(var->id(), var);
  continuous_variables_->append(qMakePair(var->id(), var));
  if (vp_.vOPT > 3) {
    ext_info("@[AddVariable]: " + var->ToString(), md_, cl_);
  }
}

// GET BINARY VAR ------------------------------------------
BinaryProperty* VarPropContainer::GetBinaryVariable(QUuid id) const {
  for(int ii=0; ii < binary_variables_->size(); ii++) {
    if (binary_variables_->at(ii).first==id) {
      return binary_variables_->at(ii).second;
    }
  }
  string tm = "Binary variable not found.";
  throw VariableIdDoesNotExistException(tm);
}

// GET DISCRETE VAR ----------------------------------------
DiscreteProperty *VarPropContainer::GetDiscreteVariable(QUuid id) const {
  for(int ii=0; ii < discrete_variables_->size(); ii++) {
    if (discrete_variables_->at(ii).first==id) {
      return discrete_variables_->at(ii).second;
    }
  }
  string tm = "Discrete variable not found.";
  throw VariableIdDoesNotExistException(tm);
}

// GET CONTINUOUS VAR --------------------------------------
ContinuousProperty *VarPropContainer::GetContinuousVariable(QUuid id) const {
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    if (continuous_variables_->at(ii).first==id) {
      return continuous_variables_->at(ii).second;
    }
  }
  string tm = "Continuous variable not found.";
  throw VariableIdDoesNotExistException(tm);
}

// SET BINARY VALUE ----------------------------------------
void VarPropContainer::SetBinaryVariableValue(QUuid id, bool val) {
  for(int ii=0; ii < binary_variables_->size(); ii++) {
    if (binary_variables_->at(ii).first==id) {
      binary_variables_->at(ii).second->setValue(val);
      return;
    }
  }
  string tm = "Binary variable not found.";
  throw VariableIdDoesNotExistException(tm);
}

// SET DISCRETE VALUE --------------------------------------
void VarPropContainer::SetDiscreteVariableValue(QUuid id, int val) {
  for(int ii=0; ii < discrete_variables_->size(); ii++) {
    if (discrete_variables_->at(ii).first==id) {
      discrete_variables_->at(ii).second->setValue(val);
      return;
    }
  }
  string tm = "Integer variable not found.";
  throw VariableIdDoesNotExistException(tm);
}

// SET CONTINUOUS VALUE ------------------------------------
void VarPropContainer::SetContinuousVariableValue(QUuid id, double val) {
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    if (continuous_variables_->at(ii).first==id) {
      continuous_variables_->at(ii).second->setValueSc(val);
      // continuous_variables_->at(ii).second->setValue(val);
      return;
    }
  }
  string tm = "Continuous variable not found.";
  throw VariableIdDoesNotExistException(tm);
}

// BINARY VALUES -------------------------------------------
//QHash<QUuid, bool> VarPropContainer::GetBinVarValues() const {
//  QHash<QUuid, bool> binary_values = QHash<QUuid, bool>();
//  for (QUuid key : binary_variables_->keys())
//    binary_values[key] = binary_variables_->value(key)->value();
//  return binary_values;
//}

// QMap<QUuid, bool> VarPropContainer::GetBinVarValues() const {
//   QMap<QUuid, bool> binary_values = QMap<QUuid, bool>();
//   for (QUuid key : binary_variables_->keys())
//     binary_values[key] = binary_variables_->value(key)->value();
//   return binary_values;
// }

QList<QPair<QUuid, bool>> VarPropContainer::GetBinVarValues() const {
  QList<QPair<QUuid, bool>> bin_vals = QList<QPair<QUuid, bool>>();
  for(int ii=0; ii < binary_variables_->size(); ii++) {
    bin_vals.append(qMakePair(binary_variables_->at(ii).first,
                              binary_variables_->at(ii).second->value()));
  }
  return bin_vals;
}

// DISCRETEVALUES ---------------------------------------
//QHash<QUuid, int> VarPropContainer::GetDiscVarValues() const {
//  QHash<QUuid, int> discrete_values = QHash<QUuid, int>();
//  for (QUuid key : discrete_variables_->keys()) {
//    discrete_values[key] = discrete_variables_->value(key)->value();
//  }
//  return discrete_values;
//}

//QMap<QUuid, int> VarPropContainer::GetDiscVarValues() const {
//  QMap<QUuid, int> discrete_values = QMap<QUuid, int>();
//  for (QUuid key : discrete_variables_->keys()) {
//    discrete_values[key] = discrete_variables_->value(key)->value();
//  }
//  return discrete_values;
//}

QList<QPair<QUuid, int>> VarPropContainer::GetDiscVarValues() const {
  QList<QPair<QUuid, int>> int_vals = QList<QPair<QUuid, int>>();
  for(int ii=0; ii < discrete_variables_->size(); ii++) {
    int_vals.append(qMakePair(discrete_variables_->at(ii).first,
                              discrete_variables_->at(ii).second->value()));
  }
  return int_vals;
}

// CONTINUOUS VALUES ---------------------------------------
//QHash<QUuid, double> VarPropContainer::GetContVarValues() const {
//  QHash<QUuid, double> continous_values = QHash<QUuid, double>();
//  for (QUuid key : continuous_variables_->keys()) {
//    continous_values[key] = continuous_variables_->value(key)->value();
//  }
//  return continous_values;
//}

//QMap<QUuid, double> VarPropContainer::GetContVarValues() const {
//  QMap<QUuid, double> continous_values = QMap<QUuid, double>();
//  for (QUuid key : continuous_variables_->keys()) {
//    continous_values[key] = continuous_variables_->value(key)->value();
//  }
//  return continous_values;
//}

QList<QPair<QUuid, double>> VarPropContainer::GetContVarValues() const {
  QList<QPair<QUuid, double>> cont_vals = QList<QPair<QUuid, double>>();
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    cont_vals.append(qMakePair(continuous_variables_->at(ii).first,
                               continuous_variables_->at(ii).second->valueSc()));
    // cont_vals.append(qMakePair(continuous_variables_->at(ii).first,
    //                            continuous_variables_->at(ii).second->value()));
  }
  return cont_vals;
}

// CHECK VAR.NAME UNIQUENESS--------------------------------
void VarPropContainer::CheckVariableNameUniqueness() {
  QList<QString> names = QList<QString>();
  string tm = "Encountered non-unique variable name: ";

  for(int ii=0; ii < discrete_variables_->size(); ii++) {
    auto var = discrete_variables_->at(ii).second;
    if (var->name().size() > 0) {
      if (names.contains(var->name()))
        throw std::runtime_error(tm + var->name().toStdString());
      else
        names.append(var->name());
    }
  }

  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    auto var = continuous_variables_->at(ii).second;
    if (var->name().size() > 0) {
      if (names.contains(var->name()))
        throw std::runtime_error(tm + var->name().toStdString());
      else
        names.append(var->name());
    }
  }

  for(int ii=0; ii < binary_variables_->size(); ii++) {
    auto var = binary_variables_->at(ii).second;
    if (var->name().size() > 0) {
      if (names.contains(var->name()))
        throw std::runtime_error(tm + var->name().toStdString());
      else
        names.append(var->name());
    }
  }
}

// SPLINEVARS ----------------------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetWSplineVars(QString well_name) const {
  QList<ContinuousProperty *> spline_vars;
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    auto var = continuous_variables_->at(ii).second;
    if (var->propertyInfo().prop_type == Property::PropertyType::SplinePoint &&
        QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      spline_vars.append(var);
    }
  }
  return spline_vars;
}

// POLARSPLINE ---------------------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetPolarSplineVariables(const QString well_name) const {
  QList<ContinuousProperty *> polar_vars;
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    auto var = continuous_variables_->at(ii).second;
    if (var->propertyInfo().prop_type == Property::PropertyType::PolarSpline &&
        well_name == var->propertyInfo().parent_well_name) {
      polar_vars.append(var);
    }
  }
  return polar_vars;
}

// PSEUDOCONT ----------------------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetPseudoContVertVariables(const QString well_name) const {
  QList<ContinuousProperty *> pcv_vars;
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    auto var = continuous_variables_->at(ii).second;
    if (var->propertyInfo().prop_type == Property::PropertyType::PseudoContVert &&
        QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      pcv_vars.append(var);
    }
  }

  stringstream es;
  if (pcv_vars.size() != 2) {
    es << "Error: Incorrect number (" << pcv_vars.size();
    es << ") found for well " << well_name.toStdString() << endl;
    std::cerr << es.str();
    string em = "Incorrect number of pseudocontinuous variables found.";
    throw std::runtime_error(em);
  }

  if (pcv_vars[0]->propertyInfo().coord != Property::Coordinate::x)
    pcv_vars.swap(0, 1);
  return pcv_vars;
}

// CONTROL -------------------------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetWellControlVariables() const {
  QList<ContinuousProperty *> well_controls;
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    auto var = continuous_variables_->at(ii).second;
    if (var->propertyInfo().prop_type == Property::PropertyType::BHP
        || var->propertyInfo().prop_type == Property::PropertyType::Rate
        || var->propertyInfo().prop_type == Property::PropertyType::ICD) {
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

QList<ContinuousProperty *>
VarPropContainer::GetWellControlVariables(const QString well_name) const {
  QList<ContinuousProperty *> well_controls;
  for (auto var : GetWellControlVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      well_controls.append(var);
    }
  }
  return well_controls;
}

// BHP -----------------------------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetWellBHPVariables() const {
  QList<ContinuousProperty *> bhp_variables;
  for (auto var : GetWellControlVariables()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::BHP)
      bhp_variables.append(var);
  }
  return bhp_variables;
}

QList<ContinuousProperty *>
VarPropContainer::GetWellBHPVars(QString const &well_name) const {
  QList<ContinuousProperty *> well_bhp_variables;
  for (auto var : GetWellBHPVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      well_bhp_variables.append(var);
    }
  }
  return well_bhp_variables;
}

// ICD -----------------------------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetWellICDVariables() const {
  QList<ContinuousProperty *> icd_variables;
  for (auto var : GetWellControlVariables()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::ICD)
      icd_variables.append(var);
  }
  return icd_variables;
}

QList<ContinuousProperty *>
VarPropContainer::GetWellICDVars(QString const &well_name) const {
  QList<ContinuousProperty *> well_icd_variables;
  for (auto var : GetWellICDVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      well_icd_variables.append(var);
    }
  }
  return well_icd_variables;
}

// RATE ----------------------------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetWellRateVariables() const {
  QList<ContinuousProperty *> rate_variables;
  for (auto var : GetWellControlVariables()) {
    if (var->propertyInfo().prop_type == Property::PropertyType::Rate)
      rate_variables.append(var);
  }
  return rate_variables;
}

QList<ContinuousProperty *>
VarPropContainer::GetWellRateVars(QString well_name) const {
  QList<ContinuousProperty *> well_rate_variables;
  for (auto var : GetWellRateVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0) {
      well_rate_variables.append(var);
    }
  }
  return well_rate_variables;
}

// WELL BLOCKS ---------------------------------------------
QList<DiscreteProperty *>
VarPropContainer::GetWellBlockVariables() const {
  QList<DiscreteProperty *> wb_vars;
  for(int ii=0; ii < discrete_variables_->size(); ii++) {
    auto var = discrete_variables_->at(ii).second;
    if (var->propertyInfo().prop_type == Property::PropertyType::WellBlock)
      wb_vars.append(var);
  }
  return wb_vars;
}

QList<DiscreteProperty *>
VarPropContainer::GetWellBlockVariables(const QString well_name) const {
  QList<DiscreteProperty *> wb_vars;
  for (auto var : GetWellBlockVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0)
      wb_vars.append(var);
  }
  return wb_vars;
}

// TRANSMISSIBILITY VARS -----------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetTransmissibilityVariables() const {
  QList<ContinuousProperty *> trans_vars;
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    auto var = continuous_variables_->at(ii).second;
    if (var->propertyInfo().prop_type == Property::PropertyType::Transmissibility)
      trans_vars.append(var);
  }
  return trans_vars;
}

QList<ContinuousProperty *>
VarPropContainer::GetTransmissibilityVariables(const QString well_name) const {
  QList<ContinuousProperty *> trans_vars;
  for (auto var : GetTransmissibilityVariables()) {
    if (QString::compare(well_name, var->propertyInfo().parent_well_name) == 0)
      trans_vars.append(var);
  }
  return trans_vars;
}

// GET BINARY VAR ------------------------------------------
BinaryProperty*
VarPropContainer::GetBinaryVariable(QString name) const {
  for(int ii=0; ii < binary_variables_->size(); ii++) {
    auto var = binary_variables_->at(ii).second;
    if (QString::compare(var->name(), name) == 0)
      return var;
  }
  throw std::runtime_error("Unable to find binary variable with name " + name.toStdString());
}

// GET DISCRETE VAR ----------------------------------------
DiscreteProperty*
VarPropContainer::GetDiscreteVariable(QString name) const {
  string tm = "Unable to find discrete variable with name ";
  for(int ii=0; ii < discrete_variables_->size(); ii++) {
    auto var = discrete_variables_->at(ii).second;
    if (QString::compare(var->name(), name) == 0)
      return var;
  }
  throw std::runtime_error(tm + name.toStdString());
}

// GET CONTINUOUS VAR --------------------------------------
ContinuousProperty*
VarPropContainer::GetContinousVariable(QString name) const {
  string tm = "Unable to find continuous variable with name ";
  for(int ii=0; ii < continuous_variables_->size(); ii++) {
    auto var = continuous_variables_->at(ii).second;
    if (QString::compare(var->name(), name) == 0)
      return var;
  }
  throw std::runtime_error(tm + name.toStdString());
}

// GET CONTINUOUS PROPS ------------------------------------
QList<ContinuousProperty *>
VarPropContainer::GetContinuousProperties(QList<QUuid> ids) const {
  QList<ContinuousProperty *> props;
  for (int ii = 0; ii < ids.size(); ++ii) {
    for(int jj=0; jj < continuous_variables_->size(); jj++) {
      if (continuous_variables_->at(jj).first == ids[ii]) {
        props.append(continuous_variables_->at(ii).second);
      }
    }
  }
  return props;
}

// SCALE VARIABLES -----------------------------------------
void VarPropContainer::ScaleVariables(){
  stringstream ss;
  ss << "@[VarPropContainer::ScaleVariables]|";

  // constraint_handler constructor inserts the bounds
  // of each variable in the variable container

  double D = 2.0, c = 0.0;
  // x = y * bnds_b_mns_a_/2 + bnds_a_pls_b_/2
  // y = 2*x/bnds_b_mns_a_ - bnds_a_pls_b_/bnds_b_mns_a_

  for(int ii=0; ii < continuous_variables_->size(); ii++) {

    // set up scaling multiplier (D) and constant (c) terms
    auto bounds = continuous_variables_->value(ii).second->bounds();

    D = bounds(1) - bounds(0); // bnds_b_mns_a_
    c = bounds(0) + bounds(1); // bnds_a_pls_b_

    continuous_variables_->value(ii).second->setScalingCoeffs(D, c);
    continuous_variables_->value(ii).second->scaleValue();

    if (vp_.vMOD >= 3) {
      ss << continuous_variables_->value(ii).second->ToString();
    }
  }

  ext_info(ss.str(), md_, cl_);

  // ext_warn(wm_, md_, cl_);
  // em_ = "";
  // throw runtime_error(em_);
}

}
}
