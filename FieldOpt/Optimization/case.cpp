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

#include <boost/serialization/map.hpp>
#include <QtCore/QUuid>
#include <Utilities/time.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <Utilities/math.hpp>
#include "case.h"
#include <Utilities/printer.hpp>

namespace Optimization {

Case::Case() {
  id_ = QUuid::createUuid();

  // binary_variables_ = QHash<QUuid, bool>();
  // integer_variables_ = QHash<QUuid, int>();
  // real_variables_ = QHash<QUuid, double>();

  //binary_variables_ = QMap<QUuid, bool>();
  //integer_variables_ = QMap<QUuid, int>();
  //real_variables_ = QMap<QUuid, double>();

  binary_variables_ = QList<QPair<QUuid, bool>>();
  integer_variables_ = QList<QPair<QUuid, int>>();
  real_variables_ = QList<QPair<QUuid, double>>();
  rvar_grads_ = QList<QPair<QUuid, double>>();

  objf_val_ = std::numeric_limits<double>::max();
  sim_time_sec_ = 0;
  wic_time_sec_ = 0;
  ensemble_realization_ = "";
  ensemble_ofvs_ = QHash<QString, double>();
}

Case::Case(const QList<QPair<QUuid, bool>> &binary_variables,
           const QList<QPair<QUuid, int>> &integer_variables,
           const QList<QPair<QUuid, double>> &real_variables) {

  id_ = QUuid::createUuid();
  binary_variables_ = binary_variables;
  integer_variables_ = integer_variables;
  real_variables_ = real_variables;
  // rvar_grads_ = real_variables;

  for(int ii=0; ii < binary_variables_.size(); ii++) {
    binary_id_index_map_.append(binary_variables_.at(ii).first);
  }

  for(int ii=0; ii < integer_variables_.size(); ii++) {
    integer_id_index_map_.append(integer_variables_.at(ii).first);
  }

  for(int ii=0; ii < real_variables_.size(); ii++) {
    real_id_index_map_.append(real_variables_.at(ii).first);
    rvar_grads_.append(qMakePair(real_variables_.at(ii).first, 0.0));
  }

  objf_val_ = std::numeric_limits<double>::max();
  sim_time_sec_ = 0;
  wic_time_sec_ = 0;
  ensemble_realization_ = "";
  ensemble_ofvs_ = QHash<QString, double>();
}

Case::Case(const Case *c) {

  id_ = QUuid::createUuid();
  binary_variables_ = QList<QPair<QUuid, bool>> (c->binary_variables());
  integer_variables_ = QList<QPair<QUuid, int>> (c->integer_variables());
  real_variables_ = QList<QPair<QUuid, double>> (c->real_variables());
  rvar_grads_ = QList<QPair<QUuid, double>>(c->real_var_grads());

  binary_id_index_map_ = c->binary_id_index_map_;
  integer_id_index_map_ = c->integer_id_index_map_;
  real_id_index_map_ = c->real_id_index_map_;

  // integer_id_index_map_ = c->integer_variables_.keys();
  // for(int ii=0; ii < integer_variables_.size(); ii++) {
  //   integer_id_index_map_.append(integer_variables_.at(ii).first);
  // }

  objf_val_ = c->objf_val_;
  sim_time_sec_ = 0;
  wic_time_sec_ = 0;
  ensemble_realization_ = "";
  ensemble_ofvs_ = c->ensemble_ofvs_;
  vp_ = c->vp_;
}

void Case::CopyCaseVals(const Case *c) {
  binary_variables_old_ = binary_variables_;
  integer_variables_old_ = integer_variables_;
  real_variables_old_ = real_variables_;
  rvar_grads_old_ = rvar_grads_;

  binary_variables_ = QList<QPair<QUuid, bool>> (c->binary_variables());
  integer_variables_ = QList<QPair<QUuid, int>> (c->integer_variables());
  real_variables_ = QList<QPair<QUuid, double>> (c->real_variables());
  rvar_grads_ = QList<QPair<QUuid, double>>(c->real_var_grads());

  binary_id_index_map_ = c->binary_id_index_map_;
  integer_id_index_map_ = c->integer_id_index_map_;
  real_id_index_map_ = c->real_id_index_map_;

  objf_val_ = c->objf_val_;
  sim_time_sec_ = c->sim_time_sec_;
  wic_time_sec_ = c->wic_time_sec_;
}

bool Case::Equals(const Case *other, double tolerance) const {

  // Check if number of variables are equal
  if (this->binary_variables().size() != other->binary_variables().size()
    || this->integer_variables().size() != other->integer_variables().size()
    || this->real_variables().size() != other->real_variables().size())
    return false;

  // for (QUuid key : this->binary_variables().keys()) {
  //   if (std::abs(this->binary_variables()[key] - other->binary_variables()[key]) > tolerance)
  //     return false;
  // }

  for(int ii=0; ii < binary_variables_.size(); ii++) {
    auto lhs = this->binary_variables_.at(ii).second;
    auto rhs = other->binary_variables_.at(ii).second;
    if (std::abs(lhs - rhs) > tolerance) {
      return false;
    }
  }

  // for (QUuid key : this->integer_variables().keys()) {
  //   if (std::abs(this->integer_variables()[key] - other->integer_variables()[key]) > tolerance)
  //     return false;
  // }

  for(int ii=0; ii < integer_variables_.size(); ii++) {
    auto lhs = this->integer_variables_.at(ii).second;
    auto rhs = other->integer_variables_.at(ii).second;
    if (std::abs(lhs - rhs) > tolerance) {
      return false;
    }
  }

  // for (QUuid key : this->real_variables().keys()) {
  //   if (std::abs(this->real_variables()[key] - other->real_variables()[key]) > tolerance)
  //     return false;
  // }

  for(int ii=0; ii < real_variables_.size(); ii++) {
    auto lhs = this->real_variables_.at(ii).second;
    auto rhs = other->real_variables_.at(ii).second;
    if (std::abs(lhs - rhs) > tolerance) {
      return false;
    }
  }

  return true; // All variable values are equal if we reach this point.
}

double Case::objf_value() const {
  if (objf_val_ == std::numeric_limits<double>::max()) {
    string em = "The objective function value has not been set in this Case.";
    throw ObjectiveFunctionException(em);
  } else {
    return objf_val_;
  }
}

void Case::set_integer_variable_value(const QUuid id, const int val) {
  for(int ii=0; ii < integer_variables_.size(); ii++) {
    if (integer_variables_.at(ii).first==id) {
      auto pair_int = qMakePair(integer_variables_.at(ii).first, val);
      integer_variables_.replace(ii, pair_int);
      return;
    }
  }
  em_ = "Unable to set value of variable " + id.toString().toStdString();
  throw VariableException(em_);
}

void Case::set_binary_variable_value(const QUuid id, const bool val) {
  for(int ii=0; ii < binary_variables_.size(); ii++) {
    if (binary_variables_.at(ii).first==id) {
      auto pair_int = qMakePair(binary_variables_.at(ii).first, val);
      binary_variables_.replace(ii, pair_int);
      return;
    }
  }
  em_ = "Unable to set value of variable " + id.toString().toStdString();
  throw VariableException(em_);
}

void Case::set_real_variable_value(const QUuid id, const double val) {
  for(int ii=0; ii < real_variables_.size(); ii++) {
    if (real_variables_.at(ii).first==id) {
      auto pair_real = qMakePair(real_variables_.at(ii).first, val);
      real_variables_.replace(ii, pair_real);
      return;
    }
  }
  em_ = "Unable to set value of variable " + id.toString().toStdString();
  throw VariableException(em_);
}

void Case::set_real_variable_value(const QUuid id, const double val, const double grad) {
  for(int ii=0; ii < real_variables_.size(); ii++) {
    if (real_variables_.at(ii).first==id) {
      // add real value
      auto pair_real = qMakePair(real_variables_.at(ii).first, val);
      real_variables_.replace(ii, pair_real);
      // add corresponding gradient
      auto pair_grad = qMakePair(real_variables_.at(ii).first, grad);
      rvar_grads_.replace(ii, pair_grad);
      return;
    }
  }
  em_ = "Unable to set value of variable " + id.toString().toStdString();
  throw VariableException(em_);
}

QList<Case *> Case::Perturb(QUuid variabe_id, Case::SIGN sign, double magnitude) {
  QList<Case *> new_cases = QList<Case *>();
  int var_indx;

  // if (this->integer_variables().contains(variabe_id)) {
  if (this->integer_vars_contain(var_indx, variabe_id)) {
    if (sign == PLUS || sign == PLUSMINUS) {
      Case *new_case_p = new Case(this);

      // new_case_p->integer_variables_[variabe_id] += magnitude;
      auto pval = integer_variables_.at(var_indx).second + magnitude;
      auto pair = qMakePair(integer_variables_.at(var_indx).first, pval);
      new_case_p->integer_variables_.replace(var_indx, pair);

      new_case_p->objf_val_ = std::numeric_limits<double>::max();
      new_cases.append(new_case_p);
    }
    if (sign == MINUS || sign == PLUSMINUS) {
      Case *new_case_m = new Case(this);

      // new_case_m->integer_variables_[variabe_id] -= magnitude;
      auto pval = integer_variables_.at(var_indx).second - magnitude;
      auto pair = qMakePair(integer_variables_.at(var_indx).first, pval);
      new_case_m->integer_variables_.replace(var_indx, pair);

      new_case_m->objf_val_ = std::numeric_limits<double>::max();
      new_cases.append(new_case_m);
    }

    // } else if (real_variables_.contains(variabe_id)) {
  } else if (real_vars_contain(var_indx, variabe_id)) {
    if (sign == PLUS || sign == PLUSMINUS) {
      Case *new_case_p = new Case(this);

      // new_case_p->real_variables_[variabe_id] += magnitude;
      auto pval = real_variables_.at(var_indx).second + magnitude;
      auto pair = qMakePair(real_variables_.at(var_indx).first, pval);
      new_case_p->real_variables_.replace(var_indx, pair);

      new_case_p->objf_val_ = std::numeric_limits<double>::max();
      new_cases.append(new_case_p);
    }
    if (sign == MINUS || sign == PLUSMINUS) {
      Case *new_case_m = new Case(this);

      // new_case_m->real_variables_[variabe_id] -= magnitude;
      auto pval = real_variables_.at(var_indx).second - magnitude;
      auto pair = qMakePair(real_variables_.at(var_indx).first, pval);
      new_case_m->real_variables_.replace(var_indx, pair);

      new_case_m->objf_val_ = std::numeric_limits<double>::max();
      new_cases.append(new_case_m);
    }
  }
  return new_cases;
}

Eigen::VectorXd Case::GetVarVector(QList<QUuid> vars) {
  Eigen::VectorXd vec(vars.length());
  for (int i = 0; i < vars.length(); ++i) {
    vec[i] = get_real_variable_value(vars.at(i));
  }
  return vec;
}

Eigen::VectorXd Case::GetRealVarVector() {
  Eigen::VectorXd vec(real_id_index_map_.length());
  for (int i = 0; i < real_id_index_map_.length(); ++i) {
    // vec[i] = real_variables_.value(real_id_index_map_[i]);
    vec[i] = real_variables_.at(i).second;
  }
  return vec;
}

void Case::GetGradsOrderedByType(
  map<Property::PropertyType, VectorXd> &vt_map,
  VarPropContainer *vars) {
  // Pres for grad normalization by variable type
  // Organized grad vals by var type into map object
  auto varTypes = vars->GetVarTypes();
  stringstream ss; // dbg

  for (int ii = 0; ii < varTypes.size(); ii++) { // loop-by-var-type

    VectorXd vals = VectorXd::Zero(varTypes[ii].size());
    for (int jj = 0; jj < varTypes[ii].size(); jj++) {
      vals(jj) = get_list_value(varTypes[ii][jj]->id(), rvar_grads_scal_);
    }

    auto pt = varTypes[ii][0]->propertyInfo().prop_type;
    auto tv = pair<Property::PropertyType, VectorXd>(pt, vals);
    vt_map.insert(tv);

    // dbg
    if(vp_.vOPT > 1) {
      ss << varTypes[ii][0]->propertyInfo().getPropTypeName();
      ss << " -> var.norm() = " << vals.norm() << " | ";
      ext_info(ss.str(), md_, cl_);
    }
  }
}

vector<VectorXd> Case::GetRealVarGradx(VarPropContainer *vars) {
  VectorXd g = VectorXd(real_id_index_map_.length());
  vector<VectorXd> grads = vector<VectorXd>(2, g);
  rvar_grads_scal_ = QList<QPair<QUuid, double>>();

  // Scaling: Get grad value, then scale it w/i var container
  for (int i = 0; i < real_id_index_map_.length(); ++i) {
    auto vid = rvar_grads_.at(i).first;
    auto val = rvar_grads_.at(i).second;
    // cont var obj has bounds info to perform scaling
    vars->GetContinuousVariable(vid)->scaleGrad(val);
    // update grad construct/vector
    rvar_grads_scal_.append(QPair<QUuid, double>(vid, val));
    grads[0][i] = val;
  }

  // Normalization prep: Get grad vals (scaled) ordered by var type
  map<Property::PropertyType, VectorXd> vTypes_map;
  GetGradsOrderedByType(vTypes_map, vars);

  // Normalization: Normalize grad components according to variable types
  for (int i = 0; i < rvar_grads_scal_.length(); ++i) {
    auto vid = rvar_grads_scal_.at(i).first;
    auto val = rvar_grads_scal_.at(i).second;
    val /= vTypes_map[vars->GetContinuousVariable(vid)->propertyInfo().prop_type].norm();
    rvar_grads_norm_.append(QPair<QUuid, double>(vid, val));
    grads[1][i] = val;
  }

  return grads;
}

void Case::SetRealVarValues(Eigen::VectorXd vec) {
  for (int i = 0; i < vec.size(); ++i) {
    set_real_variable_value(real_id_index_map_[i], vec[i]);
  }
}

void Case::SetRealVarValues(Eigen::VectorXd vec, Eigen::VectorXd grad) {
  for (int i = 0; i < vec.size(); ++i) {
    set_real_variable_value(real_id_index_map_[i], vec[i], grad[i]);
  }
}

Eigen::VectorXi Case::GetIntegerVarVector() {
  Eigen::VectorXi vec(integer_id_index_map_.length());
  for (int i = 0; i < integer_id_index_map_.length(); ++i) {
    // vec[i] = integer_variables_.value(integer_id_index_map_[i]);
    vec[i] = integer_variables_.at(i).second;
  }
  return vec;
}

void Case::SetIntegerVarValues(Eigen::VectorXi vec) {
  for (int i = 0; i < vec.size(); ++i) {
    set_integer_variable_value(integer_id_index_map_[i], vec[i]);
  }
}

void Case::set_origin_data(Case *parent,
                           int direction_index,
                           double step_length) {
  parent_ = parent;
  direction_index_ = direction_index;
  step_length_ = step_length;
}

Loggable::LogTarget Case::GetLogTarget() {
  return Loggable::LogTarget::LOG_CASE;
}

map <string, string> Case::GetState() {
  map<string, string> statemap;
  switch (state.eval) {
    case CaseState::EvalStatus::E_FAILED: statemap["EvalSt"] = "FAIL"; break;
    case CaseState::EvalStatus::E_TIMEOUT: statemap["EvalSt"] = "TMOT"; break;
    case CaseState::EvalStatus::E_PENDING: statemap["EvalSt"] = "PEND"; break;
    case CaseState::EvalStatus::E_CURRENT: statemap["EvalSt"] = "CRNT"; break;
    case CaseState::EvalStatus::E_DONE: statemap["EvalSt"] = "OKAY"; break;
    case CaseState::EvalStatus::E_BOOKKEEPED: statemap["EvalSt"] = "BKPD"; break;
  }
  switch (state.cons) {
    case CaseState::ConsStatus::C_PROJ_FAILED: statemap["ConsSt"] = "PNFL"; break;
    case CaseState::ConsStatus::C_INFEASIBLE: statemap["ConsSt"] = "INFS"; break;
    case CaseState::ConsStatus::C_PENDING: statemap["ConsSt"] = "PEND"; break;
    case CaseState::ConsStatus::C_FEASIBLE: statemap["ConsSt"] = "OKAY"; break;
    case CaseState::ConsStatus::C_PROJECTED: statemap["ConsSt"] = "PROJ"; break;
    case CaseState::ConsStatus::C_PENALIZED: statemap["ConsSt"] = "PNZD"; break;
  }
  switch (state.err_msg) {
    case CaseState::ErrorMessage::ERR_SIM: statemap["ErrMsg"] = "SIML"; break;
    case CaseState::ErrorMessage::ERR_WIC: statemap["ErrMsg"] = "WLIC"; break;
    case CaseState::ErrorMessage::ERR_CONS: statemap["ErrMsg"] = "CONS"; break;
    case CaseState::ErrorMessage::ERR_UNKNOWN: statemap["ErrMsg"] = "UNWN"; break;
    case CaseState::ErrorMessage::ERR_OK: statemap["ErrMsg"] = "OKAY"; break;
  }
  return statemap;
}

QUuid Case::GetId() {
  return id();
}

map <string, vector<double>> Case::GetValues() {
  map<string, vector<double>> valmap;
  valmap["OFnVal"] = vector<double>{objf_val_};
  valmap["SimDur"] = vector<double>{sim_time_sec_};
  valmap["WicDur"] = vector<double>{wic_time_sec_};
  if (ensemble_ofvs_.size() > 1) {
    valmap["OFvSTD"] = vector<double>{GetEnsembleExpectedOfv().second};
  }
  return valmap;
}

string Case::StringRepresentation(Model::Properties::VarPropContainer *varcont) {
  stringstream str;
  str << "=========================================================|";
  str << "> fobj: " << num2str(objf_val_, 3, 1);
  str << " cs: " << id_stdstr() << " |";
  str << "---------------------------------------------------------|";

  if (!real_variables_.empty()) {
    str << "Continuous variable values:                             |";
    for (int ii=0; ii < real_variables_.size(); ii++) {
      string vn = varcont->GetContinuousVariables()->at(ii).second->name().toStdString();
      double vv = real_variables_.at(ii).second;
      str << ">  " << vn << ": " << num2str(vv,3,0,20)  << " |";

    }
  }

  if (!integer_variables_.empty()) {
    str << "| Discrete variable values:                             |" << endl;
    for (int ii=0; ii < real_variables_.size(); ii++) {
      string vn = varcont->GetDiscreteVariables()->at(ii).second->name().toStdString();
      double vv = integer_variables_.at(ii).second;
      str << "| > " << vn << ": " << std::setw (51 - vn.size()) << vv << " |" << endl;
    }
  }

  if (!binary_variables_.empty()) {
    str << "| Discrete variable values:                             |" << endl;
    for (int ii=0; ii < real_variables_.size(); ii++) {
      string vn = varcont->GetBinaryVariables()->at(ii).second->name().toStdString();
      double vv = binary_variables_.at(ii).second;
      str << "| > " << vn << ": " << std::setw (51 - vn.size()) << vv << " |" << endl;
    }
  }
  str << "=========================================================|";
  return str.str();

//    tringstream entry;
//    entry << setw(cas_log_col_widths_["TimeSt"]) << timestamp_string() << " ,";
//    entry << setw(cas_log_col_widths_["EvalSt"]) << obj->GetState()["EvalSt"] << " ,";
//    entry << setw(cas_log_col_widths_["ConsSt"]) << obj->GetState()["ConsSt"] << " ,";
//    entry << setw(cas_log_col_widths_["ErrMsg"]) << obj->GetState()["ErrMsg"] << " ,";
//    entry << setw(cas_log_col_widths_["SimDur"]) << timespan_string(obj->GetValues()["SimDur"][0]) << " , ";
//    entry << setw(cas_log_col_widths_["WicDur"]) << timespan_string(obj->GetValues()["WicDur"][0]) << " , ";
//    entry.precision(6);
//    entry << setw(cas_log_col_widths_["OFnVal"]) << scientific << obj->GetValues()["OFnVal"][0] << " ,";
//    entry << setw(cas_log_col_widths_["CaseId"]) << obj->GetId().toString().toStdString();
//    string str = entry.str();
//    Utilities::FileHandling::WriteLineToFile(QString::fromStdString(str), cas_log_path_);
}

void Case::SetRealizationOfv(const QString &alias, const double &ofv) {
  ensemble_ofvs_[alias] = ofv;
}

double Case::GetRealizationOfv(const QString &alias) {
  if (HasRealizationOfv(alias)) {
    return ensemble_ofvs_[alias];
  } else {
    wm_ = "In Case - Attempting to get unrecorded OFV for alias ";
    wm_ += alias.toStdString() + ". Returning sentinel value (-1.0).";
    ext_warn(wm_, md_, cl_);
    return -1.0;
  }
}

bool Case::HasRealizationOfv(const QString &alias) {
  return ensemble_ofvs_.count(alias) > 0;
}

double Case::GetEnsembleAverageOfv() const {
  if (ensemble_ofvs_.size() == 1) {
    string wm = "Only one realization case was successfully ";
    wm += "evaluated. You should consider tweaking well/reservoir ";
    wm += "parameters or increasing the number of realizations ";
    wm += "considered for each case.";
    ext_warn(wm, md_, cl_);
  }
  double sum = 0;
  for (double value : ensemble_ofvs_.values()) {
    sum += value;
  }
  return sum / ensemble_ofvs_.size();
}

QPair<double, double> Case::GetEnsembleExpectedOfv() const {
  if (ensemble_ofvs_.size() == 1) {
    string wm = "Only one realization case was successfully ";
    wm += "evaluated. You should consider tweaking well/reservoir ";
    wm += "parameters or increasing the number of realizations ";
    wm += "considered for each case.";
    ext_warn(wm, md_, cl_);
  }
  vector<double> list_of_ofvs;
  for (double value : ensemble_ofvs_.values()){
    list_of_ofvs.push_back(value);
  }
  auto ofvs_avg = calc_average(list_of_ofvs);
  auto ofvs_dev = calc_standard_deviation(list_of_ofvs);
  auto pair = QPair<double, double>(ofvs_avg, ofvs_dev);
  return pair;
}

void Case::set_objf_value(double fval) {
  objf_val_ = fval;
}

}
