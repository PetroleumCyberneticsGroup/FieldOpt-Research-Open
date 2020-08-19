/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#include "eclsummaryreader.h"
#include "ertwrapper_exceptions.h"
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>

namespace ERTWrapper {
namespace ECLSummary {

using std::string;
using Printer::ext_info;

ECLSummaryReader::ECLSummaryReader(string file_name) {
  file_name_ = file_name;
  ecl_sum_ = ecl_sum_fread_alloc_case(file_name_.c_str(), "");
  if (ecl_sum_ == nullptr) {
    throw SummaryFileNotFoundAtPathException(file_name);
  }
  populateKeyLists();
  initializeVectors();
}

ECLSummaryReader::~ECLSummaryReader() {
  if (ecl_sum_ != nullptr) {
    ecl_sum_free(ecl_sum_);
    ecl_sum_ = nullptr;
  }
}

double ECLSummaryReader::GetMiscVar(const string& var_name,
                                    int time_index) {
  if (!hasMiscVar(var_name)) {
    throw SummaryVariableDoesNotExistException(
      "Misc variable " + string(var_name) + " does not exist.");
  }
  if (!HasReportStep(time_index)) {
    throw SummaryTimeStepDoesNotExistException("Time step does not exist");
  }
  return ecl_sum_get_misc_var(ecl_sum_, time_index, var_name.c_str());
}

double ECLSummaryReader::GetFieldVar(const string& var_name,
                                     int time_index) {
  if (!HasReportStep(time_index)) {
    throw SummaryTimeStepDoesNotExistException("Time step does not exist");
  }
  if (!hasFieldVar(var_name)) {
    throw SummaryVariableDoesNotExistException("Field variable" + var_name + " does not exist.");
  }
  return ecl_sum_get_field_var(ecl_sum_, time_index, var_name.c_str());
}

double ECLSummaryReader::GetWellVar(const string& well_name, string var_name, int time_index) {
  if (!hasWellVar(well_name, var_name)) {
    throw SummaryVariableDoesNotExistException("Well variable " + std::string(well_name) + ":"
                                                 + std::string(var_name) + " does not exist.");
  }
  if (!HasReportStep(time_index)) {
    throw SummaryTimeStepDoesNotExistException("Time step does not exist");
  }
  return ecl_sum_get_well_var(ecl_sum_, time_index, well_name.c_str(), var_name.c_str());
}

int ECLSummaryReader::GetLastReportStep() {
  int last_step = ecl_sum_get_last_report_step(ecl_sum_);
  return ecl_sum_iget_report_end(ecl_sum_, last_step);
}

int ECLSummaryReader::GetFirstReportStep() {
  int first_step = ecl_sum_get_first_report_step(ecl_sum_);
  return ecl_sum_iget_report_start(ecl_sum_, first_step);
}

bool ECLSummaryReader::HasReportStep(int report_step) {
  return report_step <= GetLastReportStep() && report_step >= GetFirstReportStep();
}

bool ECLSummaryReader::hasWellVar(const string& well_name, string var_name) {
  return ecl_sum_has_well_var(ecl_sum_, well_name.c_str(), var_name.c_str());
}

bool ECLSummaryReader::hasGroupVar(const string& group_name, string var_name) {
  return ecl_sum_has_group_var(ecl_sum_, group_name.c_str(), var_name.c_str());
}

bool ECLSummaryReader::hasFieldVar(const string& var_name) {
  return ecl_sum_has_field_var(ecl_sum_, var_name.c_str());
}

bool ECLSummaryReader::hasBlockVar(int block_nr, const string& var_name) {
  return ecl_sum_has_block_var(ecl_sum_, var_name.c_str(), block_nr);
}

bool ECLSummaryReader::hasMiscVar(const string& var_name) {
  return ecl_sum_has_misc_var(ecl_sum_, var_name.c_str());
}

void ECLSummaryReader::populateKeyLists() {
  stringlist_type * keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, nullptr);
  stringlist_type * wells = ecl_sum_alloc_well_list(ecl_sum_, nullptr);
  stringlist_type * field_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "F*");
  stringlist_type * well_keys = ecl_sum_alloc_well_var_list(ecl_sum_);

  for (int i = 0; i < stringlist_get_size(keys); ++i) {
    keys_.insert(stringlist_safe_iget(keys, i));
  }

  for (int j = 0; j < stringlist_get_size(wells); ++j) {
    wells_.insert(stringlist_safe_iget(wells, j));
  }

  for (int k = 0; k < stringlist_get_size(field_keys); ++k) {
    field_keys_.insert(stringlist_safe_iget(field_keys, k));
  }

  for (int l = 0; l < stringlist_get_size(well_keys); ++l) {
    well_keys_.insert(stringlist_safe_iget(well_keys, l));
  }

  std::stringstream ss;
  ss << "Found summary keys: ";
  for (const auto& key : keys_) ss << key << ", ";
  for (const auto& key : field_keys_) ss << key << ", ";
  for (const auto& well : wells_) {
    for (const auto& key : well_keys_) ss << well << ":" << key << ", ";
  }
  if (VERB_SIM >= 2) {
    ext_info("Found summary keys: " + ss.str(), md_, cl_);
  }

  stringlist_free(keys);
  stringlist_free(wells);
  stringlist_free(field_keys);
  stringlist_free(well_keys);
}

void ECLSummaryReader::initializeVectors() {
  initializeTimeVector();
  initializeWellRates();
  initializeWellCumulatives();
  initFieldTotals();
  initFieldRates();

  initVectorsXd();
}

void ECLSummaryReader::initVectorsXd() {
  timeXd_ = Eigen::Map<VectorXd>(time_.data(), time_.size());
  foptXd_ = Eigen::Map<VectorXd>(fopt_.data(), fopt_.size());
  fwptXd_ = Eigen::Map<VectorXd>(fwpt_.data(), fwpt_.size());
  fgptXd_ = Eigen::Map<VectorXd>(fgpt_.data(), fgpt_.size());
  fwitXd_ = Eigen::Map<VectorXd>(fwit_.data(), fwit_.size());
  fgitXd_ = Eigen::Map<VectorXd>(fgit_.data(), fgit_.size());

  foprXd_ = Eigen::Map<VectorXd>(fopr_.data(), fopr_.size());
  fwprXd_ = Eigen::Map<VectorXd>(fwpr_.data(), fwpr_.size());
  fgprXd_ = Eigen::Map<VectorXd>(fgpr_.data(), fgpr_.size());
  fwirXd_ = Eigen::Map<VectorXd>(fwir_.data(), fwir_.size());
  fgirXd_ = Eigen::Map<VectorXd>(fgir_.data(), fgir_.size());
}

void ECLSummaryReader::initializeTimeVector() {
  int days_var_index = ecl_sum_get_misc_var_index(ecl_sum_, "TIME");
  double_vector_type * time = ecl_sum_alloc_data_vector(ecl_sum_, days_var_index, true);
  time_.resize(double_vector_size(time));
  for (int i = 0; i < double_vector_size(time); ++i) {
    time_[i] = (double_vector_safe_iget(time, i));
  }
  time_[0] = GetFirstReportStep();

  step_.resize(double_vector_size(time));
  step_[0] = time_[0];
  for (int i = 1; i < double_vector_size(time); ++i) {
    step_[i] = time_[i] - time_[i-1];
  }
  double_vector_free(time);
}

void ECLSummaryReader::initializeWellRates() {
  const ecl_smspec_type * smspec;
  smspec = ecl_sum_get_smspec(ecl_sum_);
  for(auto wname : wells_) {
    wopr_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WOPR")) {
      int index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WOPR");
      auto * data = ecl_sum_alloc_data_vector(ecl_sum_, index, true);
      assert(double_vector_size(data) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wopr_[wname][i] = double_vector_safe_iget(data, i);
      }
      wopr_[wname][0] = GetWellVar(wname, "WOPR", 0);
      double_vector_free(data);
    }

    wwpr_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WWPR")) {
      int index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WWPR");
      auto * data = ecl_sum_alloc_data_vector(ecl_sum_, index, true);
      assert(double_vector_size(data) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wwpr_[wname][i] = double_vector_safe_iget(data, i);
      }
      wwpr_[wname][0] = GetWellVar(wname, "WWPR", 0);
      double_vector_free(data);
    }

    wgpr_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WGPR")) {
      int index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WGPR");
      auto * data = ecl_sum_alloc_data_vector(ecl_sum_, index, true);
      assert(double_vector_size(data) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wgpr_[wname][i] = double_vector_safe_iget(data, i);
      }
      wgpr_[wname][0] = GetWellVar(wname, "WGPR", 0);
      double_vector_free(data);
    }

    wwir_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WWIR")) {
      int index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WWIR");
      auto * data = ecl_sum_alloc_data_vector(ecl_sum_, index, true);
      assert(double_vector_size(data) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wwir_[wname][i] = double_vector_safe_iget(data, i);
      }
      wwir_[wname][0] = GetWellVar(wname, "WWIR", 0);
      double_vector_free(data);
    }

    wgir_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WGIR")) {
      int index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WGIR");
      auto * data = ecl_sum_alloc_data_vector(ecl_sum_, index, true);
      assert(double_vector_size(data) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wgir_[wname][i] = double_vector_safe_iget(data, i);
      }
      wgir_[wname][0] = GetWellVar(wname, "WGIR", 0);
      double_vector_free(data);
    }
  }
}

void ECLSummaryReader::initializeWellCumulatives() {
  const ecl_smspec_type * smspec = ecl_sum_get_smspec(ecl_sum_);

  for (auto wname : wells_) {

    wopt_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WOPT")) {
      int wopt_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WOPT");
      double_vector_type * wopt = ecl_sum_alloc_data_vector(ecl_sum_, wopt_index, true);
      assert(double_vector_size(wopt) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wopt_[wname][i] = double_vector_safe_iget(wopt, i);
      }
      wopt_[wname][0] = 0.0;
      double_vector_free(wopt);
    } else if (hasWellVar(wname, "WOPR")) {
      if (VERB_SIM >= 2) {
        ext_info("WOPT not found, computing from WOPR.", md_, cl_);
      }
      wopt_[wname] = computeCumulativeFromRate(wopr_[wname]);
    }

    wwpt_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WWPT")) {
      int wwpt_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WWPT");
      double_vector_type * wwpt = ecl_sum_alloc_data_vector(ecl_sum_, wwpt_index, true);
      assert(double_vector_size(wwpt) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wwpt_[wname][i] = double_vector_safe_iget(wwpt, i);
      }
      wwpt_[wname][0] = 0.0;
      double_vector_free(wwpt);
    } else if (hasWellVar(wname, "WWPR")) {
      if (VERB_SIM >= 2) {
        ext_info("WWPT not found, computing from WWPR.", md_, cl_);
      }
      wwpt_[wname] = computeCumulativeFromRate(wwpr_[wname]);
    }

    wgpt_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WGPT")) {
      int wgpt_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WGPT");
      double_vector_type * wgpt = ecl_sum_alloc_data_vector(ecl_sum_, wgpt_index, true);
      assert(double_vector_size(wgpt) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wgpt_[wname][i] = double_vector_safe_iget(wgpt, i);
      }
      wgpt_[wname][0] = 0.0;
      double_vector_free(wgpt);
    } else if (hasWellVar(wname, "WGPR")) {
      if (VERB_SIM >= 2) {
        ext_info("WGPT not found, computing from WGPR.", md_, cl_);
      }
      wgpt_[wname] = computeCumulativeFromRate(wgpr_[wname]);
    }

    wwit_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WWIT")) {
      int wwit_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WWIT");
      double_vector_type * wwit = ecl_sum_alloc_data_vector(ecl_sum_, wwit_index, true);
      assert(double_vector_size(wwit) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wwit_[wname][i] = double_vector_safe_iget(wwit, i);
      }
      wwit_[wname][0] = 0.0;
      double_vector_free(wwit);
    } else if (hasWellVar(wname, "WWIR")) {
      if (VERB_SIM >= 2) {
        ext_info("WWIT not found, computing from WWIR.", md_, cl_);
      }
      wwit_[wname] = computeCumulativeFromRate(wwir_[wname]);
    }

    wgit_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WGIT")) {
      int wgit_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WGIT");
      double_vector_type * wgit = ecl_sum_alloc_data_vector(ecl_sum_, wgit_index, true);
      assert(double_vector_size(wgit) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wgit_[wname][i] = double_vector_safe_iget(wgit, i);
      }
      wgit_[wname][0] = 0.0;
      double_vector_free(wgit);
    } else if (hasWellVar(wname, "WGIR")) {
      if (VERB_SIM >= 2) {
        ext_info("WGIT not found, computing from WGIR.", md_, cl_);
      }
      wgit_[wname] = computeCumulativeFromRate(wgir_[wname]);
    }
  }
}

void ECLSummaryReader::initFieldRates() {
  const ecl_smspec_type * smspec = ecl_sum_get_smspec(ecl_sum_);

  // FOPR
  fopr_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FOPR")) {
    int fopr_index = ecl_smspec_get_field_var_params_index(smspec, "FOPR");
    double_vector_type * fopr = ecl_sum_alloc_data_vector(ecl_sum_, fopr_index, true);
    assert(double_vector_size(fopr) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fopr_[i] = double_vector_safe_iget(fopr, i);
    }
    fopr_[0] = 0.0;
  } else {
    warnPropertyNotFound("FOPR");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fopr_[i] += wopr_[wname][i];
      }
    }
  }

  // FWPR
  fwpr_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FWPR")) {
    int fwpr_index = ecl_smspec_get_field_var_params_index(smspec, "FWPR");
    double_vector_type * fwpr = ecl_sum_alloc_data_vector(ecl_sum_, fwpr_index, true);
    assert(double_vector_size(fwpr) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fwpr_[i] = double_vector_safe_iget(fwpr, i);
    }
    fwpr_[0] = 0.0;
  } else {
    warnPropertyNotFound("FWPR");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fwpr_[i] += wwpr_[wname][i];
      }
    }
  }

  // FGPR
  fgpr_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FGPR")) {
    int fgpr_index = ecl_smspec_get_field_var_params_index(smspec, "FGPR");
    double_vector_type * fgpr = ecl_sum_alloc_data_vector(ecl_sum_, fgpr_index, true);
    assert(double_vector_size(fgpr) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fgpr_[i] = double_vector_safe_iget(fgpr, i);
    }
    fgpr_[0] = 0.0;
  } else {
    warnPropertyNotFound("FGPR");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fgpr_[i] += wgpr_[wname][i];
      }
    }
  }

  // FWIR
  fwir_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FWIR")) {
    int fwir_index = ecl_smspec_get_field_var_params_index(smspec, "FWIR");
    double_vector_type * fwir = ecl_sum_alloc_data_vector(ecl_sum_, fwir_index, true);
    assert(double_vector_size(fwir) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fwir_[i] = double_vector_safe_iget(fwir, i);
    }
    fwir_[0] = 0.0;
  } else {
    warnPropertyNotFound("FWIR");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fwir_[i] += wwir_[wname][i];
      }
    }
  }

  // FGIR
  fgir_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FGIR")) {
    int fgir_index = ecl_smspec_get_field_var_params_index(smspec, "FGIR");
    double_vector_type * fgir = ecl_sum_alloc_data_vector(ecl_sum_, fgir_index, true);
    assert(double_vector_size(fgir) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fgir_[i] = double_vector_safe_iget(fgir, i);
    }
    fgir_[0] = 0.0;
  } else {
    warnPropertyNotFound("FGIR");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fgir_[i] += wgir_[wname][i];
      }
    }
  }

}

void ECLSummaryReader::initFieldTotals() {
  const ecl_smspec_type * smspec = ecl_sum_get_smspec(ecl_sum_);

  // FOPT
  fopt_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FOPT")) {
    int fopt_index = ecl_smspec_get_field_var_params_index(smspec, "FOPT");
    double_vector_type * fopt = ecl_sum_alloc_data_vector(ecl_sum_, fopt_index, true);
    assert(double_vector_size(fopt) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fopt_[i] = double_vector_safe_iget(fopt, i);
    }
    fopt_[0] = 0.0;
  } else {
    warnPropertyNotFound("FOPT");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fopt_[i] += wopt_[wname][i];
      }
    }
  }

  // FWPT
  fwpt_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FWPT")) {
    int fwpt_index = ecl_smspec_get_field_var_params_index(smspec, "FWPT");
    double_vector_type * fwpt = ecl_sum_alloc_data_vector(ecl_sum_, fwpt_index, true);
    assert(double_vector_size(fwpt) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fwpt_[i] = double_vector_safe_iget(fwpt, i);
    }
    fwpt_[0] = 0.0;
  } else {
    warnPropertyNotFound("FWPT");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fwpt_[i] += wwpt_[wname][i];
      }
    }
  }

  // FGPT
  fgpt_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FGPT")) {
    int fgpt_index = ecl_smspec_get_field_var_params_index(smspec, "FGPT");
    double_vector_type * fgpt = ecl_sum_alloc_data_vector(ecl_sum_, fgpt_index, true);
    assert(double_vector_size(fgpt) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fgpt_[i] = double_vector_safe_iget(fgpt, i);
    }
    fgpt_[0] = 0.0;
  } else {
    warnPropertyNotFound("FGPT");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fgpt_[i] += wgpt_[wname][i];
      }
    }
  }

  // FWIT
  fwit_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FWIT")) {
    int fwit_index = ecl_smspec_get_field_var_params_index(smspec, "FWIT");
    double_vector_type * fwit = ecl_sum_alloc_data_vector(ecl_sum_, fwit_index, true);
    assert(double_vector_size(fwit) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fwit_[i] = double_vector_safe_iget(fwit, i);
    }
    fwit_[0] = 0.0;
  } else {
    warnPropertyNotFound("FWIT");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fwit_[i] += wwit_[wname][i];
      }
    }
  }

  // FGIT
  fgit_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FGIT")) {
    int fgit_index = ecl_smspec_get_field_var_params_index(smspec, "FGIT");
    double_vector_type * fgit = ecl_sum_alloc_data_vector(ecl_sum_, fgit_index, true);
    assert(double_vector_size(fgit) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fgit_[i] = double_vector_safe_iget(fgit, i);
    }
    fgit_[0] = 0.0;
  } else {
    warnPropertyNotFound("FGIT");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fgit_[i] += wgit_[wname][i];
      }
    }
  }
}

const std::vector<double> ECLSummaryReader::wopt(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  if (wopt_.at(well_name).back() == 0.0)
    warnPropertyZero(well_name, "WOPT");
  return wopt_.at(well_name);
}

const std::vector<double> ECLSummaryReader::wwpt(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  if (wwpt_.at(well_name).back() == 0.0)
    warnPropertyZero(well_name, "WWPT");
  return wwpt_.at(well_name);
}

const std::vector<double> ECLSummaryReader::wgpt(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  if (wgpt_.at(well_name).back() == 0.0)
    warnPropertyZero(well_name, "WGPT");
  return wgpt_.at(well_name);
}

const std::vector<double> ECLSummaryReader::wwit(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  if (wwit_.at(well_name).back() == 0.0)
    warnPropertyZero(well_name, "WWIT");
  return wwit_.at(well_name);
}

const std::vector<double> ECLSummaryReader::wgit(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  if (wgit_.at(well_name).back() == 0.0)
    warnPropertyZero(well_name, "WGIT");
  return wgit_.at(well_name);
}

void ECLSummaryReader::warnPropertyZero(string wname, string propname) const {
  if (VERB_SIM >= 3) {
    Printer::ext_warn("Returning cumulative vector with final value 0.0 for " + propname + " for well " + wname + ".",
                      "Simulation", "ECLSummaryReader");
  }
}

void ECLSummaryReader::warnPropertyNotFound(string propname) const {
  if (VERB_SIM >= 2) {
    Printer::ext_warn("The property " + propname + " was not found in the summary. Calculating from corresponding well properties.",
                      "Simulation", "ECLSummaryReader");
  }
}

void ECLSummaryReader::warnPropertyZero(string propname) const {
  if (VERB_SIM >= 3) {
    Printer::ext_warn("Returning cumulative vector with final value 0.0 for " + propname + ".",
                      "Simulation", "ECLSummaryReader");
  }
}

const std::vector<double> &ECLSummaryReader::fopt() const {
  if (fopt_.back() == 0.0) { warnPropertyZero("FOPT"); }
  return fopt_;
}

const std::vector<double> &ECLSummaryReader::fwpt() const {
  if (fwpt_.back() == 0.0) { warnPropertyZero("FWPT"); }
  return fwpt_;
}

const std::vector<double> &ECLSummaryReader::fgpt() const {
  if (fgpt_.back() == 0.0) { warnPropertyZero("FGPT"); }
  return fgpt_;
}

const std::vector<double> &ECLSummaryReader::fwit() const {
  if (fgpt_.back() == 0.0) { warnPropertyZero("FWIT"); }
  return fwit_;
}

const std::vector<double> &ECLSummaryReader::fgit() const {
  if (fgit_.back() == 0.0) { warnPropertyZero("FGIT"); }
  return fgit_;
}

const std::vector<double> &ECLSummaryReader::fopr() const {
  if (fopr_.back() == 0.0) { warnPropertyZero("FOPR"); }
  return fopr_;
}

const std::vector<double> &ECLSummaryReader::fwpr() const {
  if (fwpr_.back() == 0.0) { warnPropertyZero("FWPR"); }
  return fwpr_;
}

const std::vector<double> &ECLSummaryReader::fgpr() const {
  if (fgpr_.back() == 0.0) { warnPropertyZero("FGPR"); }
  return fgpr_;
}

const std::vector<double> &ECLSummaryReader::fwir() const {
  if (fgpr_.back() == 0.0) { warnPropertyZero("FWIR"); }
  return fwir_;
}

const std::vector<double> &ECLSummaryReader::fgir() const {
  if (fgir_.back() == 0.0) { warnPropertyZero("FGIR"); }
  return fgir_;
}













const VectorXd &ECLSummaryReader::foptXd() const {
  if (foptXd_.tail(1).value() == 0.0) { warnPropertyZero("FOPT"); }
  return foptXd_;
}

const VectorXd &ECLSummaryReader::fwptXd() const {
  if (fwptXd_.tail(1).value() == 0.0) { warnPropertyZero("FWPT"); }
  return fwptXd_;
}

const VectorXd &ECLSummaryReader::fgptXd() const {
  if (fgptXd_.tail(1).value() == 0.0) { warnPropertyZero("FGPT"); }
  return fgptXd_;
}

const VectorXd &ECLSummaryReader::fwitXd() const {
  if (fgptXd_.tail(1).value() == 0.0) { warnPropertyZero("FWIT"); }
  return fwitXd_;
}

const VectorXd &ECLSummaryReader::fgitXd() const {
  if (fgitXd_.tail(1).value() == 0.0) { warnPropertyZero("FGIT"); }
  return fgitXd_;
}

const VectorXd &ECLSummaryReader::foprXd() const {
  if (foprXd_.tail(1).value() == 0.0) { warnPropertyZero("FOPR"); }
  return foprXd_;
}

const VectorXd &ECLSummaryReader::fwprXd() const {
  if (fwprXd_.tail(1).value() == 0.0) { warnPropertyZero("FWPR"); }
  return fwprXd_;
}

const VectorXd &ECLSummaryReader::fgprXd() const {
  if (fgprXd_.tail(1).value() == 0.0) { warnPropertyZero("FGPR"); }
  return fgprXd_;
}

const VectorXd &ECLSummaryReader::fwirXd() const {
  if (fgprXd_.tail(1).value() == 0.0) { warnPropertyZero("FWIR"); }
  return fwirXd_;
}

const VectorXd &ECLSummaryReader::fgirXd() const {
  if (fgirXd_.tail(1).value() == 0.0) { warnPropertyZero("FGIR"); }
  return fgirXd_;
}














std::vector<double> ECLSummaryReader::wopr(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  return wopr_.at(well_name);
}

std::vector<double> ECLSummaryReader::wwpr(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  return wwpr_.at(well_name);
}

std::vector<double> ECLSummaryReader::wgpr(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  return wgpr_.at(well_name);
}

std::vector<double> ECLSummaryReader::wwir(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  return wwir_.at(well_name);
}

std::vector<double> ECLSummaryReader::wgir(const string well_name) const {
  if (wells_.find(well_name) == wells_.end())
    throw SummaryVariableDoesNotExistException("The well " + well_name + " was not found in the summary.");
  return wgir_.at(well_name);
}
vector<double> ECLSummaryReader::computeCumulativeFromRate(vector<double> rate) {
  assert(time_.size() == rate.size());
  auto cumulative = vector<double>(rate.size(), 0.0);
  for (int i = 1; i < rate.size(); ++i) {
    double dt = time_[i] - time_[i-1];
    cumulative[i] = dt * rate[i-1];
  }
  return cumulative;
}

}
}
