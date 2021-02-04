/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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
using std::modf;
using Printer::info;
using Printer::ext_info;
using Printer::num2str;
using Printer::ext_warn;

ECLSummaryReader::ECLSummaryReader(const string& file_name,
                                   Settings::Settings *settings) {
  settings_ = settings;
  vp_ = settings_->global()->verbParams();
  // settings_->global()->showVerbParams(vp_, "[ @ECLSummaryReader ]");

  file_name_ = file_name;
  ecl_sum_ = ecl_sum_fread_alloc_case(file_name_.c_str(), "");
  if (ecl_sum_ == nullptr) {
    throw SmryFileNotFoundAtPathExc(file_name);
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
    string em = "Misc variable " + var_name + " does not exist.";
    throw SmryVarDoesNotExistExc(em);
  }
  if (!HasReportStep(time_index)) {
    throw SmryTimeStepDoesNotExistExc("Time step does not exist");
  }
  return ecl_sum_get_misc_var(ecl_sum_, time_index, var_name.c_str());
}

double ECLSummaryReader::GetFieldVar(const string& var_name,
                                     int time_index) {
  if (!HasReportStep(time_index)) {
    throw SmryTimeStepDoesNotExistExc("Time step does not exist");
  }
  if (!hasFieldVar(var_name)) {
    throw SmryVarDoesNotExistExc("Field variable" + var_name + " does not exist.");
  }
  return ecl_sum_get_field_var(ecl_sum_, time_index, var_name.c_str());
}

double ECLSummaryReader::GetWellVar(const string& well_name,
                                    string var_name, int time_index) {
  if (!hasWellVar(well_name, var_name)) {
    string em = "Well variable " + std::string(well_name) + ":";
    em += std::string(var_name) + " does not exist.";
    throw SmryVarDoesNotExistExc(em);
  }
  if (!HasReportStep(time_index)) {
    throw SmryTimeStepDoesNotExistExc("Time step does not exist");
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

  stringlist_type * seg_sofr_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "SOFR*");
  stringlist_type * seg_swfr_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "SWFR*");
  stringlist_type * seg_sgfr_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "SGFR*");
  stringlist_type * seg_spr_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "SPR*");
  stringlist_type * seg_sprd_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "SPRD*");
  stringlist_type * seg_swct_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "SWCT*");
  stringlist_type * seg_scsa_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "SCSA*");
  stringlist_type * comp_keys = ecl_sum_alloc_matching_general_var_list(ecl_sum_, "C*");

  for (int l = 0; l < stringlist_get_size(seg_sofr_keys); ++l) {
    seg_sofr_keys_.insert(stringlist_safe_iget(seg_sofr_keys, l));
  }

  for (int l = 0; l < stringlist_get_size(seg_swfr_keys); ++l) {
    seg_swfr_keys_.insert(stringlist_safe_iget(seg_swfr_keys, l));
  }

  for (int l = 0; l < stringlist_get_size(seg_sgfr_keys); ++l) {
    seg_sgfr_keys_.insert(stringlist_safe_iget(seg_sgfr_keys, l));
  }

  for (int l = 0; l < stringlist_get_size(seg_spr_keys); ++l) {
    seg_sprp_keys_.insert(stringlist_safe_iget(seg_spr_keys, l));
  }

  for (int l = 0; l < stringlist_get_size(seg_sprd_keys); ++l) {
    seg_sprd_keys_.insert(stringlist_safe_iget(seg_sprd_keys, l));
  }

  for (int l = 0; l < stringlist_get_size(seg_swct_keys); ++l) {
    seg_swct_keys_.insert(stringlist_safe_iget(seg_swct_keys, l));
  }

  for (int l = 0; l < stringlist_get_size(seg_scsa_keys); ++l) {
    seg_scsa_keys_.insert(stringlist_safe_iget(seg_scsa_keys, l));
  }

  for (int l = 0; l < stringlist_get_size(comp_keys); ++l) {
    comp_keys_.insert(stringlist_safe_iget(comp_keys, l));
  }

  if (vp_.vSIM >= 5) {
    std::stringstream skeys, sfield, swells, sseg, scomp;
    int ii = 0;

    skeys << "All summary keys: ";
    for (const auto& skey : keys_) { skeys << skey << ", "; }
    // ext_info(skeys.str(), md_, cl_, vp_.lnw);

    sfield << "Field keys: ";
    for (const auto& fkey : field_keys_) { sfield << fkey << ", "; }
    ext_info(sfield.str(), md_, cl_, vp_.lnw);

    swells << "Well keys: ";
    for (const auto& well : wells_) {
      for (const auto& wkey : well_keys_) {
        swells << well << ":" << wkey << ", ";
      }
    }
    ext_info(swells.str(), md_, cl_, vp_.lnw);

    sseg << "Segment SOFR keys: ";
    for (const auto& ckey : seg_sofr_keys_) {
      if(sseg.str().size() % 130 <= 10) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    }
    ext_info(sseg.str(), md_, cl_, vp_.lnw);
    sseg.str("");

    sseg << "Segment SWFR keys: ";
    for (const auto& ckey : seg_swfr_keys_) {
      if(sseg.str().size() % 130 <= 10) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    }
    ext_info(sseg.str(), md_, cl_, vp_.lnw);
    sseg.str("");

    sseg << "Segment SGFR keys: ";
    for (const auto& ckey : seg_sgfr_keys_) {
      if(sseg.str().size() % 130 <= 10) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    }
    ext_info(sseg.str(), md_, cl_, vp_.lnw);
    sseg.str("");

    sseg << "Segment SOFR keys: ";
    for (const auto& ckey : seg_sofr_keys_) {
      if(sseg.str().size() % 130 <= 10) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    }
    ext_info(sseg.str(), md_, cl_, vp_.lnw);
    sseg.str("");

    sseg << "Segment SPR keys: ";
    for (const auto& ckey : seg_sprp_keys_) {
      if(sseg.str().size() % 130 <= 20) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    }
    ext_info(sseg.str(), md_, cl_, vp_.lnw);
    sseg.str("");

    sseg << "Segment SPRD keys: ";
    for (const auto& ckey : seg_sprd_keys_) {
      if(sseg.str().size() % 130 <= 10) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    }
    ext_info(sseg.str(), md_, cl_, vp_.lnw);
    sseg.str("");

    sseg << "Segment SWCT keys: ";
    for (const auto& ckey : seg_swct_keys_) {
      if(sseg.str().size() % 130 <= 10) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    }
    ext_info(sseg.str(), md_, cl_, vp_.lnw);
    sseg.str("");

    sseg << "Segment SCSA keys: ";
    for (const auto& ckey : seg_scsa_keys_) {
      if(sseg.str().size() % 130 <= 10) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    }
    ext_info(sseg.str(), md_, cl_, vp_.lnw);
    sseg.str("");

    // sseg << "Compartment keys: ";
    // for (const auto& ckey : comp_keys_) {
    //   if(sseg.str().size() % 130 <= 15) { sseg << ckey << " "; } else { sseg << ckey << "_"; }
    // }
    // ext_info(sseg.str(), md_, cl_, vp_.lnw);
    // sseg.str("");
  }

  if (!seg_sofr_keys_.empty() && !seg_swfr_keys_.empty() && !seg_sgfr_keys_.empty()
    && !seg_sprp_keys_.empty() && !seg_sprd_keys_.empty() && !seg_swct_keys_.empty()
    && !seg_scsa_keys_.empty()) {
    read_segments_ = true;
  } else {
    ext_warn("Some or all segment result keywords empty - not read.", md_, cl_, vp_.lnw);
  }

  stringlist_free(keys);
  stringlist_free(wells);
  stringlist_free(field_keys);
  stringlist_free(well_keys);
}

void ECLSummaryReader::initializeVectors() {
  initializeTimeVector();
  initializeWellRates();
  initWellTotals();
  initFieldTotals();
  initFieldRates();

  initVectorsXd();
  if(read_segments_) {
    initWellSegRates();
    initVectorsSegXd();
  }
}

void ECLSummaryReader::initVectorsXd() {
  foptXd_ = Eigen::Map<VectorXd>(fopt_.data(), fopt_.size());
  fwptXd_ = Eigen::Map<VectorXd>(fwpt_.data(), fwpt_.size());
  fgptXd_ = Eigen::Map<VectorXd>(fgpt_.data(), fgpt_.size());
  flptXd_ = Eigen::Map<VectorXd>(flpt_.data(), flpt_.size());
  fwitXd_ = Eigen::Map<VectorXd>(fwit_.data(), fwit_.size());
  fgitXd_ = Eigen::Map<VectorXd>(fgit_.data(), fgit_.size());

  for (const auto& wname : wells_) {
    woptXd_[wname] = Eigen::Map<VectorXd>(wopt_[wname].data(), wopt_[wname].size());
    wwptXd_[wname] = Eigen::Map<VectorXd>(wwpt_[wname].data(), wwpt_[wname].size());
    wgptXd_[wname] = Eigen::Map<VectorXd>(wgpt_[wname].data(), wgpt_[wname].size());
    wlptXd_[wname] = Eigen::Map<VectorXd>(wlpt_[wname].data(), wlpt_[wname].size());
    wwitXd_[wname] = Eigen::Map<VectorXd>(wwit_[wname].data(), wwit_[wname].size());
    wgitXd_[wname] = Eigen::Map<VectorXd>(wgit_[wname].data(), wgit_[wname].size());

    // cout << "Well: " << wname << endl;
    // cout << "wopt[wn].sz(): " << wopt_[wname].size() << " -- woptXd[wn].sz(): " << woptXd_[wname].size() << endl;
    // cout << "wgpt[wn].sz(): " << wgpt_[wname].size() << " -- wgptXd[wn].sz(): " << wgptXd_[wname].size() << endl;
    // cout << "wwpt[wn].sz(): " << wwpt_[wname].size() << " -- wwptXd[wn].sz(): " << wwptXd_[wname].size() << endl;
    // cout << "wwit[wn].sz(): " << wwit_[wname].size() << " -- wwitXd[wn].sz(): " << wwitXd_[wname].size() << endl;
    // cout << "wgit[wn].sz(): " << wgit_[wname].size() << " -- wgitXd[wn].sz(): " << wgitXd_[wname].size() << endl;
  }

  timeXd_ = Eigen::Map<VectorXd>(time_.data(), time_.size());
  stepXd_ = Eigen::Map<VectorXd>(step_.data(), step_.size());

  foprXd_ = Eigen::Map<VectorXd>(fopr_.data(), fopr_.size());
  fwprXd_ = Eigen::Map<VectorXd>(fwpr_.data(), fwpr_.size());
  fgprXd_ = Eigen::Map<VectorXd>(fgpr_.data(), fgpr_.size());
  flprXd_ = Eigen::Map<VectorXd>(flpr_.data(), flpr_.size());
  fwirXd_ = Eigen::Map<VectorXd>(fwir_.data(), fwir_.size());
  fgirXd_ = Eigen::Map<VectorXd>(fgir_.data(), fgir_.size());

  for (const auto& wname : wells_) {
    woprXd_[wname] = Eigen::Map<VectorXd>(wopr_[wname].data(), wopr_[wname].size());
    wwprXd_[wname] = Eigen::Map<VectorXd>(wwpr_[wname].data(), wwpr_[wname].size());
    wgprXd_[wname] = Eigen::Map<VectorXd>(wgpr_[wname].data(), wgpr_[wname].size());
    wlprXd_[wname] = Eigen::Map<VectorXd>(wlpr_[wname].data(), wlpr_[wname].size());
    wwirXd_[wname] = Eigen::Map<VectorXd>(wwir_[wname].data(), wwir_[wname].size());
    wgirXd_[wname] = Eigen::Map<VectorXd>(wgir_[wname].data(), wgir_[wname].size());

    wbhpXd_[wname] = Eigen::Map<VectorXd>(wbhp_[wname].data(), wbhp_[wname].size());
    wwctXd_[wname] = Eigen::Map<VectorXd>(wwct_[wname].data(), wwct_[wname].size());
  }
}

void ECLSummaryReader::initVectorsSegXd() {
  for (const auto& wname : wells_) {
    for (int ii=0; ii < seg_sofr_[wname].size(); ++ii) {

      if (!seg_sofr_[wname].empty() && !seg_swfr_[wname].empty()
        &&!seg_sprp_[wname].empty() && !seg_sprd_[wname].empty()
        &&!seg_swct_[wname].empty() && !seg_scsa_[wname].empty()) {

        if (vp_.vSIM >= 5) { info("Init SegXd vectors -> Well:" + wname + ":Seg#" + num2str(ii, 0), vp_.lnw); }

        seg_sofrXd_[wname].push_back(Eigen::Map<VectorXd>(seg_sofr_[wname][ii].data(), seg_sofr_[wname][ii].size()));
        seg_swfrXd_[wname].push_back(Eigen::Map<VectorXd>(seg_swfr_[wname][ii].data(), seg_swfr_[wname][ii].size()));
        seg_sgfrXd_[wname].push_back(Eigen::Map<VectorXd>(seg_sgfr_[wname][ii].data(), seg_sgfr_[wname][ii].size()));
        seg_sprpXd_[wname].push_back(Eigen::Map<VectorXd>(seg_sprp_[wname][ii].data(), seg_sprp_[wname][ii].size()));
        seg_sprdXd_[wname].push_back(Eigen::Map<VectorXd>(seg_sprd_[wname][ii].data(), seg_sprd_[wname][ii].size()));
        seg_swctXd_[wname].push_back(Eigen::Map<VectorXd>(seg_swct_[wname][ii].data(), seg_swct_[wname][ii].size()));
        seg_scsaXd_[wname].push_back(Eigen::Map<VectorXd>(seg_scsa_[wname][ii].data(), seg_scsa_[wname][ii].size()));

        seg_softXd_[wname].push_back(Eigen::Map<VectorXd>(seg_soft_[wname][ii].data(), seg_soft_[wname][ii].size()));
        seg_swftXd_[wname].push_back(Eigen::Map<VectorXd>(seg_swft_[wname][ii].data(), seg_swft_[wname][ii].size()));
        seg_sgftXd_[wname].push_back(Eigen::Map<VectorXd>(seg_sgft_[wname][ii].data(), seg_sgft_[wname][ii].size()));
        seg_slftXd_[wname].push_back(Eigen::Map<VectorXd>(seg_slft_[wname][ii].data(), seg_slft_[wname][ii].size()));

      }
    }
  }
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

    wlpr_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WLPR")) {
      int index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WLPR");
      auto * data = ecl_sum_alloc_data_vector(ecl_sum_, index, true);
      assert(double_vector_size(data) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wlpr_[wname][i] = double_vector_safe_iget(data, i);
      }
      wlpr_[wname][0] = GetWellVar(wname, "WLPR", 0);
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

    wbhp_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WBHP")) {
      int index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WBHP");
      auto * data = ecl_sum_alloc_data_vector(ecl_sum_, index, true);
      assert(double_vector_size(data) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wbhp_[wname][i] = double_vector_safe_iget(data, i);
      }
      wbhp_[wname][0] = GetWellVar(wname, "WBHP", 0);
      double_vector_free(data);
    }

    wwct_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WWCT")) {
      int index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WWCT");
      auto * data = ecl_sum_alloc_data_vector(ecl_sum_, index, true);
      assert(double_vector_size(data) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wwct_[wname][i] = double_vector_safe_iget(data, i);
      }
      wwct_[wname][0] = GetWellVar(wname, "WWCT", 0);
      double_vector_free(data);
    }

  }
}

void ECLSummaryReader::initWellTotals() {
  const ecl_smspec_type * smspec = ecl_sum_get_smspec(ecl_sum_);

  for (auto wname : wells_) {

    // WOPT
    wopt_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WOPT")) {
      int wopt_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WOPT");
      double_vector_type * wopt = ecl_sum_alloc_data_vector(ecl_sum_, wopt_index,true);
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

    // WWPT
    wwpt_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WWPT")) {
      int wwpt_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WWPT");
      double_vector_type * wwpt = ecl_sum_alloc_data_vector(ecl_sum_, wwpt_index,true);
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

    // WGPT
    wgpt_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WGPT")) {
      int wgpt_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WGPT");
      double_vector_type * wgpt = ecl_sum_alloc_data_vector(ecl_sum_, wgpt_index,true);
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

    // WLPT
    wlpt_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WLPT")) {
      int wlpt_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WLPT");
      double_vector_type * wlpt = ecl_sum_alloc_data_vector(ecl_sum_, wlpt_index,true);
      assert(double_vector_size(wlpt) == time_.size());
      for (int i = 0; i < time_.size(); ++i) {
        wlpt_[wname][i] = double_vector_safe_iget(wlpt, i);
      }
      wlpt_[wname][0] = 0.0;
      double_vector_free(wlpt);
    } else if (hasWellVar(wname, "WLPT")) {
      if (VERB_SIM >= 2) {
        ext_info("WLPT not found, computing from WLPT.", md_, cl_);
      }
      wlpt_[wname] = computeCumulativeFromRate(wlpr_[wname]);
    }

    // WWIT
    wwit_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WWIT")) {
      int wwit_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WWIT");
      double_vector_type * wwit = ecl_sum_alloc_data_vector(ecl_sum_, wwit_index,true);
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

    // WGIT
    wgit_[wname] = std::vector<double>(time_.size(), 0.0);
    if (hasWellVar(wname, "WGIT")) {
      int wgit_index = ecl_smspec_get_well_var_params_index(smspec, wname.c_str(), "WGIT");
      double_vector_type * wgit = ecl_sum_alloc_data_vector(ecl_sum_, wgit_index,true);
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
    double_vector_type * fopr = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fopr_index,
                                                          true);
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
    double_vector_type * fwpr = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fwpr_index,
                                                          true);
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
    double_vector_type * fgpr = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fgpr_index,
                                                          true);
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

  // FLPR
  flpr_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FLPR")) {
    int flpr_index = ecl_smspec_get_field_var_params_index(smspec, "FLPR");
    double_vector_type * flpr = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          flpr_index,
                                                          true);
    assert(double_vector_size(flpr) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      flpr_[i] = double_vector_safe_iget(flpr, i);
    }
    flpr_[0] = 0.0;
  } else {
    warnPropertyNotFound("FLPR");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        flpr_[i] += wlpr_[wname][i];
      }
    }
  }

  // FWIR
  fwir_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FWIR")) {
    int fwir_index = ecl_smspec_get_field_var_params_index(smspec, "FWIR");
    double_vector_type * fwir = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fwir_index,
                                                          true);
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
    double_vector_type * fgir = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fgir_index,
                                                          true);
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
    double_vector_type * fopt = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fopt_index,
                                                          true);
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
    double_vector_type * fwpt = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fwpt_index,
                                                          true);
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
    double_vector_type * fgpt = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fgpt_index,
                                                          true);
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

  // FLPT
  flpt_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FLPT")) {
    int flpt_index = ecl_smspec_get_field_var_params_index(smspec, "FLPT");
    double_vector_type * flpt = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          flpt_index,
                                                          true);
    assert(double_vector_size(flpt) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      flpt_[i] = double_vector_safe_iget(flpt, i);
    }
    flpt_[0] = 0.0;
  } else {
    warnPropertyNotFound("FLPT");
    for (auto wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        flpt_[i] += wlpt_[wname][i];
      }
    }
  }

  // FWIT
  fwit_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FWIT")) {
    int fwit_index = ecl_smspec_get_field_var_params_index(smspec, "FWIT");
    double_vector_type * fwit = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fwit_index,
                                                          true);
    assert(double_vector_size(fwit) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fwit_[i] = double_vector_safe_iget(fwit, i);
    }
    fwit_[0] = 0.0;
  } else {
    warnPropertyNotFound("FWIT");
    for (const auto& wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fwit_[i] += wwit_[wname][i];
      }
    }
  }

  // FGIT
  fgit_ = std::vector<double>(time_.size(), 0.0);
  if (hasFieldVar("FGIT")) {
    int fgit_index = ecl_smspec_get_field_var_params_index(smspec, "FGIT");
    double_vector_type * fgit = ecl_sum_alloc_data_vector(ecl_sum_,
                                                          fgit_index,
                                                          true);
    assert(double_vector_size(fgit) == time_.size());
    for (int i = 0; i < time_.size(); ++i) {
      fgit_[i] = double_vector_safe_iget(fgit, i);
    }
    fgit_[0] = 0.0;
  } else {
    warnPropertyNotFound("FGIT");
    for (const auto& wname : wells_) {
      for (int i = 0; i < time_.size(); ++i) {
        fgit_[i] += wgit_[wname][i];
      }
    }
  }
}

// ===================================================== STD

// Field vectors TOTALS (std::vector) ----------------------
const vector<double> &ECLSummaryReader::fopt() const {
  if (fopt_.back() == 0.0) { warnPropertyZero("FOPT"); }
  return fopt_;
}

const vector<double> &ECLSummaryReader::fwpt() const {
  if (fwpt_.back() == 0.0) { warnPropertyZero("FWPT"); }
  return fwpt_;
}

const vector<double> &ECLSummaryReader::fgpt() const {
  if (fgpt_.back() == 0.0) { warnPropertyZero("FGPT"); }
  return fgpt_;
}

const vector<double> &ECLSummaryReader::flpt() const {
  if (flpt_.back() == 0.0) { warnPropertyZero("FLPT"); }
  return flpt_;
}

const vector<double> &ECLSummaryReader::fwit() const {
  if (fgpt_.back() == 0.0) { warnPropertyZero("FWIT"); }
  return fwit_;
}

const vector<double> &ECLSummaryReader::fgit() const {
  if (fgit_.back() == 0.0) { warnPropertyZero("FGIT"); }
  return fgit_;
}

// Well vectors TOTALS (std::vector) -----------------------

const vector<double> &ECLSummaryReader::wopt(string wname) {
  if (wopt_[wname].back() == 0.0) { warnPropertyZero(wname, "WOPT"); }
  return wopt_[wname];
}

const vector<double> &ECLSummaryReader::wwpt(string wname) {
  if (wwpt_[wname].back() == 0.0) { warnPropertyZero(wname, "WWPT"); }
  return wwpt_[wname];
}

const vector<double> &ECLSummaryReader::wgpt(string wname) {
  if (wgpt_[wname].back() == 0.0) { warnPropertyZero(wname, "WGPT"); }
  return wgpt_[wname];
}

const vector<double> &ECLSummaryReader::wlpt(string wname) {
  if (wlpt_[wname].back() == 0.0) { warnPropertyZero(wname, "WLPT"); }
  return wlpt_[wname];
}

const vector<double> &ECLSummaryReader::wwit(string wname) {
  if (wwit_[wname].back() == 0.0) { warnPropertyZero(wname, "WWIT"); }
  return wwit_[wname];
}

const vector<double> &ECLSummaryReader::wgit(string wname) {
  if (wgit_[wname].back() == 0.0) { warnPropertyZero(wname, "WGIT"); }
  return wgit_[wname];
}

const vector<double> &ECLSummaryReader::wbhp(string wname) {
  if (wbhp_[wname].back() == 0.0) { warnPropertyZero(wname, "WBHP"); }
  return wbhp_[wname];
}

const vector<double> &ECLSummaryReader::wwct(string wname) {
  if (wwct_[wname].back() == 0.0) { warnPropertyZero(wname, "WWCT"); }
  return wwct_[wname];
}

// Field vectors RATES (std::vector) -----------------------

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

const std::vector<double> &ECLSummaryReader::flpr() const {
  if (flpr_.back() == 0.0) { warnPropertyZero("FLPR"); }
  return flpr_;
}

const std::vector<double> &ECLSummaryReader::fwir() const {
  if (fgpr_.back() == 0.0) { warnPropertyZero("FWIR"); }
  return fwir_;
}

const std::vector<double> &ECLSummaryReader::fgir() const {
  if (fgir_.back() == 0.0) { warnPropertyZero("FGIR"); }
  return fgir_;
}

// Field vectors RATES (std::vector) -----------------------

const vector<double> &ECLSummaryReader::wopr(string wname) {
  if (wopr_[wname].back() == 0.0) { warnPropertyZero(wname, "WOPR"); }
  return wopr_[wname];
}

const vector<double> &ECLSummaryReader::wwpr(string wname) {
  if (wwpr_[wname].back() == 0.0) { warnPropertyZero(wname, "WWPR"); }
  return wwpr_[wname];
}

const vector<double> &ECLSummaryReader::wgpr(string wname) {
  if (wgpr_[wname].back() == 0.0) { warnPropertyZero(wname, "WGPR"); }
  return wgpr_[wname];
}

const vector<double> &ECLSummaryReader::wlpr(string wname) {
  if (wlpr_[wname].back() == 0.0) { warnPropertyZero(wname, "WLPR"); }
  return wlpr_[wname];
}

const vector<double> &ECLSummaryReader::wwir(string wname) {
  if (wwir_[wname].back() == 0.0) { warnPropertyZero(wname, "WWIR"); }
  return wwir_[wname];
}

const vector<double> &ECLSummaryReader::wgir(string wname) {
  if (wgir_[wname].back() == 0.0) { warnPropertyZero(wname, "WGIR"); }
  return wgir_[wname];
}

// Segment vectors RATES (std::vector) -----------------------

const vector<vector<double> > &ECLSummaryReader::seg_sofr(string wname) {
  for (int jj=0; jj < seg_sofr_[wname].size(); ++jj) {
    if (seg_sofr_[wname][jj].empty()) {
      warnPropertyZero("SOFR:" + wname + ":" + num2str(jj));
    }
  }
  return seg_sofr_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_swfr(string wname) {
  for (int ii=0; ii < seg_swfr_[wname].size(); ++ii) {
    if (seg_swfr_[wname][ii].empty()) {
      warnPropertyZero("SWFR:" + wname + ":" + num2str(ii));
    }
  }
  return seg_swfr_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_sgfr(string wname) {
  for (int ii=0; ii < seg_sgfr_[wname].size(); ++ii) {
    if (seg_sgfr_[wname][ii].empty()) {
      warnPropertyZero("SGFR:" + wname + ":" + num2str(ii));
    }
  }
  return seg_sgfr_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_slfr(string wname) {
  for (int ii=0; ii < seg_slfr_[wname].size(); ++ii) {
    if (seg_slfr_[wname][ii].empty()) {
      warnPropertyZero("SLFR:" + wname + ":" + num2str(ii));
    }
  }
  return seg_slfr_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_sprp(string wname) {
  for (int ii=0; ii < seg_sprp_[wname].size(); ++ii) {
    if (seg_sprp_[wname][ii].empty()) {
      warnPropertyZero("SPRP:" + wname + ":" + num2str(ii));
    }
  }
  return seg_sprp_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_sprd(string wname) {
  for (int ii=0; ii < seg_sprd_[wname].size(); ++ii) {
    if (seg_sprd_[wname][ii].empty()) {
      warnPropertyZero("SPRD:" + wname + ":" + num2str(ii));
    }
  }
  return seg_sprd_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_swct(string wname) {
  for (int ii=0; ii < seg_swct_[wname].size(); ++ii) {
    if (seg_swct_[wname][ii].empty()) {
      warnPropertyZero("SWCT:" + wname + ":" + num2str(ii));
    }
  }
  return seg_swct_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_scsa(string wname) {
  for (int ii=0; ii < seg_scsa_[wname].size(); ++ii) {
    if (seg_scsa_[wname][ii].empty()) {
      warnPropertyZero("SCSA:" + wname + ":" + num2str(ii));
    }
  }
  return seg_scsa_[wname];
}

// Segment vectors TOTALS (std::vector) -----------------------
const vector<vector<double> > &ECLSummaryReader::seg_soft(string wname) {
  for (int jj=0; jj < seg_soft_[wname].size(); ++jj) {
    if (seg_soft_[wname][jj].empty()) {
      warnPropertyZero("SOFT:" + wname + ":" + num2str(jj));
    }
  }
  return seg_soft_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_swft(string wname) {
  for (int jj=0; jj < seg_swft_[wname].size(); ++jj) {
    if (seg_swft_[wname][jj].empty()) {
      warnPropertyZero("SWFT:" + wname + ":" + num2str(jj));
    }
  }
  return seg_swft_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_sgft(string wname) {
  for (int jj=0; jj < seg_sgft_[wname].size(); ++jj) {
    if (seg_sgft_[wname][jj].empty()) {
      warnPropertyZero("SGFT:" + wname + ":" + num2str(jj));
    }
  }
  return seg_sgft_[wname];
}

const vector<vector<double> > &ECLSummaryReader::seg_slft(string wname) {
  for (int jj=0; jj < seg_slft_[wname].size(); ++jj) {
    if (seg_slft_[wname][jj].empty()) {
      warnPropertyZero("SLFT:" + wname + ":" + num2str(jj));
    }
  }
  return seg_slft_[wname];
}

// ================================================ VECTORXD

// Field vectors TOTALS (VectorXd) -------------------------
const VectorXd &ECLSummaryReader::foptXd() const {
  if (foptXd_.tail(1).value() == 0.0) { warnPropertyZero("FOPTXd"); }
  return foptXd_;
}

const VectorXd &ECLSummaryReader::fwptXd() const {
  if (fwptXd_.tail(1).value() == 0.0) { warnPropertyZero("FWPTXd"); }
  return fwptXd_;
}

const VectorXd &ECLSummaryReader::fgptXd() const {
  if (fgptXd_.tail(1).value() == 0.0) { warnPropertyZero("FGPTXd"); }
  return fgptXd_;
}

const VectorXd &ECLSummaryReader::flptXd() const {
  if (flptXd_.tail(1).value() == 0.0) { warnPropertyZero("FLPTXd"); }
  return flptXd_;
}

const VectorXd &ECLSummaryReader::fwitXd() const {
  if (fgptXd_.tail(1).value() == 0.0) { warnPropertyZero("FWITXd"); }
  return fwitXd_;
}

const VectorXd &ECLSummaryReader::fgitXd() const {
  if (fgitXd_.tail(1).value() == 0.0) { warnPropertyZero("FGITXd"); }
  return fgitXd_;
}

// Well vectors TOTALS (VectorXd) -------------------------

const VectorXd &ECLSummaryReader::woptXd(string wname) {
  if (wells_.find(wname) == wells_.end()) {
    string em = "Well " + wname + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return woptXd_[wname];
}

const VectorXd &ECLSummaryReader::wwptXd(string wname) {
  if (wells_.find(wname) == wells_.end()) {
    string em = "Well " + wname + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wwptXd_[wname];
}

const VectorXd &ECLSummaryReader::wgptXd(string wname) {
  if (wells_.find(wname) == wells_.end()) {
    string em = "Well " + wname + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wgptXd_[wname];
}

const VectorXd &ECLSummaryReader::wlptXd(string wname) {
  if (wells_.find(wname) == wells_.end()) {
    string em = "Well " + wname + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wlptXd_[wname];
}

const VectorXd &ECLSummaryReader::wwitXd(string wname) {
  if (wells_.find(wname) == wells_.end()) {
    string em = "Well " + wname + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wwitXd_[wname];
}

const VectorXd &ECLSummaryReader::wgitXd(string wname) {
  if (wells_.find(wname) == wells_.end()) {
    string em = "Well " + wname + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wgitXd_[wname];
}

const VectorXd &ECLSummaryReader::wbhpXd(string wname) {
  if (wells_.find(wname) == wells_.end()) {
    string em = "Well " + wname + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wbhpXd_[wname];
}

const VectorXd &ECLSummaryReader::wwctXd(string wname) {
  if (wells_.find(wname) == wells_.end()) {
    string em = "Well " + wname + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wwctXd_[wname];
}

// Field vectors RATES (VectorXd) -------------------------

const VectorXd &ECLSummaryReader::foprXd() const {
  if (foprXd_.tail(1).value() == 0.0) { warnPropertyZero("FOPRXd"); }
  return foprXd_;
}

const VectorXd &ECLSummaryReader::fwprXd() const {
  if (fwprXd_.tail(1).value() == 0.0) { warnPropertyZero("FWPRXd"); }
  return fwprXd_;
}

const VectorXd &ECLSummaryReader::fgprXd() const {
  if (fgprXd_.tail(1).value() == 0.0) { warnPropertyZero("FGPRXd"); }
  return fgprXd_;
}

const VectorXd &ECLSummaryReader::flprXd() const {
  if (flprXd_.tail(1).value() == 0.0) { warnPropertyZero("FLPRXd"); }
  return flprXd_;
}

const VectorXd &ECLSummaryReader::fwirXd() const {
  if (fgprXd_.tail(1).value() == 0.0) { warnPropertyZero("FWIRXd"); }
  return fwirXd_;
}

const VectorXd &ECLSummaryReader::fgirXd() const {
  if (fgirXd_.tail(1).value() == 0.0) { warnPropertyZero("FGIRXd"); }
  return fgirXd_;
}

// Well vectors RATES (VectorXd) -------------------------

const VectorXd &ECLSummaryReader::woprXd(string well_name) {
  if (wells_.find(well_name) == wells_.end()) {
    string em = "Well " + well_name + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return woprXd_[well_name];
}

const VectorXd &ECLSummaryReader::wwprXd(string well_name) {
  if (wells_.find(well_name) == wells_.end()) {
    string em = "The well " + well_name + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wwprXd_[well_name];
}

const VectorXd &ECLSummaryReader::wgprXd(string well_name) {
  if (wells_.find(well_name) == wells_.end()) {
    string em = "Well " + well_name + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wgprXd_[well_name];
}

const VectorXd &ECLSummaryReader::wlprXd(string well_name) {
  if (wells_.find(well_name) == wells_.end()) {
    string em = "Well " + well_name + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wlprXd_[well_name];
}

const VectorXd &ECLSummaryReader::wwirXd(string well_name) {
  if (wells_.find(well_name) == wells_.end()) {
    string em = "Well " + well_name + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wwirXd_[well_name];
}

const VectorXd &ECLSummaryReader::wgirXd(string well_name) {
  if (wells_.find(well_name) == wells_.end()) {
    string em = "Well " + well_name + " was not found in the summary.";
    throw SmryVarDoesNotExistExc(em);
  }
  return wgirXd_[well_name];
}

// Segment vectors RATES (VectorXd) ------------------------

const vector<VectorXd > &ECLSummaryReader::seg_sofrXd(string wname) {
  for (int ii=0; ii < seg_sofrXd_[wname].size(); ++ii) {
    if (seg_sofrXd_[wname][ii].size() == 0) {
      warnPropertyZero("SOFRXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_sofrXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_swfrXd(string wname) {
  for (int ii=0; ii < seg_swfrXd_[wname].size(); ++ii) {
    if (seg_swfrXd_[wname][ii].size() == 0) {
      warnPropertyZero("SWFRXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_swfrXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_sgfrXd(string wname) {
  for (int ii=0; ii < seg_sgfrXd_[wname].size(); ++ii) {
    if (seg_sgfrXd_[wname][ii].size() == 0) {
      warnPropertyZero("SGFRXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_sgfrXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_slfrXd(string wname) {
  for (int ii=0; ii < seg_slfrXd_[wname].size(); ++ii) {
    if (seg_slfrXd_[wname][ii].size() == 0) {
      warnPropertyZero("SLFRXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_slfrXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_sprpXd(string wname) {
  for (int ii=0; ii < seg_sprpXd_[wname].size(); ++ii) {
    if (seg_sprpXd_[wname][ii].size() == 0) {
      warnPropertyZero("SPRPXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_sprpXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_sprdXd(string wname) {
  for (int ii=0; ii < seg_sprdXd_[wname].size(); ++ii) {
    if (seg_sprdXd_[wname][ii].size() == 0) {
      warnPropertyZero("SPRDXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_sprdXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_swctXd(string wname) {
  for (int ii=0; ii < seg_swctXd_[wname].size(); ++ii) {
    if (seg_swctXd_[wname][ii].size() == 0) {
      warnPropertyZero("SWCTXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_swctXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_scsaXd(string wname) {
  for (int ii=0; ii < seg_scsaXd_[wname].size(); ++ii) {
    if (seg_scsaXd_[wname][ii].size() == 0) {
      warnPropertyZero("SCSAXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_scsaXd_[wname];
}

// Segment vectors TOTALS (VectorXd) -----------------------

const vector<VectorXd > &ECLSummaryReader::seg_softXd(string wname) {
  for (int ii=0; ii < seg_softXd_[wname].size(); ++ii) {
    if (seg_softXd_[wname][ii].size() == 0) {
      warnPropertyZero("SOFTXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_softXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_swftXd(string wname) {
  for (int ii=0; ii < seg_swftXd_[wname].size(); ++ii) {
    if (seg_swftXd_[wname][ii].size() == 0) {
      warnPropertyZero("SWFTXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_swftXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_sgftXd(string wname) {
  for (int ii=0; ii < seg_sgftXd_[wname].size(); ++ii) {
    if (seg_sgftXd_[wname][ii].size() == 0) {
      warnPropertyZero("SGFTXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_sgftXd_[wname];
}

const vector<VectorXd > &ECLSummaryReader::seg_slftXd(string wname) {
  for (int ii=0; ii < seg_slftXd_[wname].size(); ++ii) {
    if (seg_slftXd_[wname][ii].size() == 0) {
      warnPropertyZero("SLFTXd:" + wname + ":" + num2str(ii));
    }
  }
  return seg_slftXd_[wname];
}

// ---------------------------------------------------------

void ECLSummaryReader::warnPropertyZero(const string& wname, string propname) const {
  if (vp_.vSIM >= 3) {
    string wm = "Returning cumulative vector with final value 0.0 for ";
    wm += propname + " for well " + wname + ".";
    ext_warn(wm, md_, cl_, vp_.lnw);
  }
}

void ECLSummaryReader::warnPropertyNotFound(const string& propname) const {
  if (vp_.vSIM >= 3) {
    string wm = "Property " + propname + " was not found in the summary. ";
    wm += "Calculating from corresponding well properties.";
    ext_warn(wm, md_, cl_, vp_.lnw);
  }
}

void ECLSummaryReader::warnPropertyZero(const string& propname) const {
  if (vp_.vSIM >= 3) {
    string wm = "Returning cumulative vector with final value 0.0 for " + propname + ".";
    ext_warn(wm, md_, cl_, vp_.lnw);
  }
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

// Well segment rates --------------------------------------
void ECLSummaryReader::initWellSegRates() {
  const ecl_smspec_type *smspec = ecl_sum_get_smspec(ecl_sum_);
  auto s_obj = settings_->optimizer()->objective();
  stringstream ss;

  if (s_obj.type == Settings::Optimizer::ObjectiveType::Augmented) {

    // for each objf term, treat only terms with segment sections
    for (int ii = 0; ii < s_obj.terms.size(); ++ii) {
      if (!s_obj.terms[ii].segments.empty()) {

        // init segment data structures
        string seg_nm, pn, ln = "";
        vector<double> rseg = vector<double>(time_.size(), 0.0);
        vector<double> tseg = vector<double>(time_.size(), 0.0);
        pn = "s_obj.terms[ii=" + num2str(ii,0) + "].prop_name:_" + s_obj.terms[ii].prop_name;
        Printer::pad_text(pn, 42, '_');
        Printer::pad_text(ln, 138, '-');

        // -------------------------------------------------
        for (const auto& wname : wells_) {

          // nseg = # of segments from settings
          vector<int> nseg = s_obj.terms[ii].segments[wname];
          vector<vector<double>> wseg = vector<vector<double>>(nseg.size(), rseg);

          // loop through each specified segment
          for (int jj=0; jj < nseg.size(); ++jj) {

            // SOFR ----------------------------------------
            seg_nm = "SOFR" + wname + num2str(nseg[jj], 0);
            assert(ecl_sum_has_general_var(ecl_sum_ , seg_nm.c_str()));

            // get index for SOFR
            int sofr_index = ecl_smspec_get_general_var_params_index(smspec, seg_nm.c_str());
            double_vector_type *sofr = ecl_sum_alloc_data_vector(ecl_sum_, sofr_index, true);
            assert(double_vector_size(sofr) == time_.size());

            // loop through time at current segment
            for (int kk = 0; kk < time_.size(); ++kk) {
              rseg[kk] = double_vector_safe_iget(sofr, kk);
            }
            rseg[0] = 0.0;
            tseg[0] = rseg[0] * step_[0];
            for (int kk = 1; kk < time_.size(); ++kk) {
              tseg[kk] = tseg[kk-1] + rseg[kk] * step_[kk];
            }
            seg_sofr_[wname].push_back(rseg);
            seg_soft_[wname].push_back(tseg);
            double_vector_free(sofr);

            if (vp_.vSIM >= 5) {
              ss << pn << "seg_nm:_" << seg_nm << "___nseg[" << jj << "]:_" << nseg[jj]
                 << "___seg_sofr_[" << wname << "].sz():_" << seg_sofr_[wname].size()
                 << "___seg_sofr_[" << wname << "].[jj].sz():_" << seg_sofr_[wname][jj].size() << " ";
              // ext_info(ss.str(), md_, cl_, vp_.lnw); ss.str("");
            }

            // SWFR ----------------------------------------
            seg_nm = "SWFR" + wname + num2str(nseg[jj], 0);
            assert(ecl_sum_has_general_var(ecl_sum_ , seg_nm.c_str()));

            // get index for SWFR
            int swfr_index = ecl_smspec_get_general_var_params_index(smspec, seg_nm.c_str());
            double_vector_type *swfr = ecl_sum_alloc_data_vector(ecl_sum_, swfr_index, true);
            assert(double_vector_size(swfr) == time_.size());

            // loop through time at current segment
            for (int kk = 0; kk < time_.size(); ++kk) {
              rseg[kk] = double_vector_safe_iget(swfr, kk);
            }
            rseg[0] = 0.0;
            tseg[0] = rseg[0] * step_[0];
            for (int kk = 1; kk < time_.size(); ++kk) {
              tseg[kk] = tseg[kk-1] + rseg[kk] * step_[kk];
            }
            seg_swfr_[wname].push_back(rseg);
            seg_swft_[wname].push_back(tseg);
            double_vector_free(swfr);

            if (vp_.vSIM >= 5) {
              ss << pn << "seg_nm:_" << seg_nm << "___nseg[" << jj << "]:_" << nseg[jj]
                 << "___seg_swfr_[" << wname << "].sz():_" << seg_swfr_[wname].size()
                 << "___seg_swfr_[" << wname << "].[jj].sz():_" << seg_swfr_[wname][jj].size() << " ";
              // ext_info(ss.str(), md_, cl_, vp_.lnw); ss.str("");
            }

            // SGFR ----------------------------------------
            seg_nm = "SGFR" + wname + num2str(nseg[jj], 0);
            assert(ecl_sum_has_general_var(ecl_sum_ , seg_nm.c_str()));

            // get index for SGFR
            int sgfr_index = ecl_smspec_get_general_var_params_index(smspec, seg_nm.c_str());
            double_vector_type *sgfr = ecl_sum_alloc_data_vector(ecl_sum_, sgfr_index, true);
            assert(double_vector_size(sgfr) == time_.size());

            // loop through time at current segment
            for (int kk = 0; kk < time_.size(); ++kk) {
              rseg[kk] = double_vector_safe_iget(sgfr, kk);
            }
            rseg[0] = 0.0;
            tseg[0] = rseg[0] * step_[0];
            for (int kk = 1; kk < time_.size(); ++kk) {
              tseg[kk] = tseg[kk-1] + rseg[kk] * step_[kk];
            }

            seg_sgfr_[wname].push_back(rseg);
            seg_sgft_[wname].push_back(tseg);
            double_vector_free(sgfr);

            if (vp_.vSIM >= 5) {
              ss << pn << "seg_nm:_" << seg_nm << "___nseg[" << jj << "]:_" << nseg[jj]
                 << "___seg_sgfr_[" << wname << "].sz():_" << seg_sgfr_[wname].size()
                 << "___seg_sgfr_[" << wname << "].[jj].sz():_" << seg_sgfr_[wname][jj].size() << " ";
              // ext_info(ss.str(), md_, cl_, vp_.lnw); ss.str("");
            }

            // SLFR / SLFT ---------------------------------
            for (int kk = 0; kk < time_.size(); ++kk) {
              rseg[kk] = seg_sofr_[wname][jj][kk] + seg_swfr_[wname][jj][kk];
              tseg[kk] = seg_soft_[wname][jj][kk] + seg_swft_[wname][jj][kk];
            }
            seg_slfr_[wname].push_back(rseg);
            seg_slft_[wname].push_back(tseg);

            // SPRP ----------------------------------------
            seg_nm = "SPR" + wname + num2str(nseg[jj], 0);
            assert(ecl_sum_has_general_var(ecl_sum_ , seg_nm.c_str()));

            // get index for SPRP
            int sprp_index = ecl_smspec_get_general_var_params_index(smspec, seg_nm.c_str());
            double_vector_type *sprp = ecl_sum_alloc_data_vector(ecl_sum_, sprp_index, true);
            assert(double_vector_size(sprp) == time_.size());

            // loop through time at current segment
            for (int kk = 0; kk < time_.size(); ++kk) {
              rseg[kk] = double_vector_safe_iget(sprp, kk);
            }
            rseg[0] = 0.0;
            seg_sprp_[wname].push_back(rseg);
            double_vector_free(sprp);

            if (vp_.vSIM >= 5) {
              ss << pn << "seg_nm:_" << seg_nm << "____nseg[" << jj << "]:_" << nseg[jj]
                 << "___seg_sprp_[" << wname << "].sz():_" << seg_sprp_[wname].size()
                 << "___seg_sprp_[" << wname << "].[jj].sz():_" << seg_sprp_[wname][jj].size() << " ";
              // ext_info(ss.str(), md_, cl_, vp_.lnw); ss.str("");
            }

            // SPRD ----------------------------------------
            seg_nm = "SPRD" + wname + num2str(nseg[jj], 0);
            assert(ecl_sum_has_general_var(ecl_sum_ , seg_nm.c_str()));

            // get index for SPRD
            int sprd_index = ecl_smspec_get_general_var_params_index(smspec, seg_nm.c_str());
            double_vector_type *sprd = ecl_sum_alloc_data_vector(ecl_sum_, sprd_index, true);
            assert(double_vector_size(sprd) == time_.size());

            // loop through time at current segment
            for (int kk = 0; kk < time_.size(); ++kk) {
              rseg[kk] = double_vector_safe_iget(sprd, kk);
            }
            rseg[0] = 0.0;
            seg_sprd_[wname].push_back(rseg);
            double_vector_free(sprd);

            if (vp_.vSIM >= 5) {
              ss << pn << "seg_nm:_" << seg_nm << "___nseg[" << jj << "]:_" << nseg[jj]
                 << "___seg_sprd_[" << wname << "].sz():_" << seg_sprd_[wname].size()
                 << "___seg_sprd_[" << wname << "].[jj].sz():_" << seg_sprd_[wname][jj].size() << " ";
              // ext_info(ss.str(), md_, cl_, vp_.lnw); ss.str("");
            }

            // SWCT ----------------------------------------
            seg_nm = "SWCT" + wname + num2str(nseg[jj], 0);
            assert(ecl_sum_has_general_var(ecl_sum_ , seg_nm.c_str()));

            // get index for SWCT
            int swct_index = ecl_smspec_get_general_var_params_index(smspec, seg_nm.c_str());
            double_vector_type *swct = ecl_sum_alloc_data_vector(ecl_sum_, swct_index, true);
            assert(double_vector_size(swct) == time_.size());

            // loop through time at current segment
            for (int kk = 0; kk < time_.size(); ++kk) {
              rseg[kk] = double_vector_safe_iget(swct, kk);
            }
            rseg[0] = 0.0;
            seg_swct_[wname].push_back(rseg);
            double_vector_free(swct);

            if (vp_.vSIM >= 5) {
              ss << pn << "seg_nm:_" << seg_nm << "___nseg[" << jj << "]:_" << nseg[jj]
                 << "___seg_swct_[" << wname << "].sz():_" << seg_swct_[wname].size()
                 << "___seg_swct_[" << wname << "].[jj].sz():_" << seg_swct_[wname][jj].size() << " ";
              // ext_info(ss.str(), md_, cl_, vp_.lnw); ss.str("");
            }

            // SCSA ----------------------------------------
            seg_nm = "SCSA" + wname + num2str(nseg[jj], 0);
            assert(ecl_sum_has_general_var(ecl_sum_ , seg_nm.c_str()));

            // get index for SCSA
            int scsa_index = ecl_smspec_get_general_var_params_index(smspec, seg_nm.c_str());
            double_vector_type *scsa = ecl_sum_alloc_data_vector(ecl_sum_, scsa_index, true);
            assert(double_vector_size(scsa) == time_.size());

            // loop through time at current segment
            for (int kk = 0; kk < time_.size(); ++kk) {
              rseg[kk] = double_vector_safe_iget(scsa, kk);
            }
            rseg[0] = 0.0;
            seg_scsa_[wname].push_back(rseg);
            double_vector_free(scsa);

            if (vp_.vSIM >= 5) {
              ss << pn << "seg_nm:_" << seg_nm << "___nseg[" << jj << "]:_" << nseg[jj]
                 << "___seg_scsa_[" << wname << "].sz():_" << seg_scsa_[wname].size()
                 << "___seg_scsa_[" << wname << "].[jj].sz():_" << seg_scsa_[wname][jj].size() << " ";
              // ext_info(ss.str(), md_, cl_, vp_.lnw); ss.str("");
              ss << ln << " ";
            }

            // dbg -- 1
            if (vp_.vSIM >= 6) {
              cout << "[ECLSummaryReader::initWellSegRates() -- 1] Seg#" << jj << endl;
              for (int kk = 0; kk < time_.size(); ++kk) {
                cout <<   "[j:" << num2str(jj, 0, 0, 2);
                cout << "|i:" << num2str(kk, 0, 0, 3) << "] ";
                cout <<  " -- [ sofr: " << num2str(seg_sofr_[wname][jj][kk],3, 0, 8);
                cout << "] -- [ swfr: " << num2str(seg_swfr_[wname][jj][kk],3, 0, 8);
                cout << "] -- [ sprp: " << num2str(seg_sprp_[wname][jj][kk],3, 0, 8);
                cout << "] -- [ sprd: " << num2str(seg_sprd_[wname][jj][kk],3, 0, 8);
                cout << "] -- [ swct: " << num2str(seg_swct_[wname][jj][kk],3, 0, 8);
                cout << "] -- [ scsa: " << num2str(seg_scsa_[wname][jj][kk],3, 1, 8);
                cout << endl;
              }
            }

          } // for each segment
          // ext_info(ss.str(), md_, cl_, 140, 2); ss.str("");

        } // for each well

        // dbg -- 2
        if (vp_.vSIM >= 6) {
          cout << ln << endl;
          for (const auto &wname : wells_) {
            for (int jj = 0; jj < seg_sofr_[wname].size(); ++jj) {
              cout << "[ECLSummaryReader::initializeVectors() -- 2] Seg#" << jj << endl;
              for (int kk = 0; kk < time_.size(); ++kk) {
                cout <<   "[j:" << num2str(jj, 0, 0, 2);
                cout << "|i:" << num2str(kk, 0, 0, 3) << "] ";
                cout <<  " -- [ sofr: " << num2str(seg_sofr_[wname][jj][kk], 3, 0, 8);
                cout << "] -- [ swfr: " << num2str(seg_swfr_[wname][jj][kk], 3, 0, 8);
                cout << "] -- [ sprp: " << num2str(seg_sprp_[wname][jj][kk], 3, 0, 8);
                cout << "] -- [ sprd: " << num2str(seg_sprd_[wname][jj][kk], 3, 0, 8);
                cout << "] -- [ swct: " << num2str(seg_swct_[wname][jj][kk], 3, 0, 8);
                cout << "] -- [ scsa: " << num2str(seg_scsa_[wname][jj][kk], 3, 1, 8);
                cout << endl;
              }
            }
          }
          cout << ln << endl;
        }

        return;
      } // if segments
    } // for each objf term
  } // if augmented objf type
}


}
}
