/***********************************************************
Copyright (C) 2015-2018
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

#ifndef ECLSUMMARYREADER_H
#define ECLSUMMARYREADER_H

#include <ert/ecl/ecl_sum.h>
#include <ert/ecl/ecl_smspec.h>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <Eigen/Core>

#include "Settings/settings.h"

namespace ERTWrapper {
namespace ECLSummary {

using namespace std;
using namespace Eigen;

/*!
 * \brief ECLSummaryReader class is a wrapper for ecl_sum in
 * ERT. It lets you retrieve information from summary files
 * generated by Eclipse.
 */
class ECLSummaryReader
{
 public:
  /*!
   * \brief ECLSummaryReader Reads the Eclipse summary file
   * \param file_name Path to the eclipse summary, with or
   * without file suffix.
   */
  explicit ECLSummaryReader(const string& file_name,
                            Settings::Settings *settings);
  ~ECLSummaryReader();

  /*!
   * \brief GetMiscVar Get a Misc var. Calls ecl_sum_get_misc_var.
   * \param var_name Name of the variable to get, e.g. TIME or YEARS.
   * \param time_index The time index (0 and up).
   * \return The value of the variable at the specified time index.
   */
  double GetMiscVar(const string& var_name, int time_index);

  /*!
   * \brief GetFieldVar Get a Field variable. Calls ecl_sum_get_field_var.
   * \param var_name Name of the variable to get, e.g. FOPT or FWCT.
   * \param time_index The time index (0 and up).
   * \return The value of the variable at the specified time index.
   */
  double GetFieldVar(const string& var_name, int time_index);

  /*!
   * \brief GetWellVar Get a Well variable. Calls ecl_sum_get_well_var.
   * \param well_name Name of the variable to get, e.g. WOPR or WBHP.
   * \param var_name Name of the well to which the variable belongs, e.g. "PROD" or "INJ" (must be upper-case).
   * \param time_index The time index (0 and up).
   * \return The value of the variable at the specified time index.
   */
  double GetWellVar(const string& well_name,
                    string var_name,
                    int time_index);

  //!< Get the last report step, i.e. the highest possible time index.
  int GetLastReportStep();

  //!< Get the first report step, i.e. the lowest possible time index (usually 0).
  int GetFirstReportStep();

  //!< Check whether the report step is valid, i.e. < last and > first.
  bool HasReportStep(int report_step);

  //!< Get the list of all the keys contained in the summary.
  const set<string> &keys() const { return keys_; }

  //!< Get the list of all wells found in the summary.
  const set<string> &wells() const { return wells_; }

  //!< Get the list of all field-level keys contained in the summary.
  const set<string> &field_keys() const { return field_keys_; }

  //!< Get the list of all well-level keys contained in the summary.
  const set<string> &well_keys() const { return well_keys_; }

  //!< Get the time vector (days).
  const vector<double> &time() const { return time_; }
  const vector<double> &step() const { return step_; }
  const VectorXd &timeXd() const { return timeXd_; }
  const VectorXd &stepXd() const { return stepXd_; }

  const vector<double> &fopt() const;
  const vector<double> &fwpt() const;
  const vector<double> &fgpt() const;
  const vector<double> &fwit() const;
  const vector<double> &fgit() const;

  const vector<double> &fopr() const;
  const vector<double> &fwpr() const;
  const vector<double> &fgpr() const;
  const vector<double> &fwir() const;
  const vector<double> &fgir() const;

  const VectorXd &foptXd() const;
  const VectorXd &fwptXd() const;
  const VectorXd &fgptXd() const;
  const VectorXd &fwitXd() const;
  const VectorXd &fgitXd() const;

  const VectorXd &foprXd() const;
  const VectorXd &fwprXd() const;
  const VectorXd &fgprXd() const;
  const VectorXd &fwirXd() const;
  const VectorXd &fgirXd() const;

  vector<double> wopt(string well_name) const;
  vector<double> wwpt(string well_name) const;
  vector<double> wgpt(string well_name) const;
  vector<double> wwit(string well_name) const;
  vector<double> wgit(const string& well_name) const;

  vector<double> wopr(string well_name) const;
  vector<double> wwpr(string well_name) const;
  vector<double> wgpr(string well_name) const;
  vector<double> wwir(string well_name) const;
  vector<double> wgir(string well_name) const;

 private:
  Settings::Settings* settings_;
  Settings::VerbParams vp_;
  string md_ = "ERTWrapper";

  string cl_ = "ECLSummaryReader";
  string file_name_;

  ecl_sum_type *ecl_sum_;

  set<string> keys_; //!< List of all keys in smry.
  set<string> wells_; //!< List of all wells in smry.
  set<string> field_keys_; //!< List of all field keys in smry.
  set<string> well_keys_; //!< List of all well keys in smry.

  set<string> seg_sofr_keys_;
  set<string> seg_swfr_keys_;
  set<string> seg_spr_keys_;
  set<string> seg_sprd_keys_;
  set<string> seg_swct_keys_;
  set<string> seg_scsa_keys_;
  set<string> comp_keys_;


  //!< Populate key lists using the ecl_sum_select_matching_general_var_list function.
  void populateKeyLists();

  vector<double> time_;
  vector<double> step_;

  vector<double> fopt_;
  vector<double> fwpt_;
  vector<double> fgpt_;
  vector<double> fwit_;
  vector<double> fgit_;

  vector<double> fopr_;
  vector<double> fwpr_;
  vector<double> fgpr_;
  vector<double> fwir_;
  vector<double> fgir_;

  VectorXd timeXd_;
  VectorXd stepXd_;

  VectorXd foptXd_;
  VectorXd fwptXd_;
  VectorXd fgptXd_;
  VectorXd fwitXd_;
  VectorXd fgitXd_;

  VectorXd foprXd_;
  VectorXd fwprXd_;
  VectorXd fgprXd_;
  VectorXd fwirXd_;
  VectorXd fgirXd_;

  map<string, vector<double> > wopt_;
  map<string, vector<double> > wwpt_;
  map<string, vector<double> > wgpt_;
  map<string, vector<double> > wwit_;
  map<string, vector<double> > wgit_;

  map<string, vector<double> > wopr_;
  map<string, vector<double> > wwpr_;
  map<string, vector<double> > wgpr_;
  map<string, vector<double> > wwir_;
  map<string, vector<double> > wgir_;

  struct segSet {
    int nsegs = 0;
    vector<vector<double>> seg_data;
  };

  map<string, segSet> seg_sofr_;
  map<string, segSet> seg_swfr_;
  map<string, segSet> seg_spr_;
  map<string, segSet> seg_sprd_;
  map<string, segSet> seg_swct_;
  map<string, segSet> seg_scsa_;

  void initializeVectors();
  void initVectorsXd();
  void initializeTimeVector();

  void initializeWellRates();
  void initWellTotals();

  void initFieldTotals();
  void initFieldRates();

  void initWellSegRates();

  void warnPropertyZero(const string& wname, string propname) const;
  void warnPropertyNotFound(const string& propname) const;
  void warnPropertyZero(const string& propname) const;

  bool hasWellVar(const string& well_name, string var_name);
  bool hasGroupVar(const string& group_name, string var_name);
  bool hasFieldVar(const string& var_name);
  bool hasBlockVar(int block_nr, const string& var_name);
  bool hasMiscVar(const string& var_name);

  /*!
   * Compute a cumulative vector from a rate vector and time_.
   *
   * The first value of the cumulative vector is set to 0.0,
   * and for the remaining:
   *    cml[i] = (time[i] - time[i-1]) * rate[i-1]
   */
  vector<double> computeCumulativeFromRate(vector<double> rate);
};

}
}

#endif // ECLSUMMARYREADER_H
