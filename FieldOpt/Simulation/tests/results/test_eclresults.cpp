/***********************************************************
Copyright (C) 2015-2018
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

#include <gtest/gtest.h>
#include <QString>
#include <qvector.h>
#include "Simulation/results/eclresults.h"
#include "Simulation/tests/results/test_olympr37_data.hpp"
#include "Settings/tests/test_resource_example_file_paths.hpp"
#include "Settings/tests/test_resource_settings.hpp"

#include "Simulation/tests/test_resource_results.h"

using namespace Simulation::Results;
using namespace TestResources::ExampleFilePaths;
using namespace Printer;

namespace {

using ECLProp = Simulation::Results::ECLResults::Property;

class ECLResultsTest : public TestResources::TestResourceResults {

 protected:
  ECLResultsTest() {
    results_ = new ECLResults(settings_full_);
  }

  virtual ~ECLResultsTest() { }

  virtual void SetUp() { }

  virtual void TearDown() {
    results_->DumpResults();
  }

  Results *results_;
};

TEST_F(ECLResultsTest, ReadSummary) {
  EXPECT_THROW(results_->ReadResults("a"),
               RsltFileNotFoundExc);
  auto fp = QString::fromStdString(ecl_base_horzwell);
  EXPECT_NO_THROW(results_->ReadResults(fp));
}

TEST_F(ECLResultsTest, DumpingAndAvailability) {
  // Make results available
  auto fp = QString::fromStdString(ecl_base_horzwell);
  results_->ReadResults(fp);
  EXPECT_TRUE(results_->isAvailable());
  EXPECT_NO_THROW(results_->GetValue(ECLProp::Time));

  // Make results unavailable
  results_->DumpResults();
  EXPECT_FALSE(results_->isAvailable());
  EXPECT_THROW(results_->GetValue(ECLProp::Time),
               RsltsNotAvailExc);

  // Make results available again
  results_->ReadResults(fp);
  EXPECT_TRUE(results_->isAvailable());
  EXPECT_NO_THROW(results_->GetValue(ECLProp::Time));
}

TEST_F(ECLResultsTest, FieldVariables) {
  auto fp = QString::fromStdString(TestResources::ExampleFilePaths::ecl_base_horzwell);
  results_->ReadResults(fp);
  EXPECT_FLOAT_EQ(187866.44, results_->GetValue(ECLProp::FieldOilProdTotal));
  EXPECT_FLOAT_EQ(115870.73, results_->GetValue(ECLProp::FieldOilProdTotal, 10));
  EXPECT_THROW(results_->GetValue(ECLProp::FieldOilProdTotal, 30),
               std::runtime_error);

  std::vector<double> fopt_vec = results_->GetValueVector(ECLProp::FieldOilProdTotal);
  EXPECT_EQ(21, fopt_vec.size());
  EXPECT_FLOAT_EQ(0, fopt_vec.front());
  EXPECT_FLOAT_EQ(187866.44, fopt_vec.back());
}

TEST_F(ECLResultsTest, MiscVariables) {
  auto fp = QString::fromStdString(TestResources::ExampleFilePaths::ecl_base_horzwell);
  results_->ReadResults(fp);
  EXPECT_FLOAT_EQ(200, results_->GetValue(ECLProp::Time));
  EXPECT_FLOAT_EQ(100, results_->GetValue(ECLProp::Time, 10));
  EXPECT_THROW(results_->GetValue(ECLProp::Time, 30),
               std::runtime_error);

  std::vector<double> time_vec = results_->GetValueVector(ECLProp::Time);
  EXPECT_EQ(21, time_vec.size());
  EXPECT_EQ(0, time_vec.front());
  EXPECT_EQ(200, time_vec.back());
}

TEST_F(ECLResultsTest, WellVariables) {
  auto fp = QString::fromStdString(ecl_base_horzwell);
  results_->ReadResults(fp);
  EXPECT_FLOAT_EQ(1116.8876, results_->GetValue(ECLProp::WellWatProdTotal, "PROD"));
  EXPECT_FLOAT_EQ(524.5061, results_->GetValue(ECLProp::WellWatProdTotal, "PROD", 10));
}

TEST_F(ECLResultsTest, ReadFieldWellSegData) {
  string ln;
  Printer::pad_text(ln, 135, '-');
  cout << ln << endl;

  // =======================================================
  auto time = results_olympr37_->GetValueVector(ECLProp::Time);
  auto step = results_olympr37_->GetValueVector(ECLProp::Step);

  auto fopt = results_olympr37_->GetValueVector(ECLProp::FieldOilProdTotal);
  auto fwpt = results_olympr37_->GetValueVector(ECLProp::FieldWatProdTotal);
  auto fgpt = results_olympr37_->GetValueVector(ECLProp::FieldGasProdTotal);
  auto flpt = results_olympr37_->GetValueVector(ECLProp::FieldLiqProdTotal);
  auto fwit = results_olympr37_->GetValueVector(ECLProp::FieldWatInjTotal);
  auto fgit = results_olympr37_->GetValueVector(ECLProp::FieldGasInjTotal);

  auto fopr = results_olympr37_->GetValueVector(ECLProp::FieldOilProdRate);
  auto fwpr = results_olympr37_->GetValueVector(ECLProp::FieldWatProdRate);
  auto fgpr = results_olympr37_->GetValueVector(ECLProp::FieldGasProdRate);
  auto flpr = results_olympr37_->GetValueVector(ECLProp::FieldLiqProdRate);
  auto fwir = results_olympr37_->GetValueVector(ECLProp::FieldWatInjRate);
  auto fgir = results_olympr37_->GetValueVector(ECLProp::FieldGasInjRate);

  auto timeXd = results_olympr37_->GetValueVector(ECLProp::Time);
  auto stepXd = results_olympr37_->GetValueVector(ECLProp::Step);

  auto foptXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldOilProdTotal);
  auto fwptXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldWatProdTotal);
  auto fgptXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldGasProdTotal);
  auto flptXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldLiqProdTotal);
  auto fwitXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldWatInjTotal);
  auto fgitXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldGasInjTotal);

  auto foprXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldOilProdRate);
  auto fwprXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldWatProdRate);
  auto fgprXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldGasProdRate);
  auto flprXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldLiqProdRate);
  auto fwirXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldWatInjRate);
  auto fgirXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldGasInjRate);

  for (int ii = 0; ii < time.size(); ++ii) {
    // cout <<   "fopt[i:" << ii << "]: " << num2str(fopt[ii], 3, 0, 9);
    // cout << "\tfwpt[i:" << ii << "]: " << num2str(fwpt[ii], 3, 0, 9);
    // cout << "\tfgpt[i:" << ii << "]: " << num2str(fgpt[ii], 3, 0, 9);
    // cout << "\tflpt[i:" << ii << "]: " << num2str(flpt[ii], 3, 0, 9);
    // cout << "\tfwit[i:" << ii << "]: " << num2str(fwit[ii], 3, 0, 9);
    // cout << "\tfgit[i:" << ii << "]: " << num2str(fgit[ii], 3, 0, 9);
    // cout << endl;

    EXPECT_FLOAT_EQ(fopt[ii], foptXd(ii));
    EXPECT_FLOAT_EQ(fwpt[ii], fwptXd(ii));
    EXPECT_FLOAT_EQ(fgpt[ii], fgptXd(ii));
    EXPECT_FLOAT_EQ(flpt[ii], flptXd(ii));
    EXPECT_FLOAT_EQ(fwit[ii], fwitXd(ii));
    EXPECT_FLOAT_EQ(fgit[ii], fgitXd(ii));
  }
  cout << ln << endl;

  for (int ii = 0; ii < time.size(); ++ii) {
    // cout <<   "fopr[i:" << ii << "]: " << num2str(fopr[ii], 3, 0, 9);
    // cout << "\tfgpr[i:" << ii << "]: " << num2str(fgpr[ii], 3, 0, 9);
    // cout << "\tfwpr[i:" << ii << "]: " << num2str(fwpr[ii], 3, 0, 9);
    // cout << "\tfwir[i:" << ii << "]: " << num2str(fwir[ii], 3, 0, 9);
    // cout << "\tfgir[i:" << ii << "]: " << num2str(fgir[ii], 3, 0, 9);
    // cout << endl;

    EXPECT_FLOAT_EQ(fopr[ii], foprXd(ii));
    EXPECT_FLOAT_EQ(fwpr[ii], fwprXd(ii));
    EXPECT_FLOAT_EQ(fgpr[ii], fgprXd(ii));
    EXPECT_FLOAT_EQ(flpr[ii], flprXd(ii));
    EXPECT_FLOAT_EQ(fwir[ii], fwirXd(ii));
    EXPECT_FLOAT_EQ(fgir[ii], fgirXd(ii));
  }
  cout << ln << endl;

  // -------------------------------------------------------
  for (const auto &wn : results_olympr37_->getWells()) {
    auto wopt = results_olympr37_->GetValueVector(ECLProp::WellOilProdTotal, wn);
    auto wwpt = results_olympr37_->GetValueVector(ECLProp::WellWatProdTotal, wn);
    auto wgpt = results_olympr37_->GetValueVector(ECLProp::WellGasProdTotal, wn);
    auto wlpt = results_olympr37_->GetValueVector(ECLProp::WellLiqProdTotal, wn);
    auto wwit = results_olympr37_->GetValueVector(ECLProp::WellWatInjTotal, wn);
    auto wgit = results_olympr37_->GetValueVector(ECLProp::WellGasInjTotal, wn);

    auto wopr = results_olympr37_->GetValueVector(ECLProp::WellOilProdRate, wn);
    auto wwpr = results_olympr37_->GetValueVector(ECLProp::WellWatProdRate, wn);
    auto wgpr = results_olympr37_->GetValueVector(ECLProp::WellGasProdRate, wn);
    auto wlpr = results_olympr37_->GetValueVector(ECLProp::WellLiqProdRate, wn);
    auto wwir = results_olympr37_->GetValueVector(ECLProp::WellWatInjRate, wn);
    auto wgir = results_olympr37_->GetValueVector(ECLProp::WellGasInjRate, wn);

    auto woptXd = results_olympr37_->GetValueVectorXd(ECLProp::WellOilProdTotal, wn);
    auto wwptXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWatProdTotal, wn);
    auto wgptXd = results_olympr37_->GetValueVectorXd(ECLProp::WellGasProdTotal, wn);
    auto wlptXd = results_olympr37_->GetValueVectorXd(ECLProp::WellLiqProdTotal, wn);
    auto wwitXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWatInjTotal, wn);
    auto wgitXd = results_olympr37_->GetValueVectorXd(ECLProp::WellGasInjTotal, wn);

    auto woprXd = results_olympr37_->GetValueVectorXd(ECLProp::WellOilProdRate, wn);
    auto wwprXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWatProdRate, wn);
    auto wgprXd = results_olympr37_->GetValueVectorXd(ECLProp::WellGasProdRate, wn);
    auto wlprXd = results_olympr37_->GetValueVectorXd(ECLProp::WellLiqProdRate, wn);
    auto wwirXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWatInjRate, wn);
    auto wgirXd = results_olympr37_->GetValueVectorXd(ECLProp::WellGasInjRate, wn);

    // cout << "Well: " << wn << endl;
    // cout << "wopt.sz(): " << wopt.size() << " -- woptXd.sz(): " << woptXd.size() << endl;
    // cout << "wwpt.sz(): " << wwpt.size() << " -- wwptXd.sz(): " << wwptXd.size() << endl;
    // cout << "wgpt.sz(): " << wgpt.size() << " -- wgptXd.sz(): " << wgptXd.size() << endl;
    // cout << "wlpt.sz(): " << wwpt.size() << " -- wlptXd.sz(): " << wwptXd.size() << endl;
    // cout << "wwit.sz(): " << wwit.size() << " -- wwitXd.sz(): " << wwitXd.size() << endl;
    // cout << "wgit.sz(): " << wgit.size() << " -- wgitXd.sz(): " << wgitXd.size() << endl;
    // cout << ln << endl;

    for (int ii = 0; ii < time.size(); ++ii) {
      // cout <<   "wopt[i:" << ii << "]: " << num2str(wopt[ii], 3, 0, 9);
      // cout << "\twwpt[i:" << ii << "]: " << num2str(wwpt[ii], 3, 0, 9);
      // cout << "\twgpt[i:" << ii << "]: " << num2str(wgpt[ii], 3, 0, 9);
      // cout << "\twlpt[i:" << ii << "]: " << num2str(wwpt[ii], 3, 0, 9);
      // cout << "\twwit[i:" << ii << "]: " << num2str(wwit[ii], 3, 0, 9);
      // cout << "\twgit[i:" << ii << "]: " << num2str(wgit[ii], 3, 0, 9);
      // cout << endl;

      EXPECT_FLOAT_EQ(wopt[ii], woptXd(ii));
      EXPECT_FLOAT_EQ(wwpt[ii], wwptXd(ii));
      EXPECT_FLOAT_EQ(wgpt[ii], wgptXd(ii));
      EXPECT_FLOAT_EQ(wlpt[ii], wlptXd(ii));
      EXPECT_FLOAT_EQ(wwit[ii], wwitXd(ii));
      EXPECT_FLOAT_EQ(wgit[ii], wgitXd(ii));
    }
    cout << ln << endl;

    for (int ii = 0; ii < time.size(); ++ii) {
      // cout <<   "wopr[i:" << ii << "]: " << num2str(wopr[ii], 3, 0, 9);
      // cout << "\twwpr[i:" << ii << "]: " << num2str(wwpr[ii], 3, 0, 9);
      // cout << "\twgpr[i:" << ii << "]: " << num2str(wgpr[ii], 3, 0, 9);
      // cout << "\twlpr[i:" << ii << "]: " << num2str(wlpr[ii], 3, 0, 9);
      // cout << "\twwir[i:" << ii << "]: " << num2str(wwir[ii], 3, 0, 9);
      // cout << "\twgir[i:" << ii << "]: " << num2str(wgir[ii], 3, 0, 9);
      // cout << endl;

      EXPECT_FLOAT_EQ(wopr[ii], woprXd(ii));
      EXPECT_FLOAT_EQ(wwpr[ii], wwprXd(ii));
      EXPECT_FLOAT_EQ(wgpr[ii], wgprXd(ii));
      EXPECT_FLOAT_EQ(wlpr[ii], wlprXd(ii));
      EXPECT_FLOAT_EQ(wwir[ii], wwirXd(ii));
      EXPECT_FLOAT_EQ(wgir[ii], wgirXd(ii));
    }
    cout << ln << endl;
  }

  // -------------------------------------------------------
  // GetValVectorSeg(prop, wname)

  for (const auto& wn : results_olympr37_->getWells()) {
    auto sofr = results_olympr37_->GetValVectorSeg(ECLProp::WellSegOilFlowRate, wn);
    auto swfr = results_olympr37_->GetValVectorSeg(ECLProp::WellSegWatFlowRate, wn);
    auto sprp = results_olympr37_->GetValVectorSeg(ECLProp::WellSegPress, wn);
    auto sprd = results_olympr37_->GetValVectorSeg(ECLProp::WellSegPressDrop, wn);
    auto swct = results_olympr37_->GetValVectorSeg(ECLProp::WellSegWaterCut, wn);
    auto scsa = results_olympr37_->GetValVectorSeg(ECLProp::WellSegXSecArea, wn);

    auto sofrXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegOilFlowRate, wn);
    auto swfrXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegWatFlowRate, wn);
    auto sprpXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegPress, wn);
    auto sprdXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegPressDrop, wn);
    auto swctXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegWaterCut, wn);
    auto scsaXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegXSecArea, wn);

    // cout << "[ECLResultsTest] Well:" << wn << "\t"
    //      <<     "time.sz(): " << time.size()
    //      << " -- sofr.sz(): " << sofr.size()
    //      << " -- swfr.sz(): " << swfr.size()
    //      << " -- sprp.sz(): " << sprp.size()
    //      << " -- sprd.sz(): " << sprd.size()
    //      << " -- swct.sz(): " << swct.size()
    //      << " -- scsa.sz(): " << scsa.size()
    //      << endl;

    for (int jj = 0; jj < sofr.size(); ++jj) {
      // cout << "[ECLResultsTest] Seg#" << jj << endl;
      for (int ii = 0; ii < time.size(); ++ii) {

        // cout <<   "sofr[j:" << jj << "|i:" <<  ii << "]: " << num2str(sofr[jj][ii], 3, 0, 7);
        // cout << "\tswfr[j:" << jj << "|i:" <<  ii << "]: " << num2str(swfr[jj][ii], 3, 0, 7);
        // cout << "\tsprp[j:" << jj << "|i:" <<  ii << "]: " << num2str(sprp[jj][ii], 3, 0, 7);
        // cout << "\tsprd[j:" << jj << "|i:" <<  ii << "]: " << num2str(sprd[jj][ii], 3, 0, 7);
        // cout << "\tswct[j:" << jj << "|i:" <<  ii << "]: " << num2str(swct[jj][ii], 3, 0, 7);
        // cout << "\tscsa[j:" << jj << "|i:" <<  ii << "]: " << num2str(scsa[jj][ii], 3, 0, 7);
        // cout << endl;

        EXPECT_FLOAT_EQ(sofr[jj][ii], sofrXd[jj](ii));
        EXPECT_FLOAT_EQ(swfr[jj][ii], swfrXd[jj](ii));
        EXPECT_FLOAT_EQ(sprp[jj][ii], sprpXd[jj](ii));
        EXPECT_FLOAT_EQ(sprd[jj][ii], sprdXd[jj](ii));
        EXPECT_FLOAT_EQ(swct[jj][ii], swctXd[jj](ii));
        EXPECT_FLOAT_EQ(scsa[jj][ii], scsaXd[jj](ii));
      }
      // cout << endl;
    }
  }
}

TEST_F(ECLResultsTest, MatchRMSDataField) {
  auto time = results_olympr37_->GetValueVector(ECLProp::Time);
  // auto step = results_olympr37_->GetValueVector(ECLProp::Step);

  auto fopt = results_olympr37_->GetValueVector(ECLProp::FieldOilProdTotal);
  auto fwpt = results_olympr37_->GetValueVector(ECLProp::FieldWatProdTotal);
  auto fgpt = results_olympr37_->GetValueVector(ECLProp::FieldGasProdTotal);
  auto flpt = results_olympr37_->GetValueVector(ECLProp::FieldLiqProdTotal);
  auto fwit = results_olympr37_->GetValueVector(ECLProp::FieldWatInjTotal);
  auto fgit = results_olympr37_->GetValueVector(ECLProp::FieldGasInjTotal);

  auto fopr = results_olympr37_->GetValueVector(ECLProp::FieldOilProdRate);
  auto fwpr = results_olympr37_->GetValueVector(ECLProp::FieldWatProdRate);
  auto fgpr = results_olympr37_->GetValueVector(ECLProp::FieldGasProdRate);
  auto flpr = results_olympr37_->GetValueVector(ECLProp::FieldLiqProdRate);
  auto fwir = results_olympr37_->GetValueVector(ECLProp::FieldWatInjRate);
  auto fgir = results_olympr37_->GetValueVector(ECLProp::FieldGasInjRate);

  auto timeXd = results_olympr37_->GetValueVectorXd(ECLProp::Time);
  // auto stepXd = results_olympr37_->GetValueVectorXd(ECLProp::Step);

  auto foptXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldOilProdTotal);
  auto fwptXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldWatProdTotal);
  auto fgptXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldGasProdTotal);
  auto flptXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldLiqProdTotal);
  auto fwitXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldWatInjTotal);
  auto fgitXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldGasInjTotal);

  auto foprXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldOilProdRate);
  auto fwprXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldWatProdRate);
  auto fgprXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldGasProdRate);
  auto flprXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldLiqProdRate);
  auto fwirXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldWatInjRate);
  auto fgirXd = results_olympr37_->GetValueVectorXd(ECLProp::FieldGasInjRate);

  // -------------------------------------------------------
  RMSData rms = RMSData();
  auto data = rms.olympr37_.data;
  string wn = "";
  int sn = 0;

  // Time
  VectorXd RMStimeXd;
  vector<double> RMStime;
  rms.getRMSVectors(ECLProp::Time, wn, sn, RMStimeXd, RMStime, data, time);
  // cout << "RMS: " << RMStimeXd << endl;

  // Step
  VectorXd RMSstepXd;
  vector<double> RMSstep;
  rms.getRMSVectors(ECLProp::Step, wn, sn, RMSstepXd, RMSstep, data, time);

  // FOPT
  VectorXd RMSfoptXd;
  vector<double> RMSfopt;
  rms.getRMSVectors(ECLProp::FieldOilProdTotal, wn, sn, RMSfoptXd, RMSfopt, data, time);
  // cout << "RMSfoptXd: " << RMSfoptXd << endl;

  // FWPT
  VectorXd RMSfwptXd;
  vector<double> RMSfwpt;
  rms.getRMSVectors(ECLProp::FieldWatProdTotal, wn, sn, RMSfwptXd, RMSfwpt, data, time);
  // cout << "RMSfwptXd: " << RMSfwptXd << endl;

  // FGPT
  VectorXd RMSfgptXd;
  vector<double> RMSfgpt;
  rms.getRMSVectors(ECLProp::FieldGasProdTotal, wn, sn, RMSfgptXd, RMSfgpt, data, time);
  // cout << "RMSfgptXd: " << RMSfgptXd << endl;

  // FLPT
  VectorXd RMSflptXd;
  vector<double> RMSflpt;
  rms.getRMSVectors(ECLProp::FieldLiqProdTotal, wn, sn, RMSflptXd, RMSflpt, data, time);
  // cout << "RMSflptXd: " << RMSflptXd << endl;

  // FWIT
  VectorXd RMSfwitXd;
  vector<double> RMSfwit;
  rms.getRMSVectors(ECLProp::FieldWatInjTotal, wn, sn, RMSfwitXd, RMSfwit, data, time);
  // cout << "RMSfwitXd: " << RMSfwitXd << endl;

  // FGIT
  VectorXd RMSfgitXd;
  vector<double> RMSfgit;
  rms.getRMSVectors(ECLProp::FieldGasInjTotal, wn, sn, RMSfgitXd, RMSfgit, data, time);
  // cout << "RMSfgitXd: " << RMSfgitXd << endl;

  for (int ii = 0; ii < time.size(); ++ii) {
    // VectorXd data
    EXPECT_FLOAT_EQ(timeXd(ii), RMStimeXd(ii));
    // EXPECT_FLOAT_EQ(stepXd(ii), RMSstepXd(ii));

    EXPECT_FLOAT_EQ(foptXd(ii), RMSfoptXd(ii));
    EXPECT_FLOAT_EQ(fwptXd(ii), RMSfwptXd(ii));
    EXPECT_FLOAT_EQ(fgptXd(ii), RMSfgptXd(ii));
    EXPECT_FLOAT_EQ(flptXd(ii), RMSflptXd(ii));
    EXPECT_FLOAT_EQ(fwitXd(ii), RMSfwitXd(ii));
    EXPECT_FLOAT_EQ(fgitXd(ii), RMSfgitXd(ii));

    // std data
    EXPECT_FLOAT_EQ(time[ii], RMStime[ii]);
    // EXPECT_FLOAT_EQ(step[ii], RMSstep[ii]);

    EXPECT_FLOAT_EQ(fopt[ii], RMSfopt[ii]);
    EXPECT_FLOAT_EQ(fwpt[ii], RMSfwpt[ii]);
    EXPECT_FLOAT_EQ(fgpt[ii], RMSfgpt[ii]);
    EXPECT_FLOAT_EQ(flpt[ii], RMSflpt[ii]);
    EXPECT_FLOAT_EQ(fwit[ii], RMSfwit[ii]);
    EXPECT_FLOAT_EQ(fgit[ii], RMSfgit[ii]);
  }
}

TEST_F(ECLResultsTest, MatchRMSDataWells) {
  auto time = results_olympr37_->GetValueVector(ECLProp::Time);

  // -------------------------------------------------------
  RMSData rms = RMSData();
  auto data = rms.olympr37_.data;
  int sn = 0;

  // Time
  VectorXd RMStimeXd;
  vector<double> RMStime;
  rms.getRMSVectors(ECLProp::Time, "", sn, RMStimeXd, RMStime, data, time);

  // Step
  VectorXd RMSstepXd;
  vector<double> RMSstep;
  rms.getRMSVectors(ECLProp::Step, "", sn, RMSstepXd, RMSstep, data, time);

  // -------------------------------------------------------
  VectorXd RMSwoptXd;
  vector<double> RMSwopt;
  VectorXd RMSwwptXd;
  vector<double> RMSwwpt;
  VectorXd RMSwgptXd;
  vector<double> RMSwgpt;
  VectorXd RMSwlptXd;
  vector<double> RMSwlpt;
  VectorXd RMSwwitXd;
  vector<double> RMSwwit;
  VectorXd RMSwgitXd;
  vector<double> RMSwgit;

  VectorXd RMSwbhpXd;
  vector<double> RMSwbhp;
  VectorXd RMSwwctXd;
  vector<double> RMSwwct;

  for (const auto &wn : results_olympr37_->getWells()) {

    // auto time = results_olympr37_->GetValueVector(ECLProp::Time);
    // auto step = results_olympr37_->GetValueVector(ECLProp::Step);

    auto wopt = results_olympr37_->GetValueVector(ECLProp::WellOilProdTotal, wn);
    auto wwpt = results_olympr37_->GetValueVector(ECLProp::WellWatProdTotal, wn);
    auto wgpt = results_olympr37_->GetValueVector(ECLProp::WellGasProdTotal, wn);
    auto wlpt = results_olympr37_->GetValueVector(ECLProp::WellLiqProdTotal, wn);
    auto wwit = results_olympr37_->GetValueVector(ECLProp::WellWatInjTotal, wn);
    auto wgit = results_olympr37_->GetValueVector(ECLProp::WellGasInjTotal, wn);

    auto wopr = results_olympr37_->GetValueVector(ECLProp::WellOilProdRate, wn);
    auto wwpr = results_olympr37_->GetValueVector(ECLProp::WellWatProdRate, wn);
    auto wgpr = results_olympr37_->GetValueVector(ECLProp::WellGasProdRate, wn);
    auto wlpr = results_olympr37_->GetValueVector(ECLProp::WellLiqProdRate, wn);
    auto wwir = results_olympr37_->GetValueVector(ECLProp::WellWatInjRate, wn);
    auto wgir = results_olympr37_->GetValueVector(ECLProp::WellGasInjRate, wn);

    auto wbhp = results_olympr37_->GetValueVector(ECLProp::WellBHP, wn);
    auto wwct = results_olympr37_->GetValueVector(ECLProp::WellWaterCut, wn);

    auto timeXd = results_olympr37_->GetValueVector(ECLProp::Time);
    // auto stepXd = results_olympr37_->GetValueVector(ECLProp::Step);

    auto woptXd = results_olympr37_->GetValueVectorXd(ECLProp::WellOilProdTotal, wn);
    auto wwptXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWatProdTotal, wn);
    auto wgptXd = results_olympr37_->GetValueVectorXd(ECLProp::WellGasProdTotal, wn);
    auto wlptXd = results_olympr37_->GetValueVectorXd(ECLProp::WellLiqProdTotal, wn);
    auto wwitXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWatInjTotal, wn);
    auto wgitXd = results_olympr37_->GetValueVectorXd(ECLProp::WellGasInjTotal, wn);

    auto woprXd = results_olympr37_->GetValueVectorXd(ECLProp::WellOilProdRate, wn);
    auto wwprXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWatProdRate, wn);
    auto wgprXd = results_olympr37_->GetValueVectorXd(ECLProp::WellGasProdRate, wn);
    auto wlprXd = results_olympr37_->GetValueVectorXd(ECLProp::WellLiqProdRate, wn);
    auto wwirXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWatInjRate, wn);
    auto wgirXd = results_olympr37_->GetValueVectorXd(ECLProp::WellGasInjRate, wn);

    auto wbhpXd = results_olympr37_->GetValueVectorXd(ECLProp::WellBHP, wn);
    auto wwctXd = results_olympr37_->GetValueVectorXd(ECLProp::WellWaterCut, wn);

    rms.getRMSVectors(ECLProp::WellOilProdTotal, wn, sn, RMSwoptXd, RMSwopt, data, time);
    rms.getRMSVectors(ECLProp::WellWatProdTotal, wn, sn, RMSwwptXd, RMSwwpt, data, time);
    rms.getRMSVectors(ECLProp::WellGasProdTotal, wn, sn, RMSwgptXd, RMSwgpt, data, time);
    rms.getRMSVectors(ECLProp::WellLiqProdTotal, wn, sn, RMSwlptXd, RMSwlpt, data, time);
    rms.getRMSVectors(ECLProp::WellWatInjTotal, wn, sn, RMSwwitXd, RMSwwit, data, time);
    rms.getRMSVectors(ECLProp::WellGasInjTotal, wn, sn, RMSwgitXd, RMSwgit, data, time);

    rms.getRMSVectors(ECLProp::WellBHP, wn, sn, RMSwbhpXd, RMSwbhp, data, time);
    rms.getRMSVectors(ECLProp::WellWaterCut, wn, sn, RMSwwctXd, RMSwwct, data, time);

    for (int ii = 0; ii < time.size(); ++ii) {
      // VectorXd data
      EXPECT_FLOAT_EQ(woptXd(ii), RMSwoptXd(ii));
      EXPECT_FLOAT_EQ(wwptXd(ii), RMSwwptXd(ii));
      EXPECT_FLOAT_EQ(wgptXd(ii), RMSwgptXd(ii));
      EXPECT_FLOAT_EQ(wlptXd(ii), RMSwlptXd(ii));
      EXPECT_FLOAT_EQ(wwitXd(ii), RMSwwitXd(ii));
      EXPECT_FLOAT_EQ(wgitXd(ii), RMSwgitXd(ii));

      EXPECT_FLOAT_EQ(wbhpXd(ii), RMSwbhpXd(ii));
      EXPECT_NEAR(wwctXd(ii), RMSwwctXd(ii), 1E-6);

      // std data
      EXPECT_FLOAT_EQ(wopt[ii], RMSwopt[ii]);
      EXPECT_FLOAT_EQ(wwpt[ii], RMSwwpt[ii]);
      EXPECT_FLOAT_EQ(wgpt[ii], RMSwgpt[ii]);
      EXPECT_FLOAT_EQ(wlpt[ii], RMSwlpt[ii]);
      EXPECT_FLOAT_EQ(wwit[ii], RMSwwit[ii]);
      EXPECT_FLOAT_EQ(wgit[ii], RMSwgit[ii]);

      EXPECT_FLOAT_EQ(wbhp[ii], RMSwbhp[ii]);
      EXPECT_NEAR(wwct[ii], RMSwwct[ii], 1E-6);
    }
  }
}

TEST_F(ECLResultsTest, MatchRMSDataSegs) {
  auto time = results_olympr37_->GetValueVector(ECLProp::Time);

  // -------------------------------------------------------
  RMSData rms = RMSData();
  auto data = rms.olympr37_.data;
  // string wn = "";
  int sn = 0;

  // Time
  VectorXd RMStimeXd;
  vector<double> RMStime;
  rms.getRMSVectors(ECLProp::Time, "", sn, RMStimeXd, RMStime, data, time);

  // Step
  VectorXd RMSstepXd;
  vector<double> RMSstep;
  rms.getRMSVectors(ECLProp::Step, "", sn, RMSstepXd, RMSstep, data, time);

  // -------------------------------------------------------
  VectorXd RMSsofrXd;
  vector<double> RMSsofr;
  VectorXd RMSswfrXd;
  vector<double> RMSswfr;
  VectorXd RMSsprpXd;
  vector<double> RMSsprp;
  VectorXd RMSsprdXd;
  vector<double> RMSsprd;
  VectorXd RMSwwctXd;
  vector<double> RMSwwct;
  VectorXd RMSscsaXd;
  vector<double> RMSscsa;

  VectorXd RMSsgfrXd;
  vector<double> RMSsgfr;

  for (const auto& wn : results_olympr37_->getWells()) {
    auto sofr = results_olympr37_->GetValVectorSeg(ECLProp::WellSegOilFlowRate, wn);
    auto swfr = results_olympr37_->GetValVectorSeg(ECLProp::WellSegWatFlowRate, wn);
    auto sgfr = results_olympr37_->GetValVectorSeg(ECLProp::WellSegGasFlowRate, wn);
    auto slfr = results_olympr37_->GetValVectorSeg(ECLProp::WellSegLiqFlowRate, wn);

    auto sprp = results_olympr37_->GetValVectorSeg(ECLProp::WellSegPress, wn);
    auto sprd = results_olympr37_->GetValVectorSeg(ECLProp::WellSegPressDrop, wn);
    auto swct = results_olympr37_->GetValVectorSeg(ECLProp::WellSegWaterCut, wn);
    auto scsa = results_olympr37_->GetValVectorSeg(ECLProp::WellSegXSecArea, wn);

    auto sofrXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegOilFlowRate, wn);
    auto swfrXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegWatFlowRate, wn);
    auto sgfrXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegGasFlowRate, wn);
    auto slfrXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegLiqFlowRate, wn);

    auto sprpXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegPress, wn);
    auto sprdXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegPressDrop, wn);
    auto swctXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegWaterCut, wn);
    auto scsaXd = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegXSecArea, wn);

    for (int jj = 0; jj < (int)sofr.size(); ++jj) {
      cout << "[ECLResultsTest] Seg#" << jj << " [nsegs=" << sofr.size() << "] " << endl;

      rms.getRMSVectors(ECLProp::WellSegOilFlowRate, wn, jj, RMSsofrXd, RMSsofr, data, time);
      rms.getRMSVectors(ECLProp::WellSegWatFlowRate, wn, jj, RMSswfrXd, RMSswfr, data, time);
      rms.getRMSVectors(ECLProp::WellSegGasFlowRate, wn, jj, RMSsgfrXd, RMSsgfr, data, time);
      rms.getRMSVectors(ECLProp::WellSegPress, wn, jj, RMSsprpXd, RMSsprp, data, time);
      rms.getRMSVectors(ECLProp::WellSegPressDrop, wn, jj, RMSsprdXd, RMSsprd, data, time);
      rms.getRMSVectors(ECLProp::WellSegWaterCut, wn, jj, RMSwwctXd, RMSwwct, data, time);
      rms.getRMSVectors(ECLProp::WellSegXSecArea, wn, jj, RMSscsaXd, RMSscsa, data, time);

      for (int ii = 0; ii < (int)time.size(); ++ii) {

        // cout <<   "sofr[j:" << jj << "|i:" <<  ii << "]: " << num2str(sofr[jj][ii], 3, 0, 7);
        // cout << "\tswfr[j:" << jj << "|i:" <<  ii << "]: " << num2str(swfr[jj][ii], 3, 0, 7);
        // cout << "\tsprp[j:" << jj << "|i:" <<  ii << "]: " << num2str(sprp[jj][ii], 3, 0, 7);
        // cout << "\tsprd[j:" << jj << "|i:" <<  ii << "]: " << num2str(sprd[jj][ii], 3, 0, 7);
        // cout << "\tswct[j:" << jj << "|i:" <<  ii << "]: " << num2str(swct[jj][ii], 3, 0, 7);
        // cout << "\tscsa[j:" << jj << "|i:" <<  ii << "]: " << num2str(scsa[jj][ii], 3, 0, 7);
        // cout << endl;

        EXPECT_FLOAT_EQ(sofr[jj][ii], sofrXd[jj](ii));
        EXPECT_FLOAT_EQ(swfr[jj][ii], swfrXd[jj](ii));
        EXPECT_FLOAT_EQ(sgfr[jj][ii], sgfrXd[jj](ii));
        EXPECT_FLOAT_EQ(sprp[jj][ii], sprpXd[jj](ii));
        EXPECT_FLOAT_EQ(sprd[jj][ii], sprdXd[jj](ii));
        EXPECT_FLOAT_EQ(swct[jj][ii], swctXd[jj](ii));
        EXPECT_FLOAT_EQ(scsa[jj][ii], scsaXd[jj](ii));
      }
      // cout << endl;
    }
  }
}

TEST_F(ECLResultsTest, TestTotalComputation) {
  // TODO Check for deletion 
  // auto fp = QString::fromStdString(olympr37_T02_base_);
  // results_olympr37_->ReadResults(fp);

  string wn = "PRODX2";
  auto time = results_olympr37_->GetValueVectorXd(ECLProp::Time);
  auto step = results_olympr37_->GetValueVectorXd(ECLProp::Step);
  auto wopr = results_olympr37_->GetValueVectorXd(ECLProp::WellOilProdRate, wn);
  auto wopt = results_olympr37_->GetValueVectorXd(ECLProp::WellOilProdTotal, wn);

  vector<double> rseg = vector<double>(time.size(), 0.0);
  vector<double> tseg = vector<double>(time.size(), 0.0);

  // loop through time at current segment
  for (int kk = 0; kk < time.size(); ++kk) {
    rseg[kk] = wopr(kk);
  }
  rseg[0] = 0.0;
  tseg[0] = rseg[0] * step[0];
  for (int kk = 1; kk < time.size(); ++kk) {
    tseg[kk] = tseg[kk-1] + rseg[kk] * step[kk];
  }

  for (int kk = 1; kk < wopt.size(); ++kk) { // skip first component b/c set to 0
    cout << "kk: " << kk << " (wopt(kk)-tseg[kk])/wopt(kk) = ";
    cout << num2str((wopt(kk)-tseg[kk])/wopt(kk),2) << endl;
    if (kk < 11) {
      EXPECT_LT((wopt(kk)-tseg[kk])/wopt(kk), 0.41);
    } else {
      EXPECT_LT((wopt(kk)-tseg[kk])/wopt(kk), 0.10);
    }
  }
}

// TEST_F(ECLResultsTest, TestFieldValue) {
//   auto time = results_olympr37_->GetValueVectorXd(ECLProp::Time);
//   auto step = results_olympr37_->GetValueVectorXd(ECLProp::Step);
//
//   string wn = "PRODX2";
//   auto sofr = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegOilFlowRate, wn);
//   auto swfr = results_olympr37_->GetValVectorSegXd(ECLProp::WellSegWatFlowRate, wn);
//
//   auto wopr = results_olympr37_->GetValueVectorXd(ECLProp::WellOilProdRate, wn);
//   auto wopt = results_olympr37_->GetValueVectorXd(ECLProp::WellOilProdTotal, wn);
//
//   for (int jj = 0; jj < time.size(); ++jj) {
//
//   }
//
// }

}
