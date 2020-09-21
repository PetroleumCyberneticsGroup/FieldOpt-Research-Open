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

#include "eclresults.h"
#include "ERTWrapper/ertwrapper_exceptions.h"
#include <QVector>
#include <boost/lexical_cast.hpp>
#include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>

namespace Simulation {
namespace Results {

using namespace Eigen;
using Printer::ext_info;
using std::runtime_error;
using std::to_string;

ECLResults::ECLResults(Settings::Settings *settings)
  : Results(settings) {
  smry_reader_ = nullptr;
  settings_ = settings;
  vp_ = settings_->global()->verbParams();
  // settings_->global()->showVerbParams("[ @ECLResults ]");
}

void ECLResults::ReadResults(QString file_path) {
  file_path_ = file_path;
  auto fp = file_path_.toStdString();
  if (vp_.vSIM >= 2) {
    ext_info("[ECLResults::ReadResults] Reading results from " + fp,
      md_, cl_, vp_.lnw);
  }

  if (smry_reader_ != nullptr) delete smry_reader_;
  try {
    smry_reader_ = new ERTWrapper::ECLSummary::ECLSummaryReader(fp, settings_);
    wells_ = smry_reader_->wells();
  } catch (ERTWrapper::SmryFileNotFoundAtPathExc &e) {
    throw RsltFileNotFoundExc(file_path.toLatin1().constData());
  }
  setAvailable();
}

void ECLResults::ReadResults(const QString& file_path,
                             const QString build_dir){ }

void ECLResults::DumpResults() {
  if (smry_reader_ != nullptr) delete smry_reader_;
  smry_reader_ = nullptr;
  setUnavailable();
}

double ECLResults::GetValue(Results::Property prop) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  return GetValueVector(prop).back();
}

double ECLResults::GetValue(Results::Property prop, int time_index) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));

  if (time_index < 0 || time_index >= smry_reader_->time().size()) {
    E("Time index " + to_string(time_index) + " is outside range of smry.");
  }
  return GetValueVector(prop)[time_index];
}

double ECLResults::GetValue(Results::Property prop, QString well) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  return GetValueVector(prop, well.toStdString()).back();
}

double ECLResults::GetValue(Results::Property prop, QString well, int time_index) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  if (time_index < 0 || time_index >= smry_reader_->time().size()) {
    E("Time index " + to_string(time_index) + " is outside range of smry.");
  }
  return GetValueVector(prop, well.toStdString())[time_index];
}

std::vector<double> ECLResults::GetValueVector(Results::Property prop) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  switch (prop) {

    case FieldOilProdTotal: return smry_reader_->fopt();
    case FieldWatProdTotal: return smry_reader_->fwpt();
    case FieldGasProdTotal: return smry_reader_->fgpt();
    case FieldLiqProdTotal: return smry_reader_->flpt();
    case FieldWatInjTotal:  return smry_reader_->fwit();
    case FieldGasInjTotal:  return smry_reader_->fgit();

    case Time:              return smry_reader_->time();
    case Step:              return smry_reader_->step();

    case FieldOilProdRate: return smry_reader_->fopr();
    case FieldWatProdRate: return smry_reader_->fwpr();
    case FieldGasProdRate: return smry_reader_->fgpr();
    case FieldLiqProdRate: return smry_reader_->flpr();
    case FieldWatInjRate:  return smry_reader_->fwir();
    case FieldGasInjRate:  return smry_reader_->fgir();

    default: {
      string p = GetPropertyKey(prop);
      E("[ECLResults] Requested property" + p + " is not a field or misc property.");
    }
  }
}

vector<double> ECLResults::GetValueVector(Results::Property prop,
                                          string wn) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  switch (prop) {
    case WellOilProdTotal: return smry_reader_->wopt(wn);
    case WellWatProdTotal: return smry_reader_->wwpt(wn);
    case WellGasProdTotal: return smry_reader_->wgpt(wn);
    case WellLiqProdTotal: return smry_reader_->wlpt(wn);
    case WellWatInjTotal:  return smry_reader_->wwit(wn);
    case WellGasInjTotal:  return smry_reader_->wgit(wn);

    case WellBHP:         return smry_reader_->wbhp(wn);
    case WellWaterCut:    return smry_reader_->wwct(wn);

    case WellOilProdRate: return smry_reader_->wopr(wn);
    case WellWatProdRate: return smry_reader_->wwpr(wn);
    case WellGasProdRate: return smry_reader_->wgpr(wn);
    case WellLiqProdRate: return smry_reader_->wlpr(wn);
    case WellWatInjRate:  return smry_reader_->wwir(wn);
    case WellGasInjRate:  return smry_reader_->wgir(wn);

    default: {
      string p = GetPropertyKey(prop);
      E("[ECLResults->GetValueVector] Req. prop. " + p + " is not a well property.");
    }
  }
}

VectorXd ECLResults::GetValueVectorXd(Results::Property prop) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  switch (prop) {

    case FieldOilProdTotal: return smry_reader_->foptXd();
    case FieldWatProdTotal: return smry_reader_->fwptXd();
    case FieldGasProdTotal: return smry_reader_->fgptXd();
    case FieldLiqProdTotal: return smry_reader_->flptXd();
    case FieldWatInjTotal:  return smry_reader_->fwitXd();
    case FieldGasInjTotal:  return smry_reader_->fgitXd();

    case Time:              return smry_reader_->timeXd();
    case Step:              return smry_reader_->stepXd();

    case FieldOilProdRate: return smry_reader_->foprXd();
    case FieldWatProdRate: return smry_reader_->fwprXd();
    case FieldGasProdRate: return smry_reader_->fgprXd();
    case FieldLiqProdRate: return smry_reader_->flprXd();
    case FieldWatInjRate:  return smry_reader_->fwirXd();
    case FieldGasInjRate:  return smry_reader_->fgirXd();

    default: {
      string p = GetPropertyKey(prop);
      E("[ECLResults->GetValueVectorXd] Requ. prop " + p + " is not a field or misc property.");
    }
  }
}

VectorXd ECLResults::GetValueVectorXd(Results::Property prop, string wn) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  switch (prop) {

    case WellOilProdTotal: return smry_reader_->woptXd(wn);
    case WellWatProdTotal: return smry_reader_->wwptXd(wn);
    case WellGasProdTotal: return smry_reader_->wgptXd(wn);
    case WellLiqProdTotal: return smry_reader_->wlptXd(wn);
    case WellWatInjTotal:  return smry_reader_->wwitXd(wn);
    case WellGasInjTotal:  return smry_reader_->wgitXd(wn);

    case WellBHP:         return smry_reader_->wbhpXd(wn);
    case WellWaterCut:    return smry_reader_->wwctXd(wn);

    case WellOilProdRate: return smry_reader_->woprXd(wn);
    case WellWatProdRate: return smry_reader_->wwprXd(wn);
    case WellGasProdRate: return smry_reader_->wgprXd(wn);
    case WellLiqProdRate: return smry_reader_->wlprXd(wn);
    case WellWatInjRate:  return smry_reader_->wwirXd(wn);
    case WellGasInjRate:  return smry_reader_->wgirXd(wn);

    default: {
      string p = GetPropertyKey(prop);
      E("[ECLResults->GetValueVectorXd] Req. prop. " + p + " is not a field or misc property.");
    }
  }
}

vector<vector<double>> ECLResults::GetValVectorSeg(Results::Property prop,
                                                   string wn) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  switch (prop) {

    case WellSegOilFlowRate: return smry_reader_->seg_sofr(wn);
    case WellSegWatFlowRate: return smry_reader_->seg_swfr(wn);
    case WellSegGasFlowRate: return smry_reader_->seg_sgfr(wn);
    case WellSegLiqFlowRate: return smry_reader_->seg_slfr(wn);

    case WellSegPress:       return smry_reader_->seg_sprp(wn);
    case WellSegPressDrop:   return smry_reader_->seg_sprd(wn);
    case WellSegWaterCut:    return smry_reader_->seg_swct(wn);
    case WellSegXSecArea:    return smry_reader_->seg_scsa(wn);

    case WellSegOilFlowTotal: return smry_reader_->seg_soft(wn);
    case WellSegWatFlowTotal: return smry_reader_->seg_swft(wn);
    case WellSegGasFlowTotal: return smry_reader_->seg_sgft(wn);
    case WellSegLiqFlowTotal: return smry_reader_->seg_slft(wn);

    default: {
      string p = GetPropertyKey(prop);
      E("[ECLResults->GetValVectorSeg] Req. pro " + p + " is not a well segment property.");
    }
  }
}

vector<VectorXd> ECLResults::GetValVectorSegXd(Results::Property prop,
                                               string wn) {
  if (!isAvailable()) throw RsltsNotAvailExc(GetPropertyKey(prop));
  switch (prop) {

    case WellSegOilFlowRate: return smry_reader_->seg_sofrXd(wn);
    case WellSegWatFlowRate: return smry_reader_->seg_swfrXd(wn);
    case WellSegGasFlowRate: return smry_reader_->seg_sgfrXd(wn);
    case WellSegLiqFlowRate: return smry_reader_->seg_slfrXd(wn);

    case WellSegPress:       return smry_reader_->seg_sprpXd(wn);
    case WellSegPressDrop:   return smry_reader_->seg_sprdXd(wn);
    case WellSegWaterCut:    return smry_reader_->seg_swctXd(wn);
    case WellSegXSecArea:    return smry_reader_->seg_scsaXd(wn);

    case WellSegOilFlowTotal: return smry_reader_->seg_softXd(wn);
    case WellSegWatFlowTotal: return smry_reader_->seg_swftXd(wn);
    case WellSegGasFlowTotal: return smry_reader_->seg_sgftXd(wn);
    case WellSegLiqFlowTotal: return smry_reader_->seg_slftXd(wn);

    default: {
      string p = GetPropertyKey(prop);
      E("[ECLResults] Requested property " + p + " is not a well segment property.");
    }
  }
}

}}
