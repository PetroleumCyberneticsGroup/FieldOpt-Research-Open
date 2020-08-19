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

ECLResults::ECLResults(Settings::Simulator *settings)
  : Results(settings) {
  smry_reader_ = 0;
}

void ECLResults::ReadResults(QString file_path) {
  if (vp_.vSIM >= 2) {
    ext_info("Attempting to read results from" + file_path.toStdString(), md_, cl_);
  }
  file_path_ = file_path;
  if (smry_reader_ != 0) delete smry_reader_;
  try {
    smry_reader_ = new ERTWrapper::ECLSummary::ECLSummaryReader(file_path.toStdString());
  } catch (ERTWrapper::SummaryFileNotFoundAtPathException) {
    throw ResultFileNotFoundException(file_path.toLatin1().constData());
  }

  setAvailable();
}

void ECLResults::ReadResults(QString file_path, QString build_dir){ }

void ECLResults::DumpResults() {
  if (smry_reader_ != 0) delete smry_reader_;
  smry_reader_ = 0;
  setUnavailable();
}

double ECLResults::GetValue(Results::Property prop) {
  if (!isAvailable()) throw RsltsNotAvailExc();
  return GetValueVector(prop).back();
}

double ECLResults::GetValue(Results::Property prop, int time_index) {
  if (!isAvailable()) throw RsltsNotAvailExc();
  if (time_index < 0 || time_index >= smry_reader_->time().size())
    throw std::runtime_error("The time index " + boost::lexical_cast<std::string>(time_index) + " is outside the range of the summary.");
  return GetValueVector(prop)[time_index];
}

double ECLResults::GetValue(Results::Property prop, QString well) {
  if (!isAvailable()) throw RsltsNotAvailExc();
  return GetValueVector(prop, well).back();
}

double ECLResults::GetValue(Results::Property prop, QString well, int time_index) {
  if (!isAvailable()) throw RsltsNotAvailExc();
  if (time_index < 0 || time_index >= smry_reader_->time().size())
    throw std::runtime_error("The time index " + boost::lexical_cast<std::string>(time_index) + " is outside the range of the summary.");
  return GetValueVector(prop, well)[time_index];
}

std::vector<double> ECLResults::GetValueVector(Results::Property prop) {
  if (!isAvailable()) throw RsltsNotAvailExc();
  switch (prop) {
    case FieldOilProdTotal: return smry_reader_->fopt();
    case FieldGasProdTotal: return smry_reader_->fgpt();
    case FieldWatProdTotal: return smry_reader_->fwpt();
    case FieldWatInjTotal:  return smry_reader_->fwit();
    case FieldGasInjTotal:  return smry_reader_->fgit();
    case Time:              return smry_reader_->time();
    case Step:              return smry_reader_->step();

    case FieldOilProdRate: return smry_reader_->fopr();
    case FieldGasProdRate: return smry_reader_->fgpr();
    case FieldWatProdRate: return smry_reader_->fwpr();
    case FieldWatInjRate:  return smry_reader_->fwir();
    case FieldGasInjRate:  return smry_reader_->fgir();

    default: {
      string em = "[ECLResults] requested property is not a field or misc property.";
      throw std::runtime_error(em);
    }
  }
}

std::vector<double> ECLResults::GetValueVector(Results::Property prop, QString well_name) {
  if (!isAvailable()) throw RsltsNotAvailExc();
  switch (prop) {
    case WellOilProdTotal: return smry_reader_->wopt(well_name.toStdString());
    case WellGasProdTotal: return smry_reader_->wgpt(well_name.toStdString());
    case WellWatProdTotal: return smry_reader_->wwpt(well_name.toStdString());
    case WellWatInjTotal:  return smry_reader_->wwit(well_name.toStdString());
    case WellGasInjTotal:  return smry_reader_->wgit(well_name.toStdString());
    default: {
      throw runtime_error("In ECLResults: The requested property is not a well property.");
    }
  }
}

VectorXd ECLResults::GetValueVectorXd(Results::Property prop) {
  if (!isAvailable()) throw RsltsNotAvailExc();
  switch (prop) {
    case FieldOilProdTotal: return smry_reader_->foptXd();
    case FieldGasProdTotal: return smry_reader_->fgptXd();
    case FieldWatProdTotal: return smry_reader_->fwptXd();
    case FieldWatInjTotal:  return smry_reader_->fwitXd();
    case FieldGasInjTotal:  return smry_reader_->fgitXd();
    case Time:              return smry_reader_->timeXd();
    case Step:              return smry_reader_->stepXd();

    case FieldOilProdRate: return smry_reader_->foprXd();
    case FieldGasProdRate: return smry_reader_->fgprXd();
    case FieldWatProdRate: return smry_reader_->fwprXd();
    case FieldWatInjRate:  return smry_reader_->fwirXd();
    case FieldGasInjRate:  return smry_reader_->fgirXd();

    default: {
      throw std::runtime_error("ECLResults: Requested property is not a field or misc property.");
    }
  }
}

}}
