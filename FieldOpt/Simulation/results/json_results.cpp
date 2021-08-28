/***********************************************************
Created by einar on 6/11/19.
Copyright (C) 2015-2019
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020 Mathias Bellout
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


#include "Simulation/results/json_results.h"
#include "Utilities/filehandling.hpp"
#include "QFile"
#include "QByteArray"
#include "QJsonDocument"
#include "QJsonObject"
#include "QJsonArray"
#include "Utilities/printer.hpp"
#include "Utilities/verbosity.h"

namespace Simulation {
namespace Results {

JsonResults::JsonResults(std::string file_path, Settings::Simulator sim_settings) {
  vp_ = sim_settings.verbParams();

  if (vp_.vSIM >= 2) { ext_info("Reading JSON results from " + file_path, md_, cl_); }
  assert(Utilities::FileHandling::FileExists(file_path, vp_));

  if (vp_.vSIM >= 2) { info("JSON results file exists."); }
  QFile *file = new QFile(QString::fromStdString(file_path));
  if (!file->open(QIODevice::ReadOnly)) {
    throw std::runtime_error("Unable to open the JSON results file");
  }
  if (vp_.vSIM >= 2) { info("JSON results file successfully opened."); }
  QByteArray data = file->readAll();
  QJsonDocument json = QJsonDocument::fromJson(data);

  if (json.isNull())
    throw std::runtime_error("Unable to parse the JSON results file.");

  if (!json.isObject())
    throw std::runtime_error("JSON results file format incorrect.");

  QJsonObject json_results = QJsonObject(json.object());
  if (vp_.vSIM >= 2) { info("JSON results file successfully read."); }

  if (sim_settings.read_external_json_results()) {
    if (!json_results.contains("Components") || !json_results["Components"].isArray()) {
      throw std::runtime_error("JSON results file must contain an array named Components.");
    }

    if (vp_.vSIM >= 2) Printer::info("Parsing JSON results file contents.");
    for (auto comp : json_results["Components"].toArray()) {
      QJsonObject component = comp.toObject();
      if (component["Type"] == "Single") {
        singles_[component["Name"].toString().toStdString()] = component["Value"].toDouble();

      } else if (component["Type"] == "Monthly") {
        monthlies_[component["Name"].toString().toStdString()] = std::vector<double>();
        for (auto value : component["Values"].toArray()) {
          monthlies_[component["Name"].toString().toStdString()].push_back(value.toDouble());
        }

      } else if (component["Type"] == "Yearly") {
        yearlies_[component["Name"].toString().toStdString()] = std::vector<double>();
        for (auto value : component["Values"].toArray()) {
          yearlies_[component["Name"].toString().toStdString()].push_back(value.toDouble());
        }

      } else {
        throw std::runtime_error("Parsing only Single, Monthly and Yearly type results.");
      }
    }
  }

  if(sim_settings.read_adj_grad_data()) {
    for (auto adjg : json_results["AdjG"].toArray()) {
      QJsonObject gdata = adjg.toObject();
      grads_map_[gdata["DOMAIN"].toString().toStdString()] = gdata["GRAD"].toDouble();
      grads_vec_.push_back(gdata["GRAD"].toDouble());
      norms_vec_.push_back(gdata["NORM"].toDouble());
      fval_vec_.push_back(gdata["FVAL"].toDouble());
    }
    gradsXd_ = Eigen::Map<VectorXd>(grads_vec_.data(), grads_vec_.size());

    stringstream ss;
    if (vp_.vSIM >= 2 && !grads_vec_.empty()) {
      ss << "Grad data loaded into results object.";
    }
    if (vp_.vSIM >= 4 && !grads_vec_.empty()) {
      ss << "grads_ = { ";
      for (int ii=0; ii < grads_vec_.size(); ii++) {
        ss << "[" << ii << ": " << grads_vec_[ii] << "] ";
      }
      ss << "}; grads_map_ = { ";
      for (auto const& it : grads_map_) {
        ss <<  "[" << it.first << ": " << it.second << "] ";
      }
      ss << "}";
    }
    if (vp_.vSIM >= 2 && !grads_vec_.empty()) {
      ext_info(ss.str(), md_, cl_, vp_.lnw);
    }
  }

  file->close();
  if (vp_.vSIM > 1) {
    ext_info("Done reading additional JSON results.", md_, cl_);
  }
}

double JsonResults::GetSingleValue(std::string name) {
  return singles_[name];
}

std::vector<double> JsonResults::GetMonthlyValues(std::string name) {
  return monthlies_[name];
}

std::vector<double> JsonResults::GetYearlyValues(std::string name) {
  return yearlies_[name];
}

}
}
