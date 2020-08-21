/***********************************************************
Copyright (C) 2016
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020- Mathias Bellout
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

#ifndef RESULTS_H
#define RESULTS_H

#include <QString>
#include "results_exceptions.h"
#include <vector>
#include "Simulation/results/json_results.h"
#include <Eigen/Core>

namespace Simulation {
namespace Results {

/*!
 * \brief Results class handles access to the results of an
 * evaluated model. It acts as an interface to the utilities
 * doing the actual reading of simulator summary/result files.
 */
class Results
{
 public:
  /*!
   * \brief Property enum Lists the properties that can be
   * retrieved from the results.
   *
   * If results other than these are to be retrieved, they
   * must be added here and functionality to retrieve them
   * must be implemented in the underlying classes.
   */
  enum Property {
    FieldOilProdTotal,
    FieldGasProdTotal,
    FieldWatProdTotal,
    FieldWatInjTotal,
    FieldGasInjTotal,
    WellOilProdTotal,
    WellGasProdTotal,
    WellWatProdTotal,
    WellWatInjTotal,
    WellGasInjTotal,
    Time,
    Step,
    FieldOilProdRate,
    FieldGasProdRate,
    FieldWatProdRate,
    FieldWatInjRate,
    FieldGasInjRate,
    WellOilProdRate,
    WellGasProdRate,
    WellWatProdRate,
    WellWatInjRate,
    WellGasInjRate,
    //
    FieldTotal,
    GnGFunction,
    WellWatBrkTmTotal,
    AuxProp
//    EclAjdGradients
  };

  static string GetPropertyKey(Property prop) {
    string pstr = "???";
    if (prop == FieldOilProdTotal)      pstr = "FieldOilProdTotal";
    else if (prop == FieldGasProdTotal) pstr = "FieldGasProdTotal";
    else if (prop == FieldWatProdTotal) pstr = "FieldWatProdTotal";
    else if (prop == FieldWatInjTotal)  pstr = "FieldWatInjTotal";
    else if (prop == WellOilProdTotal)  pstr = "WellOilProdTotal";
    else if (prop == WellGasProdTotal)  pstr = "WellGasProdTotal";
    else if (prop == WellWatProdTotal)  pstr = "WellWatProdTotal";
    else if (prop == Time)              pstr = "Time";
    else if (prop == Step)              pstr = "Step";
    else if (prop == FieldOilProdRate)  pstr = "FieldOilProdRate";
    else if (prop == FieldGasProdRate)  pstr = "FieldGasProdRate";
    else if (prop == FieldWatProdRate)  pstr = "FieldWatProdRate";
    else if (prop == FieldWatInjRate)   pstr = "FieldWatInjRate";
    else if (prop == FieldGasInjRate)   pstr = "FieldGasInjRate";
    else if (prop == WellOilProdRate)   pstr = "WellOilProdRate";
    else if (prop == WellGasProdRate)   pstr = "WellGasProdRate";
    else if (prop == WellWatProdRate)   pstr = "WellWatProdRate";
    else if (prop == WellWatInjRate)    pstr = "WellWatInjRate";
    else if (prop == WellGasInjRate)    pstr = "WellGasInjRate";

    else if (prop == FieldTotal)        pstr = "FieldTotal";
    else if (prop == GnGFunction)       pstr = "GnGFunction";
    else if (prop == WellWatBrkTmTotal) pstr = "WellWatBrkTmTotal";
    else if (prop == AuxProp)           pstr = "AuxProp";
    else throw RsltPropKeyDoesNotExistExc("", pstr);
    return pstr;
  }

  static Property GetPropertyKey(QString prop) {
    if (prop == "FieldOilProdTotal")      return FieldOilProdTotal;
    else if (prop == "FieldGasProdTotal") return FieldGasProdTotal;
    else if (prop == "FieldWatProdTotal") return FieldWatProdTotal;
    else if (prop == "FieldWatInjTotal")  return FieldWatInjTotal;
    else if (prop == "WellOilProdTotal")  return WellOilProdTotal;
    else if (prop == "WellGasProdTotal")  return WellGasProdTotal;
    else if (prop == "WellWatProdTotal")  return WellWatProdTotal;
    else if (prop == "Time")              return Time;
    else if (prop == "Step")              return Step;
    else if (prop == "FieldOilProdRate")  return FieldOilProdRate;
    else if (prop == "FieldGasProdRate")  return FieldGasProdRate;
    else if (prop == "FieldWatProdRate")  return FieldWatProdRate;
    else if (prop == "FieldWatInjRate")   return FieldWatInjRate;
    else if (prop == "FieldGasInjRate")   return FieldGasInjRate;
    else if (prop == "WellOilProdRate")   return WellOilProdRate;
    else if (prop == "WellGasProdRate")   return WellGasProdRate;
    else if (prop == "WellWatProdRate")   return WellWatProdRate;
    else if (prop == "WellWatInjRate")    return WellWatInjRate;
    else if (prop == "WellGasInjRate")    return WellGasInjRate;

    else if (prop == "FieldTotal")        return FieldTotal;
    else if (prop == "GnGFunction")       return GnGFunction;
    else if (prop == "WellWatBrkTmTotal") return WellWatBrkTmTotal;
    else if (prop == "AuxProp")           return AuxProp;
    else throw RsltPropKeyDoesNotExistExc("", prop.toStdString());
  }

  static Property GetPropertyKey(std::string p) {
    auto prop = QString::fromStdString(p);
    return GetPropertyKey(prop);
  }

  struct EclAjdGData {
    std::vector<string> id;
    std::vector<string> domain;
    std::vector<string> type;
    std::vector<double> grad;
    std::vector<double> time;
    std::vector<double> norm;
  };

  /*!
   * \brief ReadResults Read the summary data from file.
   * \param file_path path to summary file without suffixes.
   */
  virtual void ReadResults(QString file_path) = 0;

  /*!
   * \brief DumpResults Dump the results that have been read.
   * Should be called when the results are no longer valid.
   *
   * Should call setUnavailable().
   */
  virtual void DumpResults() = 0;

  /*!
   * \brief GetFinalValue Gets the value of the given _field_
   * or _misc_ property at the final time index.
   * \param prop The property to be retrieved.
   */
  virtual double GetValue(Property prop) = 0;

  /*!
   * \brief GetValueVector Get the vector containing all
   * values for the specified property.
   * \param prop The property to be retrieved.
   */
  virtual std::vector<double> GetValueVector(Property prop) = 0;
  virtual VectorXd GetValueVectorXd(Property prop) = 0;

  /*!
   * \brief GetFinalValue Gets the value of the given property
   * for the given well at the final time index.
   * \param prop The property to be retireved.
   * \param well The well to get the value from.
   */
  virtual double GetValue(Property prop, QString well) = 0;

  /*!
   * \brief GetValue Get the value of the given _field_
   * or _misc_ property at the given time index.
   * \param prop The property to be retrieved.
   * \param time_index The time index to get the value of the property at.
   */
  virtual double GetValue(Property prop, int time_index) = 0;

  /*!
   * \brief GetValue Get the value of the given property for the given well at the
   * given time index.
   * \param prop The property to be retrieved.
   * \param well The well to get the value from.
   * \param time_index The time index to get the property at.
   * \return
   */
  virtual double GetValue(Property prop, QString well,
                          int time_index) = 0;

  /*!
   * \brief available Indicates whether results are available.
   *
   * The constructor sets the underlying value to false. setAvailable() should be called
   * when results are available to indicate it.
   */
  bool isAvailable() const { return available_; }

  JsonResults GetJsonResults() { return json_results_; }
  void SetJsonResults(JsonResults results) { json_results_ = results; }

 protected:
  /*!
   * \brief Results Default constructor. A Results object is
   * not useful on its own; One of the subclasses' constructor
   * must be called.
   */
  explicit Results(Settings::Settings *settings) {
    available_ = false;
    vp_ = settings->global()->verbParams();
  }

  /*!
   * \brief setAvailable Sets the availability to true. This shuld be called when results
   * are available.
   */
  void setAvailable() { available_ = true; }

  /*!
   * \brief setUnavailable Sets the availability to false. Should be called when results
   * are unavailable.
   */
  void setUnavailable() { available_ = false; }

  Settings::VerbParams vp_;
  string md_ = "Simulation";
  string cl_ = "results";

 private:
  bool available_;
  JsonResults json_results_;


};

}}

#endif // RESULTS_H
