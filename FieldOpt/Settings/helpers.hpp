/***********************************************************
Created by einar on 6/11/19.
Copyright (C) 2019
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

#ifndef SETTINGS_HELPERS
#define SETTINGS_HELPERS

#include <QString>
#include <QList>
#include <QJsonArray>
#include "Utilities/printer.hpp"
#include "Utilities/verbosity.h"

namespace Settings {

using Printer::ext_info;
using Printer::ext_warn;
using Printer::num2str;

/*!
 * @brief Check if container property is to be variable.
 *
 * Done by checking if "IsVariable" field is present,
 * and checking whether its value is true.
 */
inline bool is_prop_variable(const QJsonObject &container) {
  return container.contains("IsVariable")
    && container["IsVariable"].isBool()
    && container["IsVariable"].toBool();
}

/*!
* @brief Sets optional double property.
 * Returns true if property is found, otherwise false.
* @param prop Property that should have its value set.
* @param container JSON object that should contain the property.
* @param prop_name Name of property to find in the container.
*/
inline bool set_opt_prop_double(double &prop,
                                const QJsonObject &container,
                                const QString &prop_name,
                                VerbParams vp) {
  string im, md = "Settings", cl = "Helpers";

  if (container.contains(prop_name) && (container[prop_name].isDouble())) {
    prop = container[prop_name].toDouble();
    return true;

  } else if (container.contains(prop_name) && (container[prop_name].isString())) {
    if (vp.vSET >= 3) { ext_info(prop_name.toStdString() + "given as string.", md, cl); }
    prop = std::stod(container[prop_name].toString().toStdString());
    return true;

  } else {
    if (vp.vSET >= 3) {
      im = "Property " + prop_name.toStdString();
      im += " not found. Using default ("+ num2str(prop) + ").";
      ext_info(im, md, cl);
    }
    return false;
  }
}

/*!
* @brief Set a required double property. Will
 * throw an exception if the property is not found.
* @param prop The property that should have its value set.
* @param container The JSON object that should contain the property.
* @param prop_name The name of the property to find in the container.
*/
inline void set_req_prop_double(double &prop, const QJsonObject &container,
                                const QString &prop_name, VerbParams vp) {
  if (!set_opt_prop_double(prop, container, prop_name, vp)) {
    throw std::runtime_error("Required property " + prop_name.toStdString() + " not found.");
  }
}

/*!
* @brief Set an optional integer property.
 * Will return true if property is found, otherwise false.
* @param prop The property that should have its value set.
* @param container The JSON object that should contain the property.
* @param prop_name The name of the property to find in the container.
*/
inline bool set_opt_prop_int(int &prop, const QJsonObject &container,
                             const QString &prop_name, VerbParams vp) {
  string im, md = "Settings", cl = "Helpers";

  if (container.contains(prop_name) && (container[prop_name].isDouble())) {
    prop = container[prop_name].toInt();
    return true;

  } else if (container.contains(prop_name) && (container[prop_name].isString())) {
    if (vp.vSET >= 3) { ext_info(prop_name.toStdString() + "given as string.", md, cl); }
    prop = std::stoi(container[prop_name].toString().toStdString());
    return true;

  } else {
    if (vp.vSET >= 6) {
      im = "Property " + prop_name.toStdString();
      im += " not found. Using default (" + num2str(prop) + ").";
      ext_info(im, md, cl);
    }
    return false;
  }
}

inline void set_req_prop_int(int &prop, const QJsonObject &container,
                             const QString &prop_name, VerbParams vp) {
  if (!set_opt_prop_int(prop, container, prop_name, vp)) {
    throw std::runtime_error("Required property " + prop_name.toStdString() + " not found.");
  }
}

inline bool set_opt_prop_bool(bool &prop, const QJsonObject &container,
                              const QString &prop_name, VerbParams vp) {
  string im, md = "Settings", cl = "Helpers";

  if (container.contains(prop_name) && container[prop_name].isBool()) {
    prop = container[prop_name].toBool();
    return true;

  } else {
    if (vp.vSET >= 6) {
      im = "Property " + prop_name.toStdString();
      im += " not found. Using default (" + num2str(prop) + ").";
      ext_info(im, md, cl);
    }
    return false;
  }
}

inline void set_req_prop_bool(bool &prop, const QJsonObject &container,
                              const QString &prop_name, VerbParams vp) {
  if (!set_opt_prop_bool(prop, container, prop_name, vp)) {
    throw std::runtime_error("Required property " + prop_name.toStdString() + " not found.");
  }
}

/*!
* @brief Set an optional string property.
 * Will return true if property is found, otherwise false.
* @param prop The property that should have its value set.
* @param container The JSON object that should contain the property.
* @param prop_name The name of the property to find in the container.
*/
inline bool set_opt_prop_string(std::string &prop, const QJsonObject &container,
                                const QString &prop_name, VerbParams vp){
  if (container.contains(prop_name) && container[prop_name].isString()) {
    prop = container[prop_name].toString().toStdString();
    return true;
  } else {
    return false;
  }
}
inline void set_req_prop_string(std::string &prop, const QJsonObject &container,
                                const QString &prop_name, VerbParams vp) {
  if (!set_opt_prop_string(prop, container, prop_name, vp)) {
    throw std::runtime_error("Required property " + prop_name.toStdString() + " not found.");
  }
}

/*!
* @brief Set an optional string array property.
 * Will return true if property is found, otherwise false.
* @param prop The property that should have its value set.
* @param container The JSON object that should contain the property.
* @param prop_name The name of the property to find in the container.
*/
inline bool set_opt_prop_string_array(std::vector<std::string> &prop, const QJsonObject &container,
                                      const QString &prop_name, VerbParams vp){
  if (container.contains(prop_name) && container[prop_name].isArray() && container[prop_name].toArray()[0].isString()) {
    for (auto elem : container[prop_name].toArray()) {
      prop.push_back(elem.toString().toStdString());
    }
    return true;
  } else {
    return false;
  }
}

inline void set_req_prop_string_array(std::vector<std::string> &prop, const QJsonObject &container,
                                      const QString &prop_name, VerbParams vp){
  if (!set_opt_prop_string_array(prop, container, prop_name, vp)) {
    throw std::runtime_error("Required property " + prop_name.toStdString() + " not found.");
  }
}

/*!
* @brief Set an optional int array property.
 * Will return true if property is found, otherwise false.
* @param prop The property that should have its value set.
* @param container The JSON object that should contain the property.
* @param prop_name The name of the property to find in the container.
*/
inline bool set_opt_prop_int_array(std::vector<int> &prop, const QJsonObject &container,
                                   const QString &prop_name, VerbParams vp) {
  if (container.contains(prop_name) && container[prop_name].isArray() && container[prop_name].toArray()[0].isDouble()) {
    for (auto elem : container[prop_name].toArray()) {
      prop.push_back(elem.toInt());
    }
    return true;
  } else {
    return false;
  }
}

inline void set_req_prop_int_array(std::vector<int> &prop, const QJsonObject &container,
  const QString &prop_name, VerbParams vp) {
  if (!set_opt_prop_int_array(prop, container, prop_name, vp)) {
    throw std::runtime_error("Required property " + prop_name.toStdString() + " not found.");
  }
}

/*!
* @brief Set an optional double array property. Will return true if property is found, otherwise false.
* @param prop The property that should have its value set.
* @param container The JSON object that should contain the property.
* @param prop_name The name of the property to find in the container.
*/
inline bool set_opt_prop_double_array(std::vector<double> &prop, const QJsonObject &container, const QString &prop_name) {
  if (container.contains(prop_name) && container[prop_name].isArray() && container[prop_name].toArray()[0].isDouble()) {
    for (auto elem : container[prop_name].toArray()) {
      prop.push_back(elem.toDouble());
    }
    return true;
  }
  else {
    return false;
  }
}

inline void set_req_prop_double_array(std::vector<double> &prop, const QJsonObject &container, const QString &prop_name) {
  if (!set_opt_prop_double_array(prop, container, prop_name)) {
    throw std::runtime_error("Required property " + prop_name.toStdString() + " not found.");
  }
}


}

#endif // SETTINGS_HELPERS
