/***********************************************************
Created: 22.09.2015 2015 by einar

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


#ifndef PROPERTY_H
#define PROPERTY_H

#include "property_exceptions.h"
#include <QString>
#include <QUuid>

#include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>

namespace Model {
namespace Properties {

using Printer::info;
using Printer::idbg;
using Printer::ext_warn;
using Printer::ext_info;
using Printer::num2str;
using Printer::num2strQ;

/*!
 * \brief The Property class is an abstract class implemented
 * by specific property types, i.e. integer, real and binary.
 * It holds and to some extends describes the value of a
 * property in the model.
 */
class Property
{
 public:
  //!< Underlying datatype of the property's value.
  enum Type { Discrete, Continuous, Binary };

  //!< Get the underlying datatype of the property's value.
  Type type() const { return type_; }

  //!< Get the name set for the variable. Returns
  //!< an empty string if name has not been set.
  QString name() const { return name_;}

  //!< Set the name of the variable.
  void setName(QString name) { name_ = name; }

  //!< Check if the property is locked.
  bool IsLocked() const { return locked_; }

  bool checkLock(bool locked) {
    string m =  "Can not change locked real variable.";
    if (locked) { throw PropertyLockedException(m); }
    return true;
  }

  //!< Set the property to locked. The value
  //!< of a locked property cannot be changed.
  void Lock() { locked_ = true; }

  //!< Unlock a property.
  void Unlock() { locked_ = false; }

  /*!
   * \brief //!< Set the property to variable.
   *
   * This will also populate the info data structure for the
   * property. The property name is used to create the Info
   * datastructure, thus the name of the property must be set
   * before setting it to variable.
   */
  void SetVariable();

  //!< Check if the property is variable.
  bool isVariable() { return is_variable_; }

  //!< Unique ID for the property.
  QUuid id() const { return id_; }

  void setBounds(double min, double max);

  /*!
   * @brief Update the UUID for this object. This is used
   * when synchronizing models between process instances
   * when using MPI.
   * @param new_id The new ID for this object.
   */
  void UpdateId(QUuid new_id);

  /*!
   * \brief The type of property represented. This type
   * decides which field in the Info datastructure is set.
   */
  enum PropertyType : int { PropertyTypeNone=2000,
    BHP=2001, Rate=2002, SplinePoint=2003, WellBlock=2004,
    Transmissibility=2005, PseudoContVert=2006, Packer=2007,
    ICD=2008, PolarSpline=2009
  };

  /*!
   * \brief For SplinePoint type properties, this is used to
   * indicate whether the property is for the heel or the toe
   * of the well.
   */
  enum SplineEnd : int { SplineEndNone=3000,
    Heel=3001, Toe=3002, Middle=3003 };

  /*!
   * \brief For PolarSpline type properties.
   */
  enum PolarProp : int { PolarPropNone=5000,
    Azimuth=5001, Length=5002, Elevation=5003, Midpoint=5004 };

  /*!
   * \brief For SplinePoint and WellBlock type variables,
   * this indicates which coordinate the property is for.
   */
  enum Coordinate : int { CoordinateNone=4000,
    x=4001, y=4002, z=4003, i=4004, j=4005, k=4006 };

  /*!
   * \brief The PropertyInfo struct contains information
   * about what the property represents. This struct is
   * only populated if the property is set variable.
   */
  struct PropertyInfo {
    //!< Indicates whether the info has been set or not.
    bool is_set_ = false;

    //!< The type of property this is (part) of.
    PropertyType prop_type = PropertyType::PropertyTypeNone;

    //!< The name of the well this property belongs to.
    QString parent_well_name = "";

    //!< The number of the well block this property belongs
    //!< to, the time step of a control, spline point index,
    //!< or packer number
    int index = -1;

    //!< Which end of a spline is held.
    SplineEnd spline_end = SplineEnd::SplineEndNone;

    //!< Which coordinate (x/y/z/i/j/k) is held.
    Coordinate coord = Coordinate::CoordinateNone;

    //!< PolarSpline properties
    PolarProp polar_prop = PolarProp::PolarPropNone;

    string getPropTypeName() {
      if (prop_type == PropertyType::BHP) {
        return "BHP";
      } else if (prop_type == PropertyType::Rate) {
        return "Rate";
      } else if (prop_type == PropertyType::SplinePoint) {
        return "SplinePoint";
      } else if (prop_type == PropertyType::WellBlock) {
        return "WellBlock";
      } else if (prop_type == PropertyType::Transmissibility) {
        return "Transmissibility";
      } else if (prop_type == PropertyType::PseudoContVert) {
        return "PseudoContVert";
      } else if (prop_type == PropertyType::Packer) {
        return "Packer";
      } else if (prop_type == PropertyType::ICD) {
        return "ICD";
      } else if (prop_type == PropertyType::PolarSpline) {
        return "PolarSpline";
      } else {
        string em = "Unable to recognize property type.";
        throw std::runtime_error(em);
      }
    }

    // TODO Obsolete functions?
    // string getSplineEndName() {
    // enum SplineEnd : int { SplineEndNone=3000,
    //   Heel=3001, Toe=3002, Middle=3003 };
    // }

    // string getPolarPropName() {
    // enum PolarProp : int { PolarPropNone=5000,
    //   Azimuth=5001, Length=5002, Elevation=5003, Midpoint=5004 };
    // }

    // string getCoordinateName() {
    // enum Coordinate : int { CoordinateNone=4000,
    //   x=4001, y=4002, z=4003, i=4004, j=4005, k=4006 };
    // }
  };

  //!< Get an Info struc containing information about the
  //!< variable. This is only set if the property is variable.
  PropertyInfo propertyInfo() const;

  Eigen::Vector2d bounds() { return bounds_; }

  void setScalingCoeffs(double D, double c) {
    bnds_b_mns_a_ = D;
    bnds_a_pls_b_ = c;
  };

 protected:
  explicit Property(Type type) {
    type_ = type;
    locked_ = false;
    is_variable_ = false;
    name_ = "";

    bounds_ << -1, 1;
    bnds_b_mns_a_ = 2.0; // D
    bnds_a_pls_b_ = 0.0; // c
  }

  double bnds_b_mns_a_;
  double bnds_a_pls_b_;

 private:
  QUuid id_ = QUuid::createUuid();
  Type type_;
  bool locked_;
  bool is_variable_;
  QString name_;
  PropertyInfo info_;

  Eigen::Vector2d bounds_;

  void set_property_info();
  PropertyType get_prop_type_name(QString prop_name) const;
  QString get_well_name(QString prop_name) const;

  SplineEnd get_spline_end(QString prop_name) const;
  int get_spline_index(QString prop_name) const;

  int get_prop_index(QString prop_name) const;
  Coordinate get_prop_coord(QString prop_name) const;
  PolarProp get_polar_prop(QString prop_name) const;

  string im_ = "", wm_ = "", em_ = "";
  string cl_ = "Property";
  string md_ = "Model::Properties";

  Settings::VerbParams vp_;

};

}
}

#endif // PROPERTY_H
