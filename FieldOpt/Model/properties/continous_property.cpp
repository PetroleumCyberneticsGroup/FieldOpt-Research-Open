/***********************************************************
Created: 23.09.2015 2015 by einar

Copyright (C) 2015
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

#include "continous_property.h"
#include <cmath>

namespace Model {
namespace Properties {

ContinuousProperty::
ContinuousProperty(double value) : Property(Continuous) {
  value_ = value;
  value_sc_ = value; // scaled values updated after bounds are set
}

double ContinuousProperty::value() const {
  return value_;
}

double ContinuousProperty::valueSc() const {
  return value_sc_;
}

void ContinuousProperty::setValue(double value) {
  if (IsLocked()) {
    throw PropertyLockedException("Can't change locked real variable.");
  } else {
    value_ = value;
    this->scaleValue();
  }
}

void ContinuousProperty::setValueSc(double value_sc) {
  if (IsLocked()) {
    throw PropertyLockedException("Can't change locked real variable.");
  } else {
    value_sc_ = value_sc;
    this->UpdateValue();
  }
}

void ContinuousProperty::scaleValue() {
  // x - > y
  // y = 2*x/bnds_b_mns_a_ - bnds_a_pls_b_/bnds_b_mns_a_
  // value_bf_scale_ = value_;
  value_sc_ = 2*value_/bnds_b_mns_a_ - bnds_a_pls_b_/bnds_b_mns_a_;

}

void ContinuousProperty::UpdateValue() {
  // y - > x
  // x = y * bnds_b_mns_a_/2 + bnds_a_pls_b_/2
  // value_bf_descale_ = value_;
  value_ = value_sc_ * bnds_b_mns_a_/2 + bnds_a_pls_b_/2;
}

void ContinuousProperty::Add(double d) {
  if (IsLocked()) {
    throw PropertyLockedException("Can't add to locked real variable");
  } else {
    // value_ += d; // obsolete, all d's are now scaled
    value_sc_ += d;
    this->UpdateValue();
  }
}

bool ContinuousProperty::EqualsValue(double other_val, double epsilon) {
  return std::abs(this->value() - other_val) <= epsilon;
}

bool ContinuousProperty::Equals(ContinuousProperty *other, double epsilon) {
  // return std::abs(this->value() - other->value()) <= epsilon;
  return std::abs(this->valueSc() - other->valueSc()) <= epsilon;
}

string ContinuousProperty::ToString() {
  string str = "Type: Continuous      ";
  str += "UUID: " + id().toString().toStdString() + "     ";
  str += "Name: " + name().toStdString() + "     ";
  str += "Value: " + num2str(value(), 5, 0) + " ---> ";
  str += "ValueSc: " + num2str(valueSc(), 5, 0) + "    |";

  auto pi = propertyInfo();
  str += "TypeName: " + pi.getPropTypeName() + "  ";
  str += "ParentWellName: " + pi.parent_well_name.toStdString() + "  ";
  str += "Index: " + num2str(pi.index, 0, 0, 2) + "    |";

  // if (pi.prop_type != PropertyType::PropertyTypeNone) {
  //   str += "TypeName: " + pi.getPropTypeName() + "  ";
  // }


  return str;
}



}
}
