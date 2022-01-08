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

#ifndef CONTINOUS_PROPERTY_H
#define CONTINOUS_PROPERTY_H

#include "property.h"

namespace Model {
namespace Properties {

/*!
 * \brief The ContinuousProperty class describes a continuous
 * property in the model. The value of the property is held
 * as a floating-point number. Continuous properties are
 * typically used for production properties, such as rates
 * and pressures.
 */
class ContinuousProperty : public Property
{
 public:
  explicit ContinuousProperty(double value);
  explicit ContinuousProperty(double value, double grad);

  double value() const { return value_; }
  double valueSc() const { return value_sc_; }
  double grad() const { return grad_; }
  double gradSc() const { return grad_sc_; }

  void setValue(double value);
  void setValueSc(double value);
  void setValue(double value, double grad);
  void setValueSc(double value, double grad);

  void scaleValue();
  void scaleGrad(double &grad);
  void UpdateValue();

  void scaleOtherValue(double &oval);

  void Add(double d); //!< Add d to the value of this property.

  /*!
   * \brief Equals checks whether the value of of this
   * property is equal to the value of another property,
   * optionally within some tolerance.
   * \param other The property to compare this to.
   * \param epsilon Optional tolerance. Default: 0.0
   * \return True if abs(this->value() - other->value()) <= epsilon; otherwise false.
   */
  bool Equals(ContinuousProperty *other, double epsilon=0.0);

  bool EqualsValue(double other_val, double epsilon=1e-6);

  string ToString();

 private:
  double value_;
  double value_sc_;

  double infd_ = std::numeric_limits<double>::infinity();
  double grad_ = -infd_;
  double grad_sc_ = -infd_;

  // double value_bf_scale_ = 0.0;
  // double value_bf_descale_ = 0.0;
};

}
}

#endif // CONTINOUS_PROPERTY_H
