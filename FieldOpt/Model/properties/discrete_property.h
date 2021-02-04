/***********************************************************
Created: 23.09.2015 2015 by einar

Copyright (C) 2015
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

#ifndef DISCRETE_PROPERTY_H
#define DISCRETE_PROPERTY_H

#include "property.h"

namespace Model {
namespace Properties {

/*!
 * \brief The DiscreteProperty class describes descibes discrete properties in the model,
 * such as locations in a discrete grid or the number of ICDs in the well. Discrete
 * values are held as integer numbers.
 */
class DiscreteProperty : public Property
{
public:
    DiscreteProperty(int value);

    int value() const { return value_; }
    void setValue(int value);

    void Add(int i); //!< Add i to the value of this property.
    bool Equals(const DiscreteProperty *other) const; // Check if the value of another variable equals the value of this variable.

    QString ToString() const;

private:
    int value_;
};

}
}
#endif // DISCRETE_PROPERTY_H
