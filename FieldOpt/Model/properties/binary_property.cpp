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

#include "binary_property.h"

namespace Model {
namespace Properties {

BinaryProperty::BinaryProperty(bool value)
    : Property(Binary)
{
    value_ = value;
}

void BinaryProperty::setValue(bool value)
{
    if (IsLocked()) throw PropertyLockedException("Cant change value of locked binary variable.");
    else value_ = value;
}

QString BinaryProperty::ToString() const
{
    return QString("Type:\tBinary\nUUID:\t%1\nName:\t%2\nValue:\t%3\n").arg(id().toString()).arg(name()).arg(value());
}


}
}
