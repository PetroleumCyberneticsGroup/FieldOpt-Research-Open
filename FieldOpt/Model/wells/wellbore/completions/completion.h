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

#ifndef COMPLETION_H
#define COMPLETION_H

#include "Settings/model.h"

namespace Model {
namespace Wells {
namespace Wellbore {
namespace Completions {

/*!
 * \brief The Completion class is the superclass for all
 * completions, e.g. perforations and ICDs.
 *
 * This class specifies members common to all completions,
 * e.g. a unique numerical ID for the completion. Completion
 * objects add themselves to their parent well block.
 */
class Completion
{
public:

    enum CompletionType : int { Perforation=1001, ICV=1002, Packer=1003, Tubing=1004, Annulus=1005 };
    CompletionType type() const { return type_; } //!< Get the type of completion, e.g. Perforation.

protected:
    Completion(::Settings::Model::Well::Completion completion_settings);
    Completion(CompletionType type);

private:
    CompletionType type_;
};

}
}
}
}

#endif // COMPLETION_H
