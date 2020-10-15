/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
<mathias.bellout@gmail.com>

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

#include "packer.h"

namespace Model {
namespace Wells {
namespace Wellbore {
namespace Completions {

Packer::Packer(const Settings::Model::Well::Completion &completion_settings,
               Properties::VarPropContainer *variable_container) : SegmentedCompletion(completion_settings,
                                                                                       variable_container) {
  if (completion_settings.variable_placement) {
    placement_->setName(completion_settings.name);
    variable_container->AddVariable(placement_);
  }

}

}
}
}
}
