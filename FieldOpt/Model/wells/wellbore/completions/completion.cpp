/***********************************************************
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

#include "completion.h"

namespace Model {
namespace Wells {
namespace Wellbore {
namespace Completions {

Completion::Completion(Settings::Model::Well::Completion completion_settings) {
  vp_ = completion_settings.verb_params_;

  switch (completion_settings.type) {
    case Settings::Model::WellCompletionType::Perforation:
      type_ = CompletionType::Perforation;
      break;

    case Settings::Model::WellCompletionType::ICV:
      type_ = CompletionType::ICV;
      break;

    case Settings::Model::WellCompletionType::Packer:
      type_ = CompletionType::Packer;
      break;

    case Settings::Model::WellCompletionType::Annulus:
      type_ = CompletionType::Annulus;
      break;

    case Settings::Model::WellCompletionType::Tubing:
      type_ = CompletionType::Tubing;
      break;
    default:
      throw std::runtime_error("Completion type not recognized.");
  }
}

Completion::Completion(Completion::CompletionType type ) {
  type_ = type;
}

}
}
}
}

