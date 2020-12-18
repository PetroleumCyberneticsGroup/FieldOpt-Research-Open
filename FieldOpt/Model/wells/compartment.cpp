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

#include "compartment.h"

namespace Model {
namespace Wells {

using Printer::num2str;
using Printer::ext_warn;

Compartment::Compartment() { }

Compartment::Compartment(const double start_md,
                         const double end_md,
                         const double start_tvd,
                         const double end_tvd,
                         const double well_length,
                         const std::vector<Wellbore::WellBlock *> block_range,
                         const std::vector<int> block_idx_range,
                         const Settings::Model::Well &well_settings,
                         Properties::VarPropContainer *variable_container,
                         std::vector<Compartment> &compartments_) {

  Settings::Model::Well::Completion comp_settings;
  block_range_ = block_range;
  block_idx_range_ = block_idx_range;

  comp_settings.verb_params_ = well_settings.verb_params_;

  if (compartments_.empty()) {
    // assert(start_md == 0.0); // -- ???
    if(start_md != 0.0) {
      wm_ = "Well MD larger than zero: " + num2str(start_md, 3, 0, 7);
      ext_warn(wm_, md_, cl_);
    }

    comp_settings.name = "Packer#" + well_settings.name + "#0";
    comp_settings.placement = start_md / well_length;
    comp_settings.true_vertical_depth = start_tvd;
    comp_settings.type = Settings::Model::WellCompletionType::Packer;
    comp_settings.variable_placement = well_settings.seg_compartment_params.variable_placement;
    start_packer = new Wellbore::Completions::Packer(comp_settings, variable_container);
  } else {
    start_packer = compartments_[compartments_.size() - 1].end_packer;
  }

  comp_settings.name = "Packer#" + well_settings.name + "#" + QString::number(compartments_.size() + 1);
  comp_settings.placement = end_md / well_length;
  comp_settings.true_vertical_depth = end_tvd;
  comp_settings.type = Settings::Model::WellCompletionType::Packer;
  comp_settings.variable_placement = well_settings.seg_compartment_params.variable_placement;

  end_packer = new Wellbore::Completions::Packer(comp_settings, variable_container);

  comp_settings.name = "ICD#" + well_settings.name + "#" + QString::number(compartments_.size());
  comp_settings.placement = start_md / well_length;
  comp_settings.true_vertical_depth = start_tvd;
  comp_settings.icd_segment_length = block_range[0]->getLength();
  comp_settings.valve_size = well_settings.seg_compartment_params.valve_size;
  comp_settings.valve_flow_coeff = well_settings.seg_compartment_params.valve_flow_coeff;
  comp_settings.variable_strength = well_settings.seg_compartment_params.variable_strength;
  comp_settings.type = Settings::Model::WellCompletionType::ICV;
  icd = new Wellbore::Completions::ICD(comp_settings, variable_container);

  if (vp_.vMOD >= 1) {
    std::stringstream ss;
    ss << "Created new compartment " << compartments_.size() + 1 << ".| ";
    ss << "Start Packer MD: " << start_packer->md(well_length) << ".| ";
    ss << "End Packer MD: " << end_packer->md(well_length) << ".| ";
    ss << "ICD MD: " << icd->md(well_length) << " ICD length: " << icd->length() << ".|";
    ext_info(ss.str(), md_, cl_);
  }
}

std::vector<Wellbore::WellBlock *> Compartment::GetBlockRange() const {
  return block_range_;
}

std::vector<int> Compartment::GetBlockIdxange() const {
  return block_idx_range_;
}

double Compartment::GetLength(const double &well_length) const {
  return end_packer->md(well_length) - start_packer->md(well_length);
}

double Compartment::GetTVDDifference() const {
  return end_packer->tvd() - start_packer->tvd();
}

void Compartment::Update() {
  icd->setPlacement(start_packer->placement());
  icd->setTvd(start_packer->tvd());
}

}
}
