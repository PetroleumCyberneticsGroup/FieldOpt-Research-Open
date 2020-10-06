/***********************************************************
Copyright (C) 2020-2021
Mathias Bellout <mathias.bellout@gmail.com>

Created:
[Tue Oct 06 2020 09:35:05 week 40 CET+0200] by bellout

Parts from original file: model.cpp
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

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

void Model::parseICVs(QJsonArray &json_icvs, Model::Well &well) {

  for (int i = 0; i < json_icvs.size(); ++i) {
    Well::Completion comp;
    comp.type = WellCompletionType::ICV;
    QJsonObject json_icv = json_icvs[i].toObject();
    comp.is_variable = is_prop_variable(json_icv);

    //bool name_set = set_opt_prop_string(comp.device_name, json_icv, "DeviceName");
    if (set_opt_prop_string(comp.device_name, json_icv, "DeviceName")) {
      comp.device_names.push_back(comp.device_name);
    } else {
      set_req_prop_string_array(comp.device_names, json_icv, "DeviceNames");
    }

    set_req_prop_double(comp.valve_size, json_icv, "ValveSize");
    set_opt_prop_double(comp.min_valve_size, json_icv, "MinValveSize");
    set_opt_prop_double(comp.max_valve_size, json_icv, "MaxValveSize");
    set_opt_prop_double(comp.valve_flow_coeff, json_icv, "FlowCoefficient");

    if (set_opt_prop_int(comp.segment_index, json_icv, "Segment")) {
      comp.segment_indexes.push_back(comp.segment_index);
    } else {
      set_req_prop_int_array(comp.segment_indexes, json_icv, "Segments");
      assert(comp.segment_indexes.size() == comp.device_names.size());
    }

    if (json_icv.contains("TimeStep")) {
      if (!controlTimeIsDeclared(json_icv["TimeStep"].toInt())) {
        string em = "All time steps must be declared in the ControlTimes array.";
        throw std::runtime_error(em);
      }
      comp.time_step = json_icv["TimeStep"].toInt();
    }

    comp.name = "ICD#" + well.name;
    well.completions.push_back(comp);

    if (vp_.vSET > 1) {
      string im = "Added ICV " + comp.name.toStdString() + " to " + well.name.toStdString();
      im += " with valve size " + num2str(comp.valve_size, 5);
      im += " and flow coefficient " + num2str(comp.valve_flow_coeff, 5);
      im += " at segment idx. " + num2str(comp.segment_index);
      ext_info(im, md_, cl_, vp_.lnw);
    }

    if (comp.is_variable && vp_.vSET > 1) {
      string im = "ICV " + comp.name.toStdString();
      im += " set as variable with name " + comp.name.toStdString();
      ext_info(im, md_, cl_, vp_.lnw);
    }
  }
}

void Model::parseICVCompartmentalization(QJsonArray &icv_compartmentalization,
                                         Well& well) {
  assert(well.completions.size() == 1);
  auto device_names = well.completions[0].device_names;

  for (auto comp : icv_compartmentalization) {
    Well::ICVGroup grp;
    grp.type = WellCompletionType::ICV;
    if (well.completions[0].is_variable) {
      grp.is_variable = true;
    }

    set_req_prop_string(grp.icv_group_name, comp.toObject(), "CompName");
    set_req_prop_string_array(grp.icvs, comp.toObject(), "ICVs");
    if (!set_opt_prop_double(grp.valve_size, comp.toObject(), "ValveSize")) {
      grp.valve_size = well.completions[0].valve_size;
    }

    grp.name = "ICD#" + well.name + "#" + QString::fromStdString(grp.icv_group_name);
    grp.min_valve_size = well.completions[0].min_valve_size;
    grp.max_valve_size = well.completions[0].max_valve_size;
    grp.valve_flow_coeff = well.completions[0].valve_flow_coeff;

    // Check that all device names listed for compartmentalization
    // have been listed as ICD device names.
    for (auto name : grp.icvs) {
      for (int i = 0; i < device_names.size(); i++) {
        if (device_names[i] == name) {
          grp.segment_indexes.push_back(well.completions[0].segment_indexes[i]);
          string im = "Added segment nr. " + num2str(grp.segment_indexes.back());
          im += " for ICV " + name + " in group " + grp.icv_group_name;
          ext_info(im, md_, cl_, vp_.lnw);
        }
      }
      if (std::find(std::begin(device_names), std::end(device_names), name) != std::end(device_names)) {
        continue; // Device was found
      }
      else {
        string em = "Unable to find compartment-device " + name + " in ICV name list";
        throw std::runtime_error(em);
      }

      // Check that none of the devices in the new compartment exist in another compartment
      for (auto other_grp : well.icv_compartments) {
        if (std::find(std::begin(other_grp.icvs), std::end(other_grp.icvs), name) != std::end(other_grp.icvs)) {
          string em = "ICV " + name + " has already been added to another compartment.";
          throw std::runtime_error(em);
        }
      }
    }

    assert(grp.segment_indexes.size() == grp.icvs.size());
    well.icv_compartments.push_back(grp);
  }

  // Check that all devices are assigned to a compartment
  for (auto icv : device_names) {
    bool icv_found = false;
    for (auto grp : well.icv_compartments) {
      if (std::find(std::begin(grp.icvs), std::end(grp.icvs), icv) != std::end(grp.icvs)) {
        icv_found = true;
        break;
      }
    }
    if (!icv_found) {
      string em = "ICV " + icv + " has not been assigned to a compartment.";
      throw std::runtime_error(em);
    }
  }
}

// SEGMENTATION --------------------------------------------

void Model::parseSegmentation(const QJsonObject& json_seg,
                              Well &well) {
  parseSegmentTubing(json_seg, well);
  parseSegmentAnnulus(json_seg, well);
  parseSegmentCompartments(json_seg, well);
}

void Model::parseSegmentTubing(const QJsonObject &json_seg,
                               Model::Well &well) const {
  if (vp_.vSET >= 2) {
    ext_info("Parsing Tubing ...", md_, cl_, vp_.lnw);
  }

  if (json_seg.contains("Tubing")) {
    try {
      well.seg_tubing.diameter = json_seg["Tubing"].toObject()["Diameter"].toDouble();
      well.seg_tubing.roughness = json_seg["Tubing"].toObject()["Roughness"].toDouble();
    } catch ( ... ) {
      throw std::runtime_error("For Tubing, both Diameter and Roughness must be defined.");
    }

  } else {
    if (vp_.vSET >= 1) {
      string im = "Tubing keyword not found in Segmentation.";
      im += " Defaulting Diameter to 0.1 and Roughness to 1.52E-5.";
      ext_info(im, md_, cl_, vp_.lnw);
    }
    well.seg_tubing.diameter = 0.1;
    well.seg_tubing.roughness = 1.52E-5;
  }
  well.seg_tubing.cross_sect_area = M_PI / 4.0 * well.seg_tubing.diameter * well.seg_tubing.diameter;
}

void Model::parseSegmentAnnulus(const QJsonObject &json_seg,
                                Model::Well &well) const {
  if (vp_.vSET >= 2) {
    ext_info("Parsing Annulus ...", md_, cl_, vp_.lnw);
  }
  if (json_seg.contains("Annulus")) {
    try {
      well.seg_annulus.diameter = json_seg["Annulus"].toObject()["Diameter"].toDouble();
      well.seg_annulus.roughness = json_seg["Annulus"].toObject()["Roughness"].toDouble();
      well.seg_annulus.cross_sect_area = json_seg["Annulus"].toObject()["CrossSectionArea"].toDouble();
    } catch ( ... ) {
      string em = "For Annulus, both Diameter, CrossSectionArea and Roughness must be defined.";
      throw std::runtime_error(em);
    }

  } else {
    if (vp_.vSET >= 1) {
      string im = "Annulus keyword not found in Segmentation. Defaulting ";
      im += "Diameter to 0.04, Ac to 8.17E-3 and Roughness to 1.52E-5.";
      ext_info(im, md_, cl_, vp_.lnw);
    }
    well.seg_annulus.diameter = 0.04;
    well.seg_annulus.roughness = 1.52E-5;
    well.seg_annulus.cross_sect_area = 8.17E-3;
  }
}

void Model::parseSegmentCompartments(const QJsonObject &json_seg,
                                     Model::Well &well) const {
  if (vp_.vSET >= 2) {
    ext_info("Parsing Compartments ...", md_, cl_, vp_.lnw);
  }

  if (json_seg.contains("Compartments")) {
    auto json_compts = json_seg["Compartments"].toObject();
    try {
      well.seg_n_compartments = json_compts["Count"].toInt();

      if (json_compts.contains("VariablePackers")
        && json_compts["VariablePackers"].toBool() == true) {
        well.seg_compartment_params.variable_placement = true;
      }

      if (json_compts.contains("VariableICDs") &&
        json_compts["VariableICDs"].toBool() == true) {
        well.seg_compartment_params.variable_strength = true;
      }

      if (json_compts.contains("ICDType")) {
        well.seg_compartment_params.type = WellCompletionType::ICV;
      }

      if (json_compts.contains("ICDValveSize")) {
        well.seg_compartment_params.valve_size = json_compts["ICDValveSize"].toDouble();
      } else {
        if (vp_.vSET >= 1) {
          string im = "ICDValveSize keyword not found in Compartments. Defaulting to 7.85E-5.";
          ext_info(im, md_, cl_, vp_.lnw);
        }
        well.seg_compartment_params.valve_size = 7.85E-5;
      }

      if (json_compts.contains("ICDValveFlowCoeff")) {
        well.seg_compartment_params.valve_flow_coeff = json_compts["ICDValveFlowCoeff"].toDouble();
      } else {
        if (vp_.vSET >= 1) {
          string im = "ICDValveFlowCoeff keyword not found in Compartments. Defaulting to 0.50.";
          ext_info(im, md_, cl_, vp_.lnw);
        }
        well.seg_compartment_params.valve_flow_coeff = 0.50;
      }

    } catch (...) {
      string em = "Something went wrong while parsing the Compartments section.";
      throw std::runtime_error(em);
    }

  } else {
    string em = "Compartments keyword must be specified when using the Segmentation keyword.";
    throw std::runtime_error(em);
  }
}