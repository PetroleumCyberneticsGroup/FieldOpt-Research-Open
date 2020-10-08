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

Model::Well Model::readSingleWell(QJsonObject json_well) {
  Well well;

  // [WELL]-> name
  well.name = json_well["Name"].toString();
  well.copyVerbParams(vp_);

  // [WELL]-> Type
  QString type = json_well["Type"].toString();

  // [WELL]-> group
  if (json_well.contains("Group")) {
    well.group = json_well["Group"].toString();
  } else {
    well.group = "";
  }

  // [WELL]-> Set well type (Producer or Injector)
  if (QString::compare(type, "Producer") == 0) {
    well.type = WellType::Producer;
  } else if (QString::compare(type, "Injector") == 0) {
    well.type = WellType::Injector;
  } else {
    string em = "Well type " + type.toStdString();
    em += " not recognized for well " + well.name.toStdString();
    throw UnableToParseWellsModelSectionException(em);
  }

  if (vp_.vSET >= 3) {
    stringstream ss;
    ss << "well.name: " + well.name.toStdString() + " | ";
    ss << "well.group: " + well.group.toStdString() + " | ";
    ss << "well.type: " + type.toStdString() + " | ";
    info(ss.str());
  }

  // [WELL]-> set ICVs
  if (json_well.contains("ICVs")) {
    auto json_icvs = json_well["ICVs"].toArray();
    parseICVs(json_icvs, well);

    // Well also contains ICVCompartmentalization
    if (json_well.contains("ICVCompartmentalization")) {
      if (!json_well["ICVCompartmentalization"].isArray()) {
        string im = "Well.ICVCompartmentalization must be an";
        im += "array of objects with the fields: CompName, ICVs.";
        ext_info(im, md_, cl_, vp_.lnw);
        throw std::runtime_error("Unable to parse ICVs.");
      }
      auto json_icv_compartmentalization = json_well["ICVCompartmentalization"].toArray();
      parseICVCompartmentalization(json_icv_compartmentalization, well);
    }
  }

  QString definition_type = json_well["DefinitionType"].toString();
  if (vp_.vSET >= 3) { info("well.def_type: " + definition_type.toStdString()); }

  // Well definition type -> WellBlocks
  if (QString::compare(definition_type, "WellBlocks") == 0) {
    well.definition_type = WellDefinitionType::WellBlocks;
    well.well_blocks = QList<Well::WellBlock>();

    QJsonArray json_well_blocks = json_well["WellBlocks"].toArray();
    for (int i = 0; i < json_well_blocks.size(); ++i) {

      QJsonObject json_well_block_i = json_well_blocks[i].toObject();
      Well::WellBlock block;
      block.i = json_well_block_i["i"].toDouble();
      block.j = json_well_block_i["j"].toDouble();
      block.k = json_well_block_i["k"].toDouble();

      if (json_well_block_i.contains("IsVariable")
        && json_well_block_i["IsVariable"].toBool() == true) {
        block.is_variable = true;
      } else {
        block.is_variable = false;
      }

      if (json_well_block_i.contains("Completion")) {
        QJsonObject json_well_block_i_completion = json_well_block_i["Completion"].toObject();
        block.completion.type = WellCompletionType::Perforation;

        if (json_well_block_i_completion.contains("TransmissibilityFactor")) {
          block.completion.transmissibility_factor = json_well_block_i_completion["TransmissibilityFactor"].toDouble();
        }

        if (json_well_block_i_completion.contains("IsVariable")) {
          block.completion.is_variable = json_well_block_i_completion["IsVariable"].toBool();
        } else {
          block.completion.is_variable = false;
        }

        block.completion.name = "Transmissibility#" + well.name + "#" + QString::number(i);
        block.has_completion = true;

      } else {
        block.has_completion = true;
        block.completion.type = WellCompletionType::Perforation;
        block.completion.transmissibility_factor = -1;
        block.completion.is_variable = false;
      }
      block.name = "WellBlock#" + well.name + "#" + QString::number(i);
      well.well_blocks.append(block);
    } // emd for loop

    // Well definition type -> PolarSpline
  } else if (QString::compare(definition_type, "PolarSpline") == 0) {

    well.definition_type = WellDefinitionType::PolarSpline;
    if (json_well.contains("UseBezierSpline")
      && json_well["UseBezierSpline"].toBool() == true) {
      well.use_bezier_spline = true;
    } else {
      well.use_bezier_spline = false;
    }

    QJsonObject json_pspline = json_well["PolarSpline"].toObject();
    if (!json_pspline.contains("Midpoint")) {
      throw UnableToParseWellsModelSectionException("No midpoint was defined for this spline-type");
    }

    QJsonObject json_midpoint = json_pspline["Midpoint"].toObject();
    Well::PolarSpline polar_spline;
    polar_spline.elevation = json_pspline["Elevation"].toDouble();
    polar_spline.azimuth = json_pspline["Azimuth"].toDouble();
    polar_spline.length = json_pspline["Length"].toDouble();
    polar_spline.midpoint.x = json_midpoint["x"].toDouble();
    polar_spline.midpoint.y = json_midpoint["y"].toDouble();
    polar_spline.midpoint.z = json_midpoint["z"].toDouble();

    if (json_pspline.contains("IsVariable")
      && json_pspline["IsVariable"].toBool() == true){
      polar_spline.is_variable = true;
    }
    well.polar_spline = polar_spline;

    // Well definition type -> WellSpline
  } else if (QString::compare(definition_type, "WellSpline") == 0) {

    well.definition_type = WellDefinitionType::WellSpline;
    if (json_well.contains("UseBezierSpline")
      && json_well["UseBezierSpline"].toBool() == true) {
      well.use_bezier_spline = true;
    } else {
      well.use_bezier_spline = false;
    }

    // SplinePointArray defined
    if (json_well.contains("SplinePointArray")) {
      auto json_sp_array = json_well["SplinePointArray"].toArray();
      assert(json_sp_array.size() >=2);
      for (int i = 0; i < json_sp_array.size(); ++i) {
        auto json_point = json_sp_array[i].toObject();
        Well::SplinePoint point;
        point.x = json_point["x"].toDouble();
        point.y = json_point["y"].toDouble();
        point.z = json_point["z"].toDouble();

        if (json_point.contains("IsVariable")
          && json_point["IsVariable"].toBool() == true) {
          point.is_variable = true;
        } else {
          point.is_variable = false;
        }

        if (i == 0) {
          point.name = "SplinePoint#" + well.name + "#heel";
        } else if (i == json_sp_array.size() - 1) {
          point.name = "SplinePoint#" + well.name + "#toe";
        } else {
          point.name = "SplinePoint#" + well.name + "#P" + QString::number(i);
        }
        well.spline_points.push_back(point);
      }
      well.spline_heel = well.spline_points.front();
      well.spline_toe = well.spline_points.back();

      // SplinePoints defined
    } else {

      QJsonObject json_points = json_well["SplinePoints"].toObject();
      if (!json_points.contains("Heel") || !json_points.contains("Toe")) {
        string em = "Both Heel and Toe must be defined for spline-type wells.";
        throw UnableToParseWellsModelSectionException(em);
      }

      QJsonObject json_heel = json_points["Heel"].toObject();
      QJsonObject json_toe = json_points["Toe"].toObject();
      well.spline_heel.x = json_heel["x"].toDouble();
      well.spline_heel.y = json_heel["y"].toDouble();
      well.spline_heel.z = json_heel["z"].toDouble();
      well.spline_toe.x = json_toe["x"].toDouble();
      well.spline_toe.y = json_toe["y"].toDouble();
      well.spline_toe.z = json_toe["z"].toDouble();

      if (json_heel.contains("IsVariable")
        && json_heel["IsVariable"].toBool()) {
        well.spline_heel.is_variable = true;
      } else {
        well.spline_heel.is_variable = false;
      }

      if (json_toe.contains("IsVariable")
        && json_toe["IsVariable"].toBool()) {
        well.spline_toe.is_variable = true;
      } else {
        well.spline_toe.is_variable = false;
      }

      if ((well.spline_heel.is_variable && well.spline_toe.is_variable)
        || (json_points.contains("IsVariable")
          && json_points["IsVariable"].toBool() == true)) {
        well.is_variable_spline = true;
      }

      well.spline_heel.name = "SplinePoint#" + well.name + "#heel";
      well.spline_toe.name = "SplinePoint#" + well.name + "#toe";

      well.spline_points.push_back(well.spline_heel);
      if (json_points.contains("AdditionalPoints")) {
        int interp_points = json_points["AdditionalPoints"].toInt(); // \todo: <- double?
        for (int p = 0; p < interp_points; ++p) {
          Well::SplinePoint point;
          point.x = well.spline_heel.x + (p+1) * (well.spline_toe.x - well.spline_heel.x) / (interp_points+1);
          point.y = well.spline_heel.y + (p+1) * (well.spline_toe.y - well.spline_heel.y) / (interp_points+1);
          point.z = well.spline_heel.z + (p+1) * (well.spline_toe.z - well.spline_heel.z) / (interp_points+1);
          point.name = "SplinePoint#" + well.name + "#P" + QString::number(p+1);
          if (well.is_variable_spline) {
            point.is_variable = true;
          }
          well.spline_points.push_back(point);
        }
      }
      well.spline_points.push_back(well.spline_toe);
    }

    // Well definition type -> PseudoContVertical2D
  } else if (QString::compare(definition_type, "PseudoContVertical2D") == 0) {

    QJsonObject json_position = json_well["Position"].toObject();
    well.definition_type = WellDefinitionType::PseudoContVertical2D;
    well.pseudo_cont_position.i = json_position["i"].toInt();
    well.pseudo_cont_position.j = json_position["j"].toInt();

    if (json_position.contains("IsVariable")
      && json_position["IsVariable"].toBool() == true) {
      well.pseudo_cont_position.is_variable = true;
    } else {
      well.spline_heel.is_variable = false;
    }

    // Well definition type not recognized
  } else {
    if (vp_.vSET > 1) {
      string im = "Well definition type not recognized.; ";
      im += "Proceeding without defining a well trajectory.";
      ext_warn(im, md_, cl_, vp_.lnw);
    }
    well.definition_type = UNDEFINED;
  }

  // Wellbore radius
  if (json_well.contains("WellboreRadius"))
    well.wellbore_radius = json_well["WellboreRadius"].toDouble();
  else {
    if (vp_.vSET > 1) {
      string im = "WellBoreRadius not set. Defaulting to 0.01905";
      ext_warn(im, md_, cl_, vp_.lnw);
    }
    well.wellbore_radius = 0.1905;
  }

  // Direction of penetration
  if (json_well.contains("Direction")) { // Direction must be specified for horizontal wells
    if (well.definition_type == WellDefinitionType::WellSpline)
      throw std::runtime_error("Direction should not be specified for spline-defined wells");
    if (QString::compare("X", json_well["Direction"].toString()) == 0) well.direction = Direction::X;
    if (QString::compare("Y", json_well["Direction"].toString()) == 0) well.direction = Direction::Y;
    if (QString::compare("Z", json_well["Direction"].toString()) == 0) well.direction = Direction::Z;
  }










  // Controls
  QJsonArray json_controls = json_well["Controls"].toArray();
  well.controls = QList<Well::ControlEntry>();
  for (int i = 0; i < json_controls.size(); ++i) {
    Well::ControlEntry control;

    if (!controlTimeIsDeclared(json_controls.at(i).toObject()["TimeStep"].toInt())) {
      string em = "All time steps must be declared in the ControlTimes array. ";
      em += "Inconsistency detected in Controls declaration.";
      throw UnableToParseWellsModelSectionException(em);
    } else {
      control.time_step = json_controls.at(i).toObject()["TimeStep"].toInt();
    }

    // State (Open or shut)
    if (json_controls[i].toObject().contains("State")
    && QString::compare("Shut", json_controls.at(i).toObject()["State"].toString()) == 0)
      control.state = WellState::WellShut;
    else
      control.state = WellState::WellOpen;


    // Control mode
    QString ctrl_mode = json_controls.at(i).toObject()["Mode"].toString();
    if (ctrl_mode == "BHP") {
      control.control_mode = ControlMode::BHPControl;
      control.name = "BHP#" + well.name + "#" + QString::number(control.time_step);
    }
    else if (ctrl_mode == "Rate" || ctrl_mode == "LRAT" ) {
      control.control_mode = ControlMode::LRATControl;
      control.name = "Rate#" + well.name + "#" + QString::number(control.time_step);
    }
    else if (ctrl_mode == "ORAT" ) {
      control.control_mode = ControlMode::ORATControl;
      control.name = "Rate#" + well.name + "#" + QString::number(control.time_step);
    }
    else if (ctrl_mode == "GRAT" ) {
      control.control_mode = ControlMode::ORATControl;
      control.name = "Rate#" + well.name + "#" + QString::number(control.time_step);
    }
    else if (ctrl_mode == "WRAT" ) {
      control.control_mode = ControlMode::WRATControl;
      control.name = "Rate#" + well.name + "#" + QString::number(control.time_step);
    }
    else if (ctrl_mode == "RESV" ) {
      control.control_mode = ControlMode::RESVControl;
      control.name = "Rate#" + well.name + "#" + QString::number(control.time_step);
    }
    else {
      string em = "Well control type " + json_controls.at(i).toObject()["Mode"].toString().toStdString();
      em += " not recognized for well " + well.name.toStdString();
      throw UnableToParseWellsModelSectionException(em);
    }

    // Control targets/limits
    set_opt_prop_double(control.liq_rate, json_controls[i].toObject(), "Rate");
    set_opt_prop_double(control.liq_rate, json_controls[i].toObject(), "LRAT");
    set_opt_prop_double(control.oil_rate, json_controls[i].toObject(), "ORAT");
    set_opt_prop_double(control.gas_rate, json_controls[i].toObject(), "GRAT");
    set_opt_prop_double(control.wat_rate, json_controls[i].toObject(), "WRAT");
    set_opt_prop_double(control.res_rate, json_controls[i].toObject(), "RESV");
    set_opt_prop_double(control.bhp, json_controls[i].toObject(), "BHP");

    // Injection type
    if (well.type == WellType::Injector) {
      if (!json_controls.at(i).toObject().contains("Type"))
        control.injection_type = InjectionType::WaterInjection;
      if (QString::compare("Water", json_controls.at(i).toObject()["Type"].toString()) == 0)
        control.injection_type = InjectionType::WaterInjection;
      else if (QString::compare("Gas", json_controls.at(i).toObject()["Type"].toString()) == 0)
        control.injection_type = InjectionType::GasInjection;
    }

    if (json_controls[i].toObject()["IsVariable"].toBool()) {
      control.is_variable = true;
    } else {
      control.is_variable = false;
    }
    well.controls.append(control);
  }

  // Preferred Phase
  if (QString::compare("Oil", json_well["PreferredPhase"].toString()) == 0) {
    well.preferred_phase = PreferredPhase::Oil;
  } else if (QString::compare("Water", json_well["PreferredPhase"].toString()) == 0) {
    well.preferred_phase = PreferredPhase::Water;
  } else if (QString::compare("Gas", json_well["PreferredPhase"].toString()) == 0) {
    well.preferred_phase = PreferredPhase::Gas;
  } else if (QString::compare("Liquid", json_well["PreferredPhase"].toString()) == 0) {
    well.preferred_phase = PreferredPhase::Liquid;
  }

  if (json_well.contains("WellSegStruct")) {
    well.wseg_structure = json_well["WellSegStruct"].toString().toStdString();
  }

  // Segmentation
  if (json_well.contains("Segmentation")) {
    well.use_segmented_model = true;
    parseSegmentation(json_well["Segmentation"].toObject(), well);
  }

  return well;
}