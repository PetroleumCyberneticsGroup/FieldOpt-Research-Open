/******************************************************************************
 *
 * model.cpp
 *
 * Created: 28.09.2015 2015 by einar
 *
 * This file is part of the FieldOpt project.
 *
 * Copyright (C) 2015-2015 Einar J.M. Baumann <einar.baumann@ntnu.no>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 *****************************************************************************/

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <boost/lexical_cast.hpp>
#include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>
#include "model.h"
#include "settings_exceptions.h"
#include "Utilities/filehandling.hpp"
#include "deck_parser.h"
#include "trajectory_importer.h"

namespace Settings {

Model::Model(QJsonObject json_model, Paths &paths)
{
    // Reservoir
    if (!paths.IsSet(Paths::ENSEMBLE_FILE)) {
        try {
            QJsonObject json_reservoir = json_model["Reservoir"].toObject();
            readReservoir(json_reservoir, paths);
        }
        catch (std::exception const &ex) {
            throw UnableToParseReservoirModelSectionException(
                "Unable to parse reservoir model section: " + std::string(ex.what()));
        }
    }

    // Control times
    if (!json_model.contains("ControlTimes") || !json_model["ControlTimes"].isArray())
        throw UnableToParseModelSectionException(
            "The ControlTimes array must be defined with at least one time for the model.");
    control_times_ = QList<int>();
    if (json_model.contains("NPVInterval")) {
        if (json_model["NPVInterval"].toString().compare("Yearly") == 0) {
            if (json_model.contains("NPVYears")) {
                for (int i = 0; i < json_model["NPVYears"].toInt(); ++i) {
                    control_times_.append(365 * i);
                }
            } else {
                throw UnableToParseModelSectionException("Unable to parse NPVYears");
            }
        } else if (json_model["NPVInterval"].toString().compare("Monthly") == 0) {
            if (json_model.contains("NPVMonths")) {
                for (int i = 0; i < json_model["NPVMonths"].toInt(); ++i) {
                    control_times_.append(30 * i);
                }
            } else {
                throw UnableToParseModelSectionException("Unable to parse NPVMonths");
            }
        }
    }
    for (int i = 0; i < json_model["ControlTimes"].toArray().size(); ++i) {
        if (!control_times_.contains(json_model["ControlTimes"].toArray().at(i).toInt())) {
            control_times_.append(json_model["ControlTimes"].toArray().at(i).toInt());
        }
    }
    qSort(control_times_);

    // Wells
    wells_ = QList<Well>();
    if (json_model.contains("Import")) {
        auto json_import = json_model["Import"].toObject();
        if (!paths.IsSet(Paths::SIM_DRIVER_FILE)) {
            throw std::runtime_error("SchedulePath must be specified (relative to DriverPath) to use the Import feature.");
        }
        if (VERB_SET >= 2) {
            Printer::ext_info("Parsing schedule ...", "Settings", "Model");
        }
        deck_parser_ = new DeckParser(paths.GetPath(Paths::SIM_DRIVER_FILE));

        if (!json_import["Keywords"].toArray().contains(QJsonValue("AllWells"))) {
            throw std::runtime_error("Unable to import simulator schedule. Import Keywords array does not contain"
                                         "any recognized keywords.");
        }
        if (VERB_SET >= 2) {
            Printer::ext_info("Importing wells ...", "Settings", "Model");
        }
        wells_.append(deck_parser_->GetWellData());
        setImportedWellDefaults(json_import);
        parseImportedWellOverrides(json_model["Wells"].toArray());

        // Spline conversion
        if (json_import.contains("SplineConversion")) {
            auto json_sconv = json_import["SplineConversion"].toObject();
            if (!json_sconv.contains("Wells")) {
                throw std::runtime_error("A list of wells must be specified when the SplineConversion is used.");
            }
            for (auto jwell : json_sconv["Wells"].toArray()) {
                QString wname = jwell.toString();
                int widx = -1;
                for (int i = 0; i < wells_.size(); ++i) {
                    if (QString::compare(wells_[i].name, wname) == 0) {
                        widx = i;
                        break;
                    }
                }
                if (widx == -1) {
                    throw std::runtime_error("Well " + wname.toStdString() +
                        " defined in SplineConversion well list not found.");
                }
                wells_[widx].convert_well_blocks_to_spline = true;
                wells_[widx].definition_type = WellDefinitionType::WellSpline;
                if (json_sconv.contains("SplinePoints") && json_sconv["SplinePoints"].toInt() >= 2) {
                    wells_[widx].n_spline_points = json_sconv["SplinePoints"].toInt();
                    if (VERB_SET >= 2) {
                        Printer::info("Setting number of spline points for "
                                              + boost::lexical_cast<std::string>(wname.toStdString()) + " to "
                                              + boost::lexical_cast<std::string>(wells_[widx].n_spline_points)
                                              + " for conversion.");
                    }
                }
                else {
                    if (VERB_SET >= 2) {
                        Printer::info("Defaulting number of spline points for "
                                          + boost::lexical_cast<std::string>(wname.toStdString()) + " to 2 "
                                          + " for conversion.");
                    }
                }
            }
        }

        // Trajectory import
        if (json_import.contains("ImportTrajectories")) {
            auto json_import_well_names = json_import["ImportTrajectories"].toArray();
            std::vector<std::string> import_well_names;
            for (auto jwn : json_import_well_names) {
                import_well_names.push_back(jwn.toString().toStdString());
            }
            std::string trajectories_path;
            if (!paths.IsSet(Paths::TRAJ_DIR)) {
                trajectories_path = paths.GetPath(Paths::SIM_DRIVER_DIR) + "/trajectories";
                assert(DirectoryExists(trajectories_path));
            }
            else {
                trajectories_path = paths.GetPath(Paths::TRAJ_DIR);
            }
            auto traj_importer = TrajectoryImporter(trajectories_path, import_well_names);

            // set list in well objects
            for (auto wname : import_well_names) {
                for (int i = 0; i < wells_.size(); ++i) {
                    if (QString::compare(wells_[i].name, QString::fromStdString(wname)) == 0) {
                        wells_[i].imported_wellblocks_ = traj_importer.GetImportedTrajectory(wname);
                        wells_[i].definition_type = WellDefinitionType::WellSpline;
                        wells_[i].convert_well_blocks_to_spline = false;
                    }
                }
            }
        }

        // Segmentation
        if (json_model.contains("Wells") && json_model["Wells"].isArray()) {
            for (auto jwell : json_model["Wells"].toArray()) { // Go through list of wells in json file
                if (jwell.toObject().contains("Segmentation")) { // Check if the well has the Segmentation keyword
                    for (int j = 0; j < wells_.size(); ++j) { // Loop through parsed wells
                        if (QString::compare(wells_[j].name, jwell.toObject()["Name"].toString()) == 0) { // Check if well names match
                            wells_[j].use_segmented_model = true;
                            parseSegmentation(jwell.toObject()["Segmentation"].toObject(), wells_[j]);
                            break;
                        }
                    }
                }
            }

        }
    }
    else {
        try {
            QJsonArray json_wells = json_model["Wells"].toArray();
            for (int i = 0; i < json_wells.size(); ++i) {
                QJsonObject json_well = json_wells[i].toObject();
                wells_.append(readSingleWell(json_well));
            }
        }
        catch (std::exception const &ex) {
            throw UnableToParseWellsModelSectionException(
                "Unable to parse wells model section: " + std::string(ex.what()));
        }
    }
}

void Model::readReservoir(QJsonObject json_reservoir, Paths &paths)
{
    // Reservoir grid source type
    QString type = json_reservoir["Type"].toString();
    if (QString::compare(type, "ECLIPSE") == 0)
        reservoir_.type = ReservoirGridSourceType::ECLIPSE;
    else throw UnableToParseReservoirModelSectionException("Grid source type " + type.toStdString() +  "not recognized.");

    // Reservoir grid path
    if (!paths.IsSet(Paths::GRID_FILE) && json_reservoir.contains("Path")) {
        paths.SetPath(Paths::GRID_FILE, json_reservoir["Path"].toString().toStdString());
    }
}

Model::Well Model::readSingleWell(QJsonObject json_well)
{
    Well well;

    // Well name
    well.name = json_well["Name"].toString();

    // Well Type
    QString type = json_well["Type"].toString();

    // Well group
    if (json_well.contains("Group"))
        well.group = json_well["Group"].toString();
    else well.group = "";

    if (QString::compare(type, "Producer") == 0)
        well.type = WellType::Producer;
    else if (QString::compare(type, "Injector") == 0)
        well.type = WellType::Injector;
    else throw UnableToParseWellsModelSectionException("Well type " + type.toStdString() + " not recognized for well " + well.name.toStdString());

    // Well definition type
    QString definition_type = json_well["DefinitionType"].toString();
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
            if (json_well_block_i.contains("IsVariable") && json_well_block_i["IsVariable"].toBool() == true)
                block.is_variable = true;
            else block.is_variable = false;

            if (json_well_block_i.contains("Completion")) {
                QJsonObject json_well_block_i_completion = json_well_block_i["Completion"].toObject();
                block.completion.type = WellCompletionType::Perforation;
                if (json_well_block_i_completion.contains("TransmissibilityFactor"))
                    block.completion.transmissibility_factor = json_well_block_i_completion["TransmissibilityFactor"].toDouble();
                if (json_well_block_i_completion.contains("IsVariable"))
                    block.completion.is_variable = json_well_block_i_completion["IsVariable"].toBool();
                else
                    block.completion.is_variable = false;
                block.completion.name = "Transmissibility#" + well.name + "#" + QString::number(i);
                block.has_completion = true;
            }
            else {
                block.has_completion = true;
                block.completion.type = WellCompletionType::Perforation;
                block.completion.transmissibility_factor = -1;
                block.completion.is_variable = false;
            }
            block.name = "WellBlock#" + well.name + "#" + QString::number(i);
            well.well_blocks.append(block);
        }
    }
    else if (QString::compare(definition_type, "WellSpline") == 0) {
        well.definition_type = WellDefinitionType::WellSpline;
        if (json_well.contains("UseBezierSpline") && json_well["UseBezierSpline"].toBool() == true) {
            well.use_bezier_spline = true;
        }
        else {
            well.use_bezier_spline = false;
        }
        if (json_well.contains("SplinePointArray")) {
            auto json_sp_array = json_well["SplinePointArray"].toArray();
            assert(json_sp_array.size() >=2);
            for (int i = 0; i < json_sp_array.size(); ++i) {
                auto json_point = json_sp_array[i].toObject();
                Well::SplinePoint point;
                point.x = json_point["x"].toDouble();
                point.y = json_point["y"].toDouble();
                point.z = json_point["z"].toDouble();
                point.is_variable = json_point["IsVariable"].toBool();
                if (i == 0) {
                    point.name = "SplinePoint#" + well.name + "#heel";
                }
                else if (i == json_sp_array.size() - 1) {
                    point.name = "SplinePoint#" + well.name + "#toe";
                }
                else {
                    point.name = "SplinePoint#" + well.name + "#P" + QString::number(i);
                }
                well.spline_points.push_back(point);
            }
            well.spline_heel = well.spline_points.front();
            well.spline_toe = well.spline_points.back();
        }
        else {
            QJsonObject json_points = json_well["SplinePoints"].toObject();
            if (!json_points.contains("Heel") || !json_points.contains("Toe"))
                throw UnableToParseWellsModelSectionException("Both Heel and Toe must be defined for spline-type wells.");

            QJsonObject json_heel = json_points["Heel"].toObject();
            QJsonObject json_toe = json_points["Toe"].toObject();
            well.spline_heel.x = json_heel["x"].toDouble();
            well.spline_heel.y = json_heel["y"].toDouble();
            well.spline_heel.z = json_heel["z"].toDouble();
            well.spline_toe.x = json_toe["x"].toDouble();
            well.spline_toe.y = json_toe["y"].toDouble();
            well.spline_toe.z = json_toe["z"].toDouble();
            if (json_heel.contains("IsVariable") && json_heel["IsVariable"].toBool())
                well.spline_heel.is_variable = true;
            else well.spline_heel.is_variable = false;
            if (json_toe.contains("IsVariable") && json_toe["IsVariable"].toBool())
                well.spline_toe.is_variable = true;
            else well.spline_toe.is_variable = false;
            if ((well.spline_heel.is_variable && well.spline_toe.is_variable)
                || (json_points.contains("IsVariable") && json_points["IsVariable"].toBool() == true)) {
                well.is_variable_spline = true;
            }
            well.spline_heel.name = "SplinePoint#" + well.name + "#heel";
            well.spline_toe.name = "SplinePoint#" + well.name + "#toe";

            well.spline_points.push_back(well.spline_heel);
            if (json_points.contains("AdditionalPoints")) {
                int interp_points = json_points["AdditionalPoints"].toInt();
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
    }
    else if (QString::compare(definition_type, "PseudoContVertical2D") == 0) {
        QJsonObject json_position = json_well["Position"].toObject();
        well.definition_type = WellDefinitionType::PseudoContVertical2D;
        well.pseudo_cont_position.i = json_position["i"].toInt();
        well.pseudo_cont_position.j = json_position["j"].toInt();
        if (json_position.contains("IsVariable") && json_position["IsVariable"].toBool() == true) {
            well.pseudo_cont_position.is_variable = true;
        }
        else well.spline_heel.is_variable = false;
    }
    else throw UnableToParseWellsModelSectionException("Well definition type " + definition_type.toStdString() + " not recognized for well " + well.name.toStdString());

    // Wellbore radius
    if (!json_well.contains("WellboreRadius"))
        throw UnableToParseWellsModelSectionException("The wellbore radius must be defined for all wells.");
    well.wellbore_radius = json_well["WellboreRadius"].toDouble();

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

        if (!controlTimeIsDeclared(json_controls.at(i).toObject()["TimeStep"].toInt()))
            throw UnableToParseWellsModelSectionException("All time steps must be declared in the ControlTimes array. Inconsistency detected in Controls declaration.");
        else control.time_step = json_controls.at(i).toObject()["TimeStep"].toInt();

        // State (Open or shut)
        if (json_controls[i].toObject().contains("State") && QString::compare("Shut", json_controls.at(i).toObject()["State"].toString()) == 0)
            control.state = WellState::WellShut;
        else
            control.state = WellState::WellOpen;


        // Control mode
        if (QString::compare("BHP", json_controls.at(i).toObject()["Mode"].toString()) == 0) {
            control.control_mode = ControlMode::BHPControl;
            control.bhp = json_controls.at(i).toObject()["BHP"].toDouble();
            control.name = "BHP#" + well.name + "#" + QString::number(control.time_step);
        }
        else if (QString::compare("Rate", json_controls.at(i).toObject()["Mode"].toString()) == 0) {
            control.control_mode = ControlMode::RateControl;
            control.rate = json_controls.at(i).toObject()["Rate"].toDouble();
            control.name = "Rate#" + well.name + "#" + QString::number(control.time_step);
        }
        else throw UnableToParseWellsModelSectionException("Well control type " + json_controls.at(i).toObject()["Mode"].toString().toStdString() + " not recognized for well " + well.name.toStdString());

        // Injection type
        if (well.type == WellType::Injector) {
            if (!json_controls.at(i).toObject().contains("Type"))
                control.injection_type = InjectionType::WaterInjection;
            if (QString::compare("Water", json_controls.at(i).toObject()["Type"].toString()) == 0)
                control.injection_type = InjectionType::WaterInjection;
            else if (QString::compare("Gas", json_controls.at(i).toObject()["Type"].toString()) == 0)
                control.injection_type = InjectionType::GasInjection;
        }
        if (json_controls[i].toObject()["IsVariable"].toBool())
            control.is_variable = true;
        else
            control.is_variable = false;
        well.controls.append(control);
    }

    // Preferred Phase
    if (QString::compare("Oil", json_well["PreferredPhase"].toString()) == 0)
        well.preferred_phase = PreferredPhase::Oil;
    else if (QString::compare("Water", json_well["PreferredPhase"].toString()) == 0)
        well.preferred_phase = PreferredPhase::Water;
    else if (QString::compare("Gas", json_well["PreferredPhase"].toString()) == 0)
        well.preferred_phase = PreferredPhase::Gas;
    else if (QString::compare("Liquid", json_well["PreferredPhase"].toString()) == 0)
        well.preferred_phase = PreferredPhase::Liquid;

    // Segmentation
    if (json_well.contains("Segmentation")) {
        well.use_segmented_model = true;
        parseSegmentation(json_well["Segmentation"].toObject(), well);
    }

    return well;
}

bool Model::controlTimeIsDeclared(int time) const
{
    return control_times_.contains(time);
}

std::string Model::Well::ControlEntry::toString() {
    std::stringstream ce;
    ce << "ControlEntry - Timestep:  " << time_step << "\n";
    ce << "               State:     " << (state == WellState::WellOpen ? "Open" : "Shut") << "\n";
    ce << "               Mode:      " << (control_mode == ControlMode::RateControl ? "Rate" : "BHP") << "\n";
    ce << "               BHP:       " << boost::lexical_cast<std::string>(bhp) << "\n";
    ce << "               Rate:      " << boost::lexical_cast<std::string>(rate) << "\n";
    ce << "               Inj. type: " << (injection_type == InjectionType::WaterInjection ? "Water" : "Gas/UNKWN") << "\n";
    ce << "               Variable:  " << (is_variable ? "Yes" : "No") << "\n";
    return ce.str();
}

std::string Model::Well::toString() {
    std::stringstream ce;
    ce << "Well - Name:           " << name.toStdString() << "\n";
    ce << "       Type:           " << (type == WellType::Injector ? "Injector" : "Producer") << "\n";
    ce << "       Group:          " << group.toStdString() << "\n";
    ce << "       Radius:         " << wellbore_radius << "\n";
    ce << "       Direction:      " << direction << "\n";
    ce << "       Pref. phase:    " << preferred_phase << "\n";
    ce << "       Def. type:      " << definition_type << "\n";
    ce << "       N. well blocks: " << well_blocks.size() << "\n";
    ce << "       N. controls:    " << controls.size() << "\n";
    return ce.str();
}

bool Model::Well::ControlEntry::isDifferent(ControlEntry other) {
    if (state != other.state)
        return true;
    if (rate != other.rate)
        return true;
    if (bhp != other.bhp)
        return true;
    if (control_mode != other.control_mode)
        return true;
    return false; // Assume they're equal if none of the above hits.
}

void Model::setImportedWellDefaults(QJsonObject json_import) {
    if (
        !json_import.contains("ProducerDefaultBHP") ||
        (!json_import.contains("InjectorDefaultRate") &&
        !json_import.contains("InjectorDefaultBHP"))
        ) {
        throw std::runtime_error("When importing from schedule, you must provide both the "
                                     "InjectorDefaultRate or InjectorDefaultBHP and "
                                     "ProducerDefaultBHP properties in the Import object.");
    }

    for (int i = 0; i < wells_.size(); ++i) {
        wells_[i].controls = QList<Well::ControlEntry>{wells_[i].controls[0]}; // Remove all but first control
        wells_[i].controls[0].time_step = getClosestControlTime(wells_[i].controls[0].time_step);
        if (wells_[i].type == Injector) {
            if (json_import.contains("InjectorDefaultRate")) {
                wells_[i].controls[0].control_mode = RateControl;
                wells_[i].controls[0].rate = json_import["InjectorDefaultRate"].toDouble();
            }
            else if (json_import.contains("InjectorDefaultBHP")) {
                wells_[i].controls[0].control_mode = BHPControl;
                wells_[i].controls[0].bhp = json_import["InjectorDefaultBHP"].toDouble();
            }
        }
        else { // Producer
            wells_[i].controls[0].control_mode = BHPControl;
            wells_[i].controls[0].bhp = json_import["ProducerDefaultBHP"].toDouble();
        }
    }
}

void Model::parseImportedWellOverrides(QJsonArray json_wells) {
    for (auto w : json_wells) {

        auto well = w.toObject();
        auto json_control = well["ControlOverrides"].toObject();
        if (json_control.contains("IsVariable") && json_control["IsVariable"].toBool() == true) {
            for (int i = 0; i < wells_.size(); ++i) { // Find the correct well
                if (QString::compare(wells_[i].name, well["Name"].toString()) ==0){
                    auto control = wells_[i].controls[0];
                    control.is_variable = true;
                    auto name_root_lst = control.name.split("#");
                    name_root_lst.removeLast();
                    auto name_root = name_root_lst.join("#") + "#";
                    auto new_control_set = QList<Well::ControlEntry>();
                    for (auto time : json_control["VariableTimeSteps"].toArray()) {
                        int t = time.toInt();
                        if (!control_times_.contains(t)) {
                            throw std::runtime_error("All entries in VariableTimeSteps must match the overall"
                                                         "ControlTimes vector. " + std::to_string(t) + " does not match");
                        }
                        auto control_copy = control;
                        control_copy.time_step = t;
                        control_copy.name = name_root + QString::number(t);
                        new_control_set.append(control_copy);
                    }
                    wells_[i].controls = new_control_set;
                    break;
                }
                else {
                    if (i == wells_.size() - 1) { // Passed the final well without finding a match.
                        throw std::runtime_error("Unable to find matching well for ControlOverrides");
                    }
                }
            }
        }
    }

}

int Model::getClosestControlTime(int deck_time) {
    auto diff = std::vector<int>(control_times_.size());
    for (int i = 0; i < control_times_.size(); ++i) { // Compute control times
        diff[i] = std::abs(deck_time - control_times_[i]);
    }
    auto idx_of_min_element = std::distance(diff.begin(), std::min_element(diff.begin(), diff.end()));
    return control_times_[idx_of_min_element];
}

void Model::parseSegmentation(QJsonObject json_seg, Well &well) {
    parseSegmentTubing(json_seg, well);
    parseSegmentAnnulus(json_seg, well);
    parseSegmentCompartments(json_seg, well);
}
void Model::parseSegmentTubing(const QJsonObject &json_seg, Model::Well &well) const {
    if (VERB_SET >= 2) {
        Printer::ext_info("Parsing Tubing ...", "Settings", "Model");
    }
    if (json_seg.contains("Tubing")) {
        try {
            well.seg_tubing.diameter = json_seg["Tubing"].toObject()["Diameter"].toDouble();
            well.seg_tubing.roughness = json_seg["Tubing"].toObject()["Roughness"].toDouble();
        }
        catch ( ... ) {
            throw std::runtime_error("For Tubing, both Diameter and Roughness must be defined.");
        }
    }
    else {
        if (VERB_SET >= 1) {
            Printer::ext_info("Tubing keyword not found in Segmentation. Defaulting Diameter to 0.1 and "
                                  "Roughness to 1.52E-5.", "Settings", "Model");
        }
        well.seg_tubing.diameter = 0.1;
        well.seg_tubing.roughness = 1.52E-5;
    }
    well.seg_tubing.cross_sect_area = M_PI / 4.0 * well.seg_tubing.diameter * well.seg_tubing.diameter;
}

void Model::parseSegmentAnnulus(const QJsonObject &json_seg, Model::Well &well) const {
    if (VERB_SET >= 2) {
        Printer::ext_info("Parsing Annulus ...", "Settings", "Model");
    }
    if (json_seg.contains("Annulus")) {
        try {
            well.seg_annulus.diameter = json_seg["Annulus"].toObject()["Diameter"].toDouble();
            well.seg_annulus.roughness = json_seg["Annulus"].toObject()["Roughness"].toDouble();
            well.seg_annulus.cross_sect_area = json_seg["Annulus"].toObject()["CrossSectionArea"].toDouble();
        }
        catch ( ... ) {
            throw std::runtime_error("For Annulus, both Diameter, CrossSectionArea and Roughness must be defined.");
        }
    }
    else {
        if (VERB_SET >= 1) {
            Printer::ext_info("Annulus keyword not found in Segmentation. Defaulting Diameter to 0.04, "
                                  "Ac to 8.17E-3 and Roughness to 1.52E-5.", "Settings", "Model");
        }
        well.seg_annulus.diameter = 0.04;
        well.seg_annulus.roughness = 1.52E-5;
        well.seg_annulus.cross_sect_area = 8.17E-3;
    }
}

void Model::parseSegmentCompartments(const QJsonObject &json_seg, Model::Well &well) const {
    if (VERB_SET >= 2) {
        Printer::ext_info("Parsing Compartments ...", "Settings", "Model");
    }
    if (json_seg.contains("Compartments")) {
        auto json_compts = json_seg["Compartments"].toObject();
        try {
            well.seg_n_compartments = json_compts["Count"].toInt();

            if (json_compts.contains("VariablePackers") && json_compts["VariablePackers"].toBool() == true) {
                well.seg_compartment_params.variable_placement = true;
            }
            if (json_compts.contains("VariableICDs") && json_compts["VariableICDs"].toBool() == true) {
                well.seg_compartment_params.variable_strength = true;
            }
            if (json_compts.contains("ICDType")) {
                well.seg_compartment_params.type = WellCompletionType::ICV;
            }
            if (json_compts.contains("ICDValveSize")) {
                well.seg_compartment_params.valve_size = json_compts["ICDValveSize"].toDouble();
            }
            else {
                if (VERB_SET >= 1) {
                    Printer::ext_info("ICDValveSize keyword not found in Compartments. Defaulting to 7.85E-5.",
                                      "Settings", "Model");
                }
                well.seg_compartment_params.valve_size = 7.85E-5;
            }
            if (json_compts.contains("ICDValveFlowCoeff")) {
                well.seg_compartment_params.valve_flow_coeff = json_compts["ICDValveFlowCoeff"].toDouble();
            }
            else {
                if (VERB_SET >= 1) {
                    Printer::ext_info("ICDValveFlowCoeff keyword not found in Compartments. Defaulting to 7.85E-5.",
                                      "Settings", "Model");
                }
                well.seg_compartment_params.valve_flow_coeff = 0.66;
            }
        }
        catch (...) {
            throw std::runtime_error("Something went wrong while parsing the Compartments section.");
        }
    }
    else {
        throw std::runtime_error("The Compartments keyword must be specified when using the Segmentation keyword.");
    }
}


}

