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

#include "wellspline.h"
#include <iostream>
#include <time.h>
#include <cmath>
#include <wells/well_exceptions.h>
#include <QtCore/QDateTime>
#include <Utilities/time.hpp>
#include <boost/algorithm/string.hpp>
#include <QList>
#include <Utilities/verbosity.h>
#include <Utilities/printer.hpp>
#include "Model/cpp-spline/src/Bezier.h"
#include "WellIndexCalculation/wicalc_rixx.h"

namespace Model {
namespace Wells {
namespace Wellbore {

using namespace Reservoir::WellIndexCalculation;
using Printer::info;
using Printer::ext_info;
using Printer::num2str;
using Printer::ext_warn;
using std::string;
using std::to_string;

WellSpline::WellSpline(Settings::Model::Well well_settings,
                       Properties::VarPropContainer *variable_container,
                       Reservoir::Grid::Grid *grid,
                       Reservoir::WellIndexCalculation::wicalc_rixx *wic) {
  vp_ = well_settings.verbParams();
  grid_ = grid;
  assert(grid_ != nullptr);
  well_settings_ = well_settings;
  is_variable_ = false;
  use_bezier_spline_ = well_settings.use_bezier_spline;

  // Initialize WIC if this is the first spline well initialized.
  if (wic == nullptr) {
    wic = new Reservoir::WellIndexCalculation::wicalc_rixx(well_settings_, grid_);
    wic_ = wic;
  } else { // If not, use existing WIC object.
    wic_ = wic;
  }

  if (!well_settings.imported_wellblocks_.empty()) { // Imported blocks present
    spline_points_from_import(well_settings);
    imported_wellblocks_ = well_settings.imported_wellblocks_;
  }

  if (vp_.vMOD >= 2) {
    auto wn = well_settings.name.toStdString();
    auto sp = num2str(well_settings.spline_points.size());
    auto sn = well_settings.spline_points[0].name.toStdString();
    im_ = "Init well spline for well " + wn + ". N points: " + sp;
    im_ += "; First spline point name: " + sn;
    ext_info(im_, md_, cl_,vp_.lnw);
  }

  for (auto point : well_settings.spline_points) {
    if (vp_.vMOD >= 2) {
      auto wn = well_settings.name.toStdString();
      im_ = "Adding new spline point for well " + wn + ": ";
      im_ += num2str(point.x) + ", " + num2str(point.y) + ", " + num2str(point.z);
      im_ += " (" + point.name.toStdString() + ")";
      ext_info(im_,md_, cl_, vp_.lnw);
    }
    auto *pt = new SplinePoint();
    pt->x = new ContinuousProperty(point.x);
    pt->y = new ContinuousProperty(point.y);
    pt->z = new ContinuousProperty(point.z);
    assert(point.name.size() > 0);
    pt->x->setName(point.name + "#x");
    pt->y->setName(point.name + "#y");
    pt->z->setName(point.name + "#z");
    if (point.is_variable) {
      is_variable_ = true;
      variable_container->AddVariable(pt->x);
      variable_container->AddVariable(pt->y);
      variable_container->AddVariable(pt->z);
    }
    spline_points_.push_back(pt);
  }
  seconds_spent_in_compute_wellblocks_ = 0;

  last_computed_grid_ = "";
  last_computed_spline_ = std::vector<Eigen::Vector3d>();
}

WellSpline::WellSpline() {
  last_computed_grid_ = "";
  last_computed_spline_ = std::vector<Eigen::Vector3d>();
}

void WellSpline::spline_points_from_import(Settings::Model::Well &well_settings) {
  QString name_base = "SplinePoint#" + well_settings.name + "#";
  well_settings_.spline_points = QList<Settings::Model::Well::SplinePoint> ();
  auto first_spline_point = Settings::Model::Well::SplinePoint();
  first_spline_point.x = well_settings.imported_wellblocks_[0].in().x();
  first_spline_point.y = well_settings.imported_wellblocks_[0].in().y();
  first_spline_point.z = well_settings.imported_wellblocks_[0].in().z();
  if (well_settings.is_variable_spline) {
    first_spline_point.is_variable = true;
  }
  first_spline_point.name =  name_base + "heel";

  well_settings.spline_points.push_back(first_spline_point);
  int i = 1;
  for (auto block : well_settings.imported_wellblocks_) {
    auto spline_point = Settings::Model::Well::SplinePoint();
    spline_point.x = block.out().x();
    spline_point.y = block.out().y();
    spline_point.z = block.out().z();
    if (well_settings.is_variable_spline) {
      spline_point.is_variable = true;
    }
    spline_point.name = name_base + "P" + QString::number(i++);
    well_settings.spline_points.push_back(spline_point);
  }
  well_settings.spline_points.back().name = name_base + "heel";
}

QList<WellBlock *> *WellSpline::computeWellBlocks() {
  assert(spline_points_.size() >= 2);
  assert(grid_ != nullptr && grid_ != 0);

  if (vp_.vMOD >= 2) {
    std::string points_str = "";
    for (auto pt : spline_points_) {
      auto point = pt->ToEigenVector();
      std::stringstream point_str;
      point_str << "(" << point.x() << ", " << point.y() << ", " << point.z() << "), ";
      points_str += point_str.str();
    }

    im_ = "Starting well index calculation. Points: ";
    im_ += points_str + "Grid: " + grid_->GetGridFilePath();
    ext_info(im_, md_, cl_, vp_.lnw);
  }

  last_computed_grid_ = grid_->GetGridFilePath();
  last_computed_spline_ = create_spline_point_vector();

  WellDefinition welldef;
  welldef.wellname = well_settings_.name.toStdString();

  auto spline_points = getPoints();
  for (int w = 0; w < spline_points.size() - 1; ++w) {
    welldef.radii.push_back(well_settings_.wellbore_radius);
    welldef.skins.push_back(0.0);
    welldef.skins.push_back(0.0);
    welldef.heels.push_back(spline_points[w]);
    welldef.toes.push_back(spline_points[w+1]);
    if (welldef.heel_md.empty()) {
      welldef.heel_md.push_back(0.0);
    }
    else {
      double prev_toe = welldef.toe_md.back();
      welldef.heel_md.push_back(prev_toe);
    }
    welldef.toe_md.push_back(
      welldef.heel_md.back() + (welldef.toes.back() - welldef.heels.back()).norm()
    );
  }

  auto start = QDateTime::currentDateTime();
  vector<IntersectedCell> block_data;
  if (imported_wellblocks_.empty() || is_variable_) {
    wic_->ComputeWellBlocks(block_data, welldef);
  }
  else {
    if (vp_.vMOD >= 1) {
      wm_ = "Well index calculation for imported ";
      wm_ += "paths is not properly implemented at this time.";
      ext_warn(wm_, md_, cl_, vp_.lnw);
    }
    block_data = convertImportedWellblocksToIntersectedCells();
    for (int i = 0; i < block_data.size(); ++i) {
      wic_->ComputeWellBlocks(block_data, welldef);
    }
  }
  auto end = QDateTime::currentDateTime();
  seconds_spent_in_compute_wellblocks_ = time_span_seconds(start, end);

  auto *blocks = new QList<WellBlock *>();
  for (auto & i : block_data) {
    blocks->append(getWellBlock(i));
    blocks->last()->setEntryPoint(i.get_segment_entry_point(0));
    blocks->last()->setExitPoint(i.get_segment_exit_point(0));
    blocks->last()->setEntryMd(i.get_segment_entry_md(0));
    blocks->last()->setExitMd(i.get_segment_exit_md(0));
    blocks->last()->setLength(i.get_segment_length(0));
  }
  if (blocks->empty()) {
    throw WellBlocksNotDefined("WIC could not compute.");
  }

  if (vp_.vMOD >= 2) {
    im_ = "Computation time for WIs [secs]: ";
    im_ += to_string(seconds_spent_in_compute_wellblocks_);
    info(im_, vp_.lnw);
  }

  if (vp_.vMOD >= 2) {
    im_ = "Computed " + num2str(blocks->size()) + " well blocks from ";
    im_ += num2str(block_data.size()) + " intersected cells for well " + welldef.wellname;
    ext_info(im_, md_, cl_, vp_.lnw);
  }
  return blocks;

}

QList<WellBlock *> *WellSpline::GetWellBlocks() {
  return computeWellBlocks();
}

WellBlock *WellSpline::getWellBlock(Reservoir::WellIndexCalculation::IntersectedCell block_data) {
  if (vp_.vMOD >= 2) {
    im_ = "Creating WellBlock for IC " + block_data.ijk_index().to_string();
    im_ += " with WI " + num2str(block_data.cell_well_index_matrix());
    cout << im_ << endl;
  }
  auto wb = new WellBlock(block_data.ijk_index().i()+1,
                          block_data.ijk_index().j()+1,
                          block_data.ijk_index().k()+1);
  auto comp = new Completions::Perforation();
  comp->setTransmissibility_factor(block_data.cell_well_index_matrix());
  wb->AddCompletion(comp);
  return wb;
}

std::vector<Eigen::Vector3d> WellSpline::create_spline_point_vector() const {
  std::vector<Eigen::Vector3d> spline_points;
  for (auto point : spline_points_) {
    spline_points.push_back(point->ToEigenVector());
  }
  return spline_points;
}

bool WellSpline::HasGridChanged() const {
  return last_computed_grid_.empty() || !boost::equals(last_computed_grid_, grid_->GetGridFilePath());
}

bool WellSpline::HasSplineChanged() const {
  if (last_computed_spline_.empty()) {
    return true;
  }

  std::vector<Eigen::Vector3d> new_spline_points;
  for (auto point : spline_points_) {
    new_spline_points.push_back(point->ToEigenVector());
  }
  assert(new_spline_points.size() == last_computed_spline_.size());

  double point_difference_sum = 0;
  for (int i = 0; i < last_computed_spline_.size(); ++i) {
    point_difference_sum += std::abs((last_computed_spline_[i] - new_spline_points[i]).norm());
  }
  return point_difference_sum > 1e-7;
}

std::vector<Reservoir::WellIndexCalculation::IntersectedCell>
  WellSpline::convertImportedWellblocksToIntersectedCells() {
  auto intersected_cells = vector<IntersectedCell>();
  double md_in = 0.0;
  double md_out = 0.0;
  double length = 0.0;
  for (auto iwb : imported_wellblocks_) {
    auto cell = grid_->GetCell(iwb.ijk().x()-1, iwb.ijk().y()-1, iwb.ijk().z()-1);
    auto ic = Reservoir::WellIndexCalculation::IntersectedCell(cell);
    md_out = md_in + (iwb.out() - iwb.in()).norm();
    length = (iwb.out() - iwb.in()).norm();
    ic.add_new_segment(iwb.in(), iwb.out(), md_in, md_out, length, well_settings_.wellbore_radius, 0.0);
    intersected_cells.push_back(ic);
    md_in = md_out;
  }
  return intersected_cells;
}

Eigen::Vector3d WellSpline::SplinePoint::ToEigenVector() const {
  return Eigen::Vector3d(this->x->value(), this->y->value(), this->z->value());
}

void WellSpline::SplinePoint::FromEigenVector(const Eigen::Vector3d vec) {
  this->x->setValue(vec.x());
  this->y->setValue(vec.y());
  this->z->setValue(vec.z());
}

vector<Eigen::Vector3d> WellSpline::getPoints() const {
  if (use_bezier_spline_) {
    return convertToBezierSpline();
  } else {
    vector<Eigen::Vector3d> points;
    for (auto point : spline_points_) {
      points.push_back(point->ToEigenVector());
    }
    return points;
  }
}

vector<Eigen::Vector3d> WellSpline::convertToBezierSpline() const {
  if (VERB_MOD >= 2) {
    string im = "Generating bezier spline for well ";
    im += well_settings_.name.toStdString() + ". N original points: ";
    im += Printer::num2str(spline_points_.size());
    ext_info(im, md_, cl_, vp_.lnw);
  }
  assert(spline_points_.size() >= 4);
  Curve *curve = new Bezier();
  curve->set_steps(50);
  for (int j = 0; j < spline_points_.size(); ++j) {
    curve->add_way_point(Vector(spline_points_[j]->x->value(),
                                spline_points_[j]->y->value(),
                                spline_points_[j]->z->value()));
  }
  vector<Eigen::Vector3d> bezier_points;
  for (int i = 0; i < curve->node_count(); ++i) {
    bezier_points.push_back(Eigen::Vector3d(curve->node(i).x,
                                            curve->node(i).y,
                                            curve->node(i).z));
  }
  delete curve;
  return bezier_points;

}
}
}
}
