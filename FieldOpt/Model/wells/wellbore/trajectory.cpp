/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#include "trajectory.h"
#include "Model/wells/well_exceptions.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "Utilities/verbosity.h"
#include "Utilities/printer.hpp"
#include "polar_spline.h"

namespace Model {
namespace Wells {
namespace Wellbore {

Trajectory::
Trajectory(Settings::Model::Well well_settings,
           Properties::VarPropContainer *variable_container,
           ::Reservoir::Grid::Grid *grid,
           Reservoir::WellIndexCalculation::wicalc_rixx *wic) {

  wsettings_ = well_settings;
  string wn = wsettings_.name.toStdString();
  vp_ = wsettings_.verbParams();

  well_blocks_ = new QList<WellBlock *>();
  well_spline_ = nullptr;
  pseudo_cont_vert_ = nullptr;
  definition_type_ = well_settings.definition_type;

  if (well_settings.definition_type == WType::WellBlocks) {
    if (vp_.vMOD >= 5) { ext_warn("BLOCK WELL: " + wn); }
    initializeWellBlocks(well_settings, variable_container);

  } else if (well_settings.definition_type == WType::WellSpline) {
    if (vp_.vMOD >= 5) { ext_warn("WELLSPLINE: " + wn); }
    if (well_settings.convert_well_blocks_to_spline) {
      convertWellBlocksToWellSpline(well_settings, grid);
    }
    well_spline_ = new WellSpline(well_settings,
                                  variable_container,
                                  grid, wic);
    well_blocks_ = well_spline_->GetWellBlocks();

  } else if (well_settings.definition_type == WType::PolarSpline) {
    if (vp_.vMOD >= 5) { ext_warn("POLARSPLINE: " + wn); }
    well_spline_ = new PolarSpline(well_settings,
                                   variable_container,
                                   grid, wic);
    well_blocks_ = well_spline_->GetWellBlocks();

  } else if (well_settings.definition_type == WType::PseudoContVertical2D) {
    if (vp_.vMOD >= 5) { ext_warn("VERT3D: " + wn); }
    pseudo_cont_vert_ = new PseudoContVert(well_settings,
                                           variable_container,
                                           grid);
    well_blocks_->append(pseudo_cont_vert_->GetWellBlock());
  }

  if (well_settings.definition_type != WType::UNDEFINED) {
    calculateDirectionOfPenetration();
  }

  if (vp_.vMOD >= 5) { printWellBlocks(); }
}

int Trajectory::GetTimeSpentInWic() const {
  if (well_spline_ != nullptr) {
    return well_spline_->GetTimeSpentInWIC();
  }
  else return 0;
}

WellBlock *Trajectory::GetWellBlock(int i, int j, int k) {
  for (int idx = 0; idx < well_blocks_->size(); ++idx) {
    if (well_blocks_->at(idx)->i() == i
      && well_blocks_->at(idx)->j() == j
      && well_blocks_->at(idx)->k() == k)
      return well_blocks_->at(idx);
  }
  throw WellBlockNotFoundException(i, j, k);
}

QList<WellBlock *> *Trajectory::GetWellBlocks() {
  return well_blocks_;
}

void Trajectory::UpdateWellBlocks() {
  // \todo This is the source of a memory leak:
  //  old well blocks are not deleted. Fix it.
  if (well_spline_ != nullptr) {
    if (well_spline_->HasGridChanged()
      || well_spline_->HasSplineChanged()) {
      well_blocks_ = well_spline_->GetWellBlocks();
    } else {
      if (vp_.vMOD >= 3) {
        im_ = "@[Trajectory::UpdateWellBlocks] ";
        im_ = "The well spline has not changed.";
        im_ += "well indices will not be recomputed.";
        info(im_, vp_.lnw);
      }
    }

  } else if (pseudo_cont_vert_ != 0) {
    well_blocks_ = new QList<WellBlock *>();
    well_blocks_->append(pseudo_cont_vert_->GetWellBlock());
  }
  calculateDirectionOfPenetration();
}

double Trajectory::GetLength() const {
  if (definition_type_ == WType::WellSpline) {
    return GetSplineLength();
  }
  else { // Block-defined
    string em = "Length of block-defined wells not yet implemented.";
    throw runtime_error(em);
  }
}

void Trajectory::
initializeWellBlocks(Settings::Model::Well well,
                     Properties::VarPropContainer *variable_container) {
  QList<Settings::Model::Well::WellBlock> blocks = well.well_blocks;
  for (int i = 0; i < blocks.size(); ++i) {
    well_blocks_->append(new WellBlock(blocks[i].i, blocks[i].j, blocks[i].k));
    if (blocks[i].is_variable) {
      well_blocks_->last()->i_->setName(blocks[i].name + "#i");
      well_blocks_->last()->j_->setName(blocks[i].name + "#j");
      well_blocks_->last()->k_->setName(blocks[i].name + "#k");
      variable_container->AddVariable(well_blocks_->last()->i_);
      variable_container->AddVariable(well_blocks_->last()->j_);
      variable_container->AddVariable(well_blocks_->last()->k_);
    }
    if (blocks[i].has_completion)
      well_blocks_->last()->AddCompletion(
        new Completions::Perforation(blocks[i].completion, variable_container));
  }
}

void Trajectory::calculateDirectionOfPenetration() {
  if (well_blocks_->size() == 1) { // Assuming that the well is vertical if it only has one block
    well_blocks_->first()->setDirOfPenetration(WellBlock::DirOfPenetration::Z);
    return;
  }
  // All but the last block use forward direction
  int delta_i, delta_j, delta_k;
  stringstream ss; ss << "@[Trajectory::calculateDirectionOfPenetration] delta ijk:" << endl;

  for (int i = 0; i < well_blocks_->size()-1; ++i) {
    delta_i = std::abs(well_blocks_->at(i)->i() - well_blocks_->at(i+1)->i());
    delta_j = std::abs(well_blocks_->at(i)->j() - well_blocks_->at(i+1)->j());
    delta_k = std::abs(well_blocks_->at(i)->k() - well_blocks_->at(i+1)->k());

    if (delta_i == 1 && delta_j == 0 && delta_k == 0) {
      well_blocks_->at(i)->setDirOfPenetration(WellBlock::DirOfPenetration::X);

    } else if (delta_i == 0 && delta_j == 1 && delta_k == 0) {
      well_blocks_->at(i)->setDirOfPenetration(WellBlock::DirOfPenetration::Y);

    } else if (delta_i == 0 && delta_j == 0 && delta_k == 1) {
      well_blocks_->at(i)->setDirOfPenetration(WellBlock::DirOfPenetration::Z);

    } else if (delta_i == 1 && delta_j == 1 && delta_k == 0) {
      well_blocks_->at(i)->setDirOfPenetration(WellBlock::DirOfPenetration::D);

    } else {
      well_blocks_->at(i)->setDirOfPenetration(WellBlock::DirOfPenetration::W);
    }
    if (vp_.vMOD >= 5) {
      ss << "delta_i = " << delta_i << "; delta_j = " << delta_j << "; delta_k = ";
      ss << delta_k << ";  ->  " << well_blocks_->at(i)->getDirPenetrationStr() << endl;
    }
  }

  // Last block uses backward direction
  delta_i = std::abs(well_blocks_->last()->i() - well_blocks_->at(well_blocks_->size()-2)->i());
  delta_j = std::abs(well_blocks_->last()->j() - well_blocks_->at(well_blocks_->size()-2)->j());
  delta_k = std::abs(well_blocks_->last()->k() - well_blocks_->at(well_blocks_->size()-2)->k());
  if (delta_i == 1 && delta_j == 0 && delta_k == 0) {
    well_blocks_->last()->setDirOfPenetration(WellBlock::DirOfPenetration::X);

  } else if (delta_i == 0 && delta_j == 1 && delta_k == 0) {
    well_blocks_->last()->setDirOfPenetration(WellBlock::DirOfPenetration::Y);

  } else if (delta_i == 0 && delta_j == 0 && delta_k == 1) {
    well_blocks_->last()->setDirOfPenetration(WellBlock::DirOfPenetration::Z);

  } else if (delta_i == 1 && delta_j == 1 && delta_k == 0) {
    well_blocks_->last()->setDirOfPenetration(WellBlock::DirOfPenetration::D);

  } else {
    well_blocks_->last()->setDirOfPenetration(WellBlock::DirOfPenetration::W);
  }
  if (vp_.vMOD >= 5) {
    ss << "delta_i = " << delta_i << "; delta_j = " << delta_j << "; delta_k = ";
    ss << delta_k << ";  ->  " << well_blocks_->last()->getDirPenetrationStr() << endl;
    cout << ss.str();
  }
}

WType Trajectory::GetDefinitionType() {
  return definition_type_;
}

void Trajectory::
convertWellBlocksToWellSpline(Settings::Model::Well &well_settings,
                              Reservoir::Grid::Grid *grid) {
  if (vp_.vMOD >= 2) {
    im_ = "Convering well " + well_settings.name.toStdString() + " to spline.";
    ext_info(im_, md_, cl_);
  }

  // std::cout << "Input blocks:" << std::endl;
  // for (auto imported_block : well_settings.well_blocks) {
  //     std::cout << imported_block.i << ",\t" << imported_block.j << ",\t" << imported_block.k << std::endl;
  // }

  QList<Settings::Model::Well::SplinePoint> points;
  for (int p = 0; p < well_settings.n_spline_points; ++p) {
    int srndg_block_idx = std::min(
      well_settings.well_blocks.size() - 1,
      int(p * std::ceil(well_settings.well_blocks.size() / (well_settings.n_spline_points-1)))
    );

    Reservoir::Grid::Cell srndg_cell;
    try {
      if (vp_.vMOD >=3) {
        cout << "Selected for spline conversion block nr " << srndg_block_idx << endl;
      }
      Settings::Model::Well::WellBlock srndg_block = well_settings.well_blocks[srndg_block_idx];
      srndg_cell = grid->GetCell(srndg_block.i-1, srndg_block.j-1, srndg_block.k-1);

    } catch (runtime_error &e) {
      wm_ = "Unable to get grid cell needed for spline conversion. Trying adjacent block.";
      ext_warn(wm_, md_, cl_);

      // std::cout << "WARNING: Unable to get grid cell needed for spline conversion: " << e.what()
      //           << " Trying adjacent block." << std::endl;
      // std::cout << "Selected for spline conversion block nr " << srndg_block_idx - 2 << std::endl;

      Settings::Model::Well::WellBlock srndg_block = well_settings.well_blocks[srndg_block_idx - 2];
      srndg_cell = grid->GetCell(srndg_block.i-1, srndg_block.j-1, srndg_block.k-1);
    }

    if (vp_.vMOD >=3) {
      std::cout << "Corresponding block: " << srndg_cell.ijk_index().i() + 1
                << ", " << srndg_cell.ijk_index().j() + 1
                << ", " << srndg_cell.ijk_index().k() + 1
                << " with center at " << srndg_cell.center().x()
                << ", " << srndg_cell.center().y()
                << ", " << srndg_cell.center().z() << std::endl;

    }

    Settings::Model::Well::SplinePoint new_point;
    new_point.name = "SplinePoint#" + well_settings.name + "#P" + QString::number(p+1);
    new_point.x = srndg_cell.center().x();
    new_point.y = srndg_cell.center().y();
    new_point.z = srndg_cell.center().z();

    if (well_settings.is_variable_spline) {
      new_point.is_variable = true;
    } else {
      new_point.is_variable = false;
    }

    points.push_back(new_point);
  }

  points.first().name = "SplinePoint#" + well_settings.name + "#heel";
  points.last().name = "SplinePoint#" + well_settings.name + "#toe";
  well_settings.spline_points = points;
  cout << "Done Converting well " << well_settings.name.toStdString() << " to spline." << endl;
}

WellBlock * Trajectory::GetWellBlockByMd(double md, string dmsg) {
  assert(well_blocks_->size() > 0);
  if (md > GetLength()) {
    E("Attempting to get well block at MD grater than well length.");
  }

  if (dmsg != "") {
    cout << dmsg << endl;
  }

  for (auto block : *well_blocks_) {
    if (vp_.vMOD >=3) { cout << block->getPropString(md) << endl; }
    if (md <= block->getExitMd() && md >= block->getEntryMd()) {
      if (vp_.vMOD >=3) {
        cout << "Returning block: [" << block->i() << " " << block->j() << " " << block->k() << "]" << endl;
      }
      return block;
    }
  }

  return dummy_block;

  // for (int ii=1; ii < well_blocks_->size()-1; ++ii) {
  //   auto curr_block = well_blocks_->at(ii);
  //   auto prev_block = well_blocks_->at(ii-1);
  //   if (md >= curr_block->getEntryMd() && md <= curr_block->getExitMd()) {
  //     return curr_block;
  //   } else if (md >= prev_block->getExitMd() && md <= curr_block->getEntryMd()) {
  //     if (vp_.vMOD >=3) {
  //       cout << "Returning block: ["
  //       << prev_block->i() << " " << prev_block->j() << " " << prev_block->k()
  //       << "] instead of block ["
  //       << curr_block->i() << " " << curr_block->j() << " " << curr_block->k()
  //       << "]" << endl;
  //     }
  //     return prev_block;
  //   }
  // }
  // E("Unable to get well block by MD.");
}

std::vector<WellBlock *> Trajectory::
GetWellBlocksByMdRange(double start_md, double end_md, string dmsg) const {
  std::vector<WellBlock *> affected_blocks;
  for (auto wb : *well_blocks_) {
    if (( start_md >= wb->getEntryMd() && start_md <= wb->getExitMd() ) || // Start md inside block
      (   end_md >= wb->getEntryMd() &&   end_md <= wb->getExitMd() ) || // End md inside block
      ( start_md <= wb->getEntryMd() &&   end_md >= wb->getExitMd() )) { // Block between start and end blocks
      affected_blocks.push_back(wb);
    }
  }
  return affected_blocks;
}

double Trajectory::GetEntryMd(const WellBlock *wb) const {
  return wb->getEntryMd();
}

double Trajectory::GetExitMd(const WellBlock *wb) const {
  return wb->getExitMd();
}

void Trajectory::printWellBlocks() {
  ext_info(wsettings_.toString("|"), md_, cl_, vp_.lnw);

  cout << endl << "well blocks -- Trajectory::printWellBlocks():" << endl;
  // cout << "I,\tJ,\tK,\tINX,\t\t\tINY,\t\t\tINZ,\t\t\tOUTX,\t\t\tOUTY,\t\t\tOUTZ" << endl;

  string sep = "         ";
  string lbl;
  lbl = "[    I    J    K ] [ Xin" + sep + "Yin" + sep + "Zin         ] ";
  lbl += "[ Xou" + sep + "You" + sep + "Zou" + sep + "] ";
  lbl +=  "[    Len EntryMd  ExitMd  DepthD ] ";
  lbl +=  "[      Dx      Dy     Dz ]";
  cout << lbl << endl;

  for (auto wb : *well_blocks_) {
    stringstream entry;
    entry << "[ "
          << num2str(wb->i(), 0, 0, 4) << " "
          << num2str(wb->j(), 0, 0, 4) << " "
          << num2str(wb->k(), 0, 0, 4) << " ] ";
    // entry.precision(3);
    // entry << std::scientific;

    entry << "[ "
          << num2str(wb->getEntryPoint().x(), 5, 1, 9) << " "
          << num2str(wb->getEntryPoint().y(), 5, 1, 9) << " "
          << num2str(wb->getEntryPoint().z(), 5, 1, 9) << " ] [ "
          << num2str(wb->getExitPoint().x(), 5, 1, 9) << " "
          << num2str(wb->getExitPoint().y(), 5, 1, 9) << " "
          << num2str(wb->getExitPoint().z(), 5, 1, 9) << " ] ";

    entry << "[ "
          << num2str(wb->getLength(), 1, 0, 6) << " "
          << num2str(wb->getEntryMd(), 1, 0, 7) << " "
          << num2str(wb->getExitMd(), 1, 0, 7) << " "
          << num2str(wb->getDepthChange(), 1, 0, 7) << " ] ";

    auto dxdydz = wb->getDxDyDz();
    if (dxdydz.size() != 0) {
      entry << "[ "
            << num2str(dxdydz(0), 1, 0, 7) << " "
            << num2str(dxdydz(1), 1, 0, 7) << " "
            << num2str(dxdydz(2), 1, 0, 6) << " ]"
            << std::endl;
    } else {
      entry << "[ "
            << num2str(0, 1, 0, 7) << " "
            << num2str(0, 1, 0, 7) << " "
            << num2str(0, 1, 0, 6) << " ]"
            << std::endl;
    }
    cout << entry.str();
  }
}

double Trajectory::GetSplineLength() const {
  return well_blocks_->back()->getExitMd();
}



}
}
}
