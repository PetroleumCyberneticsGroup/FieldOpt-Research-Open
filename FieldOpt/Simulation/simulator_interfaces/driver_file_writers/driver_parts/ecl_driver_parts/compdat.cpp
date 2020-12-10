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

#include "compdat.h"
#include <iostream>
#include <utility>

namespace Simulation {
namespace ECLDriverParts {

using Printer::ext_info;
using Printer::num2strQ;

Compdat::Compdat(QList<Model::Wells::Well *> *wells) {
  initializeBaseEntryLine(13);
  head_ = "COMPDAT";
  foot_ = "/\n\n";
  for (int i = 0; i < wells->size(); ++i) {
    if (wells->at(i)->trajectory()->GetWellBlocks()->size() > 0)
      entries_.append(createWellEntries(wells->at(i)));
  }
}

Compdat::Compdat(QList<Model::Wells::Well *> *wells, int timestep) {
  initializeBaseEntryLine(13);
  head_ = "COMPDAT";
  foot_ = "/\n\n";
  for (auto well : *wells) {
    if (well->controls()->first()->time_step() == timestep) {
      entries_.append(createWellEntries(well));
    }
  }
}

QString Compdat::GetPartString() const {
  // Return empty string if there are no entries (at this timestep)
  if (entries_.size() == 0) {
    return "";
  }

  QString entries = head_ + "\n";
  for (QStringList entry : entries_) {
    entries.append("  " + entry.join(" ") + " /\n");
  }
  entries.append("" + foot_);
  return entries;
}

QList<QStringList> Compdat::createWellEntries(Model::Wells::Well *well) {
  QList<QStringList> block_entries = QList<QStringList>();
  for (int i = 0; i < well->trajectory()->GetWellBlocks()->size(); ++i) {
    block_entries.append(createBlockEntry(well->name(),
                                          well->wellbore_radius()*2.0,
                                          well->trajectory()->GetWellBlocks()->at(i)));
  }
  return block_entries;
}

QStringList Compdat::createBlockEntry(QString well_name,
                                      double wellbore_radius,
                                      Model::Wells::Wellbore::WellBlock *well_block) {
  QStringList block_entry = QStringList(base_entry_line_);
  block_entry[0] = well_name;
  block_entry[1] = num2strQ(well_block->i(), 0, 0, 4) + "";
  block_entry[2] = num2strQ(well_block->j(), 0, 0, 4) + "";
  block_entry[3] = num2strQ(well_block->k(), 0, 0, 4) + "";
  block_entry[4] = num2strQ(well_block->k(), 0, 0, 4) + " ";

  if (well_block->HasPerforation()) {
    block_entry[5] = "OPEN";
    if (well_block->GetPerforation()->transmissibility_factor() >= 0.0)
      block_entry[7] = " " + num2strQ(well_block->GetPerforation()->transmissibility_factor(), 5, 0, 9) + "";
    // WI is defaulted if a negative value is provided.
  }
  else {
    block_entry[5] = "SHUT";
  }

  block_entry[8] = num2strQ(wellbore_radius, 4, 0, 7) + "";

  // \note By default, this is set to X. This is primarily
  // relevant if the well passes diagonally through a block.
  im_ = "@[Compdat::createBlockEntry] ";
  im_ += "Either Trajectory::calculateDirectionOfPenetration failed ";
  im_ += "or penetration dir for well: " + well_name.toStdString() + " not def. |";

  block_entry[12] = "X";
  switch (well_block->directionOfPenetration()) {
    case Model::Wells::Wellbore::WellBlock::DirOfPenetration::X:
      block_entry[12] = "X";
      break;
    case Model::Wells::Wellbore::WellBlock::DirOfPenetration::Y:
      block_entry[12] = "Y";
      break;
    case Model::Wells::Wellbore::WellBlock::DirOfPenetration::Z:
      block_entry[12] = "Z";
      break;
    case Model::Wells::Wellbore::WellBlock::DirOfPenetration::D:
      block_entry[12] = "X";
      im_ += " Well is exactly diagonal. Defaulting to \"X\"";
      im_ += "[mod: " + md_ +  "] [cls: " + cl_ + "]:\n" + im_;
      cout << im_ << endl;
      break;
    case Model::Wells::Wellbore::WellBlock::DirOfPenetration::W:
      block_entry[12] = "X";
      im_ += "DirOfPenetration::W -> overriding to X for ";
      im_ += well_block->getPropString();
      ext_info(im_, md_, cl_);
      break;
      // throw std::runtime_error(im_);
    default:
      im_ += "No DirOfPenetration defined.";
      throw std::runtime_error(im_);
  }
  return block_entry;
}

}
}
