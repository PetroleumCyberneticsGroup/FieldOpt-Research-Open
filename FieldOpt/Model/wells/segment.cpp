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
#include "segment.h"
#include "Utilities/stringhelpers.hpp"

namespace Model {
namespace Wells {

Model::Wells::Segment::Segment() { }
Segment::Segment(Segment::SegType type,
                 int index,
                 int branch,
                 int outlet,
                 double length,
                 double tvd,
                 double diameter,
                 double roughness,
                 double outlet_md) {
  type_ = type;
  index_ = index;
  branch_ = branch;
  outlet_ = outlet;
  length_ = length;
  tvd_change_ = tvd;
  diameter_ = diameter;
  roughness_ = roughness;
  md_ = outlet_md + length;
}
void Segment::AddInlet(int index) {
  inlets_.push_back(index);
}
std::vector<int> Segment::GetInlets() const {
  return inlets_;
}
void Segment::AddParentBlock(Wellbore::WellBlock *parent_block) {
  parent_block_ = parent_block;
}
void Segment::AddParentICD(Wellbore::Completions::ICD *parent_icd) {
  parent_icd_ = parent_icd;
}
bool Segment::HasParentBlock() const {
  return !parent_block_ == 0;
}
bool Segment::HasParentICD() const {
  return !parent_icd_ == 0;
}
std::string Segment::ToString() {
  std::stringstream s;
  s << "-- Segment " << index_ << " (";
  if (type_ == TUBING_SEGMENT)  s << "tubing) --\n";
  if (type_ == ICD_SEGMENT)     s << "ICD) --\n";
  if (type_ == ANNULUS_SEGMENT) s << "annulus) --\n";
  s << "   branch: " << branch_ << "\n";
  s << "   outlet: " << outlet_ << "\n";
  s << "   inlet:  " << vec_to_str(inlets_) << "\n";
  s << "   length: " << length_ << "\n";
  s << "   tvd:    " << tvd_change_ << "\n";
  s << "   diam:   " << diameter_ << "\n";
  s << "   roughness: " << roughness_ << "\n";
  if (type_ == ANNULUS_SEGMENT) {
    s << "   block i, j, k: " << parent_block_->i() << ", " << parent_block_->j() << ", " << parent_block_->k() << "\n";
  }
  return s.str();
}

}
}
