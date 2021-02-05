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

#include <Utilities/printer.hpp>
#include <boost/lexical_cast.hpp>
#include "Utilities/verbosity.h"
#include "well.h"

namespace Model {
namespace Wells {

using Printer::info;
using Printer::ext_info;
using Printer::ext_warn;
using Printer::num2str;

using WType=Settings::Model::WellDefinitionType;

Well::Well(const Settings::Model& model_settings,
           int well_number,
           Properties::VarPropContainer *variable_container,
           Reservoir::Grid::Grid *grid,
           Reservoir::WellIndexCalculation::wicalc_rixx *wic) {

  well_settings_ = model_settings.wells()[well_number];
  well_settings_.copyVerbParams(model_settings.vp_);
  vp_ = well_settings_.verbParams();

  name_ = well_settings_.name;
  type_ = well_settings_.type;

  if (well_settings_.group.length() >= 1) {
    group_ = well_settings_.group;
  } else {
    group_ = "";
  }

  preferred_phase_ = well_settings_.preferred_phase;

  wellbore_radius_ = new Properties::ContinuousProperty(well_settings_.wellbore_radius);

  controls_ = new QList<Control *>();
  for (const auto & control : well_settings_.controls) {
    controls_->append(new Control(control,well_settings_,
                                  variable_container));
  }

  trajectory_ = new Wellbore::Trajectory(well_settings_,
                                         variable_container,
                                         grid, wic);
  well_blocks_ = trajectory_->GetWellBlocks();

  trajectory_defined_ = well_settings_.definition_type != WType::UNDEFINED;
  if (trajectory_defined_) {
    heel_.i = trajectory_->GetWellBlocks()->first()->i();
    heel_.j = trajectory_->GetWellBlocks()->first()->j();
    heel_.k = trajectory_->GetWellBlocks()->first()->k();

    // if (well_settings_.wseg_structure == "ICDBranches") {
    //   initSegWellStruct(variable_container);
    //   cout << "well_settings_.wseg_structure: " << well_settings_.wseg_structure << endl;
    //   cout << "well_settings_.use_segmented_model: " << well_settings_.use_segmented_model << endl;
    // } else

    if (well_settings_.use_segmented_model) {
      is_segmented_ = true;
      // initializeSegmentedWell(variable_container);
      initSegWellStruct(variable_container);
    }

  } else {
    im_ = "No well trajectory defined for well: " + name_.toStdString();
    ext_info(im_, md_, cl_, vp_.lnw);
  }

  // Completions for wells with no defined trajectory
  // (they are specified by segment number)
  if (!well_settings_.completions.empty()
    && well_settings_.icv_compartments.empty()) {

    for (auto comp : well_settings_.completions) {
      auto base_name = comp.name;
      auto comp_num = comp.icd_names.size();
      for (int i = 0; i < comp_num; ++i) {
        comp.icd_name = comp.icd_names[i];
        comp.icd_segment = comp.icd_segments[i];
        comp.name = base_name + "#" + QString::number(comp.icd_segment);
        icds_.emplace_back(comp, variable_container);
      }
      if (vp_.vMOD >= 2) {
        im_ = "well_settings_.completions => [NOT-EMPTY] ";
        im_ += "well_settings_.icv_compartments => [EMPTY] -- ";
        im_ += "Added " + num2str(comp_num) + " ICDs.";
        ext_info(im_, md_, cl_, vp_.lnw);
      }
    }

  } else if (!well_settings_.icv_compartments.empty()) {
    for (const auto& comp : well_settings_.icv_compartments) {
      icds_.emplace_back(comp, variable_container);
    }
    if (vp_.vMOD >= 2) {
      im_ = "well_settings_.icv_compartments => [NOT-EMPTY] -- ";
      im_ += "Added " + num2str(well_settings_.icv_compartments.size()) + " ICV compartment.";
      ext_info(im_, md_, cl_, vp_.lnw);
    }
  }

}

bool Well::IsProducer() {
  return type_ == ::Settings::Model::WellType::Producer;
}

bool Well::IsInjector() {
  return type_ == ::Settings::Model::WellType::Injector;
}

void Well::Update() {
  if (trajectory_defined_) {
    trajectory_->UpdateWellBlocks();
    heel_.i = trajectory_->GetWellBlocks()->first()->i();
    heel_.j = trajectory_->GetWellBlocks()->first()->j();
    heel_.k = trajectory_->GetWellBlocks()->first()->k();

    if (is_segmented_) {
      for (auto compartment : compartments_) {
        compartment.Update();
      }
    }
  }
}

// https://stackoverflow.com/questions/6861089/how-to-split-a-vector-into-n-almost-equal-parts
template<typename T>
std::vector<std::vector<T>> Well::SplitVector(const std::vector<T>& vec, size_t n, bool dbg) {
  std::vector<std::vector<T>> outVec;

  size_t length = vec.size() / n;
  size_t remain = vec.size() % n;

  size_t begin = 0;
  size_t end = 0;

  for (size_t i = 0; i < std::min(n, vec.size()); ++i) {
    end += (remain > 0) ? (length + ((remain--) != 0)) : length;

    outVec.push_back(std::vector<T>(vec.begin() + begin, vec.begin() + end));

    begin = end;
  }

  // dbg
  if (dbg) {
    int sz_chunk;
    for(auto nn: vec) { cout << nn << " "; }
    cout << endl;

    sz_chunk=0;
    for(const auto& idx_chunk: outVec) {
      for (auto nn: idx_chunk) {
        cout << num2str(nn, 0, 0, 2) << " ";
      }
      cout << endl;
      sz_chunk += idx_chunk.size();
    }
    cout << "sz_chunk:" << sz_chunk << endl;
  }

  return outVec;
}


void Well::initSegWellStruct(Properties::VarPropContainer *variable_container) {

  // Basis: 1 segment-per-well-block, i.e., nsegs = nblocks
  // General: Orders segments into branch structure,
  // defines ncomp compartments
  // ncomp = 8 compartments w/ 1 ICD each
  // # of branches? 1 + 8 + 8 -> nbrnch = 2*ncomp + 1

  // -------------------------------------------------------
  assert(well_settings_.seg_n_compartments > 0);
  assert(trajectory_->GetDefinitionType() == Settings::Model::WellDefinitionType::WellSpline);

  tub_diam_ = well_settings_.seg_tubing.diameter;
  tub_cross_sect_area_ = well_settings_.seg_tubing.cross_sect_area;
  tub_roughness_ = well_settings_.seg_tubing.roughness;

  icd_diam_ = well_settings_.seg_compartment_params.diameter;
  icd_cross_sect_area_ = well_settings_.seg_compartment_params.cross_sect_area;
  icd_roughness_ = well_settings_.seg_compartment_params.roughness;

  ann_diam_ = well_settings_.seg_annulus.diameter;
  ann_cross_sect_area_ = well_settings_.seg_annulus.cross_sect_area;
  ann_roughness_ = well_settings_.seg_annulus.roughness;

  std::vector<Segment> segments;

  // -------------------------------------------------------
  auto wblocks = trajectory()->GetWellBlocks();
  auto ncomps = well_settings_.seg_n_compartments;

  std::vector<int> wblock_idx;
  std::vector<std::vector<int>> wblock_idx_chunks;

  wblock_idx.resize(wblocks->size());
  std::iota(wblock_idx.begin(), wblock_idx.end(), 0);

  wblock_idx_chunks = SplitVector(wblock_idx, ncomps, false);

  // -------------------------------------------------------
  // int first_idx = 0;
  std::vector<Wellbore::WellBlock *> block_range;
  std::vector<int> block_idx_range, num_vec;
  // int curr_range,
  int total_range = 0;

  // -------------------------------------------------------
  // Loop through compartments
  for (int i = 0; i < well_settings_.seg_n_compartments; ++i) {

    block_range.clear();
    for (auto range_idx: wblock_idx_chunks[i]) {
      block_range.push_back(wblocks->at(range_idx));
      block_idx_range.push_back(range_idx);
    }
    total_range += block_range.size(); // <- for dbg check

    auto first_block = wblocks->at(wblock_idx_chunks[i].front());
    auto last_block = wblocks->at(wblock_idx_chunks[i].back());

    // Make compartments
    if (first_block->i() == last_block->i() &&
      first_block->j() == last_block->j() &&
      first_block->k() == last_block->k()) {

      auto single_block = first_block;
      std::vector<Wellbore::WellBlock *> single_range;
      single_range.push_back(single_block);

      compartments_.emplace_back(single_block->getEntryMd(),
                                 single_block->getExitMd(),
                                 single_block->getEntryPoint().z(),
                                 single_block->getExitPoint().z(),
                                 trajectory()->GetLength(),
                                 single_range,
                                 block_idx_range,
                                 well_settings_,
                                 variable_container,
                                 compartments_);

    } else { // if (i < well_settings_.seg_n_compartments - 1) {

      compartments_.emplace_back(first_block->getEntryMd(),
                                 last_block->getExitMd(),
                                 first_block->getEntryPoint().z(),
                                 last_block->getExitPoint().z(),
                                 trajectory()->GetLength(),
                                 block_range,
                                 block_idx_range,
                                 well_settings_,
                                 variable_container,
                                 compartments_);

    }
  }

  if (total_range != wblocks->size()) {
    wm_ = "Sum of compartment blocks (" + num2str(total_range) +
      ") not equal to original # of blocks (" + num2str(wblocks->size()) + ").";
    ext_warn(wm_, md_, cl_);
  } else if (block_idx_range.back() + 1 != wblocks->size()) {
    wm_ = "Last block number (" + num2str(block_idx_range.back()) +
      ") not equal to total # of well blocks (" + num2str(wblocks->size()) + ").";
    ext_warn(wm_, md_, cl_);
  }
}

void Well::initializeSegmentedWell(Properties::VarPropContainer *variable_container) {
  assert(well_settings_.seg_n_compartments > 0);
  assert(trajectory_->GetDefinitionType() == Settings::Model::WellDefinitionType::WellSpline);

  tub_diam_ = well_settings_.seg_tubing.diameter;
  tub_cross_sect_area_ = well_settings_.seg_tubing.cross_sect_area;
  tub_roughness_ = well_settings_.seg_tubing.roughness;

  icd_diam_ = well_settings_.seg_compartment_params.diameter;
  icd_cross_sect_area_ = well_settings_.seg_compartment_params.cross_sect_area;
  icd_roughness_ = well_settings_.seg_compartment_params.roughness;

  ann_diam_ = well_settings_.seg_annulus.diameter;
  ann_cross_sect_area_ = well_settings_.seg_annulus.cross_sect_area;
  ann_roughness_ = well_settings_.seg_annulus.roughness;

  double compartment_length = trajectory_->GetLength() / well_settings_.seg_n_compartments;
  std::vector<double> compartment_delimiters;

  // Create approximate compartment delimiters
  compartment_delimiters.push_back(0.0);
  for (int i = 1; i < well_settings_.seg_n_compartments; ++i) {
    compartment_delimiters.push_back(i * compartment_length);
  }
  compartment_delimiters.push_back(trajectory_->GetLength());

  // Snap compartment delimiters to wellblock-intersections
  for (int l = 1; l < compartment_delimiters.size(); ++l) {
    auto wb = trajectory_->GetWellBlockByMd(compartment_delimiters[l],
                                            "Snap compartment delimiter to block:");
    while (wb->i()==0 && wb->j()==0 && wb->k()==0) {
      compartment_delimiters[l] = compartment_delimiters[l]*1.025;
      wb = trajectory_->GetWellBlockByMd(compartment_delimiters[l],
                                         "Retrying snap of compartment delimiter to block:");
    }

    compartment_delimiters[l] = wb->getExitMd() - 0.1;
    assert(compartment_delimiters[l] != compartment_delimiters[l-1]);
  }
  assert(compartment_delimiters.size() == well_settings_.seg_n_compartments + 1);

  // get trajectory wellblocks
  // auto well_blocks = trajectory_->GetWellBlocks();
  auto well_blocks = well_blocks_;
  int first_idx = 0;
  assert(compartment_delimiters.size() < trajectory_->GetWellBlocks()->size());
  std::vector<Wellbore::WellBlock *> block_range;
  std::vector<int> block_idx_range, num_vec;
  int curr_range, total_range = 0;

  // Loop through compartments
  for (int i = 0; i < well_settings_.seg_n_compartments; ++i) {

    auto first_block = well_blocks->at(first_idx);
    auto last_block = trajectory_->GetWellBlockByMd(compartment_delimiters[i+1],
                                                    "Finding last block in compartment:");

    // Get compartment blocks
    block_range = trajectory_->GetWellBlocksByMdRange(
      compartment_delimiters[i] - 0.1,
      compartment_delimiters[i+1] + 0.1);

    // Fix overlapping blocks in compartment selections
    if (i == 0) {
      block_range.erase(block_range.end()-1);
    } else if (i == well_settings_.seg_n_compartments - 1) {
      block_range.erase(block_range.begin());
    } else {
      block_range.erase(block_range.begin());
      block_range.erase(block_range.end()-1);
    }

    // Make increasing count of blocks to tag segments in each compartment
    num_vec.resize(block_range.size());
    if  (i > 0) { curr_range = block_idx_range.back() + 1; } else { curr_range = 0; }
    std::iota(num_vec.begin(), num_vec.end(), curr_range);
    block_idx_range = num_vec;
    total_range += block_idx_range.size(); // dbg check

    // Make compartments
    if (first_block->i() == last_block->i() &&
      first_block->j() == last_block->j() &&
      first_block->k() == last_block->k()) { // Compartment begins and ends in the same block

      // ext_warn("Compartment " + num2str(i) + " begins and ends in the same well block.", md_, cl_);
      auto single_block = first_block;
      std::vector<Wellbore::WellBlock *> single_range;
      single_range.push_back(single_block);

      compartments_.emplace_back(single_block->getEntryMd(),
                                 single_block->getExitMd(),
                                 single_block->getEntryPoint().z(),
                                 single_block->getExitPoint().z(),
                                 trajectory()->GetLength(),
                                 single_range,
                                 block_idx_range,
                                 well_settings_, variable_container, compartments_);
      first_idx++;

    } else if (i < well_settings_.seg_n_compartments - 1) { // Not last compartment

      compartments_.emplace_back(first_block->getEntryMd(),
                                 last_block->getExitMd(),
                                 first_block->getEntryPoint().z(),
                                 last_block->getExitPoint().z(),
                                 trajectory()->GetLength(),
                                 block_range,
                                 block_idx_range,
                                 well_settings_, variable_container, compartments_);

      // Move to after last_block
      while (! (well_blocks->at(first_idx)->i() == last_block->i()
        && well_blocks->at(first_idx)->j() == last_block->j()
        && well_blocks->at(first_idx)->k() == last_block->k())) {
        first_idx++;
      }
      first_idx++;

    } else { // Last compartment

      compartments_.emplace_back(first_block->getEntryMd(),
                                 last_block->getExitMd(),
                                 first_block->getEntryPoint().z(),
                                 last_block->getExitPoint().z(),
                                 trajectory()->GetLength(),
                                 block_range,
                                 block_idx_range,
                                 well_settings_, variable_container, compartments_);
    }
  }

  if (total_range != well_blocks->size()) {
    wm_ = "Sum of compartment blocks (" + num2str(total_range) +
      ") not equal to original # of blocks (" + num2str(well_blocks->size()) + ").";
    ext_warn(wm_, md_, cl_);
  } else if (block_idx_range.back() + 1 != well_blocks->size()) {
    wm_ = "Last block number (" + num2str(block_idx_range.back()) +
      ") not equal to total # of well blocks (" + num2str(well_blocks->size()) + ").";
    ext_warn(wm_, md_, cl_);
  }
}

// Segmented well methods
std::vector<Packer *> Well::GetPackers() const {
  assert(is_segmented_);
  std::vector<Packer *> packers;
  packers.push_back(compartments_[0].start_packer);
  for (auto comp : compartments_) {
    packers.push_back(comp.end_packer);
  }
  return packers;
}

std::vector<ICD *> Well::GetICDs() const {
  assert(is_segmented_);
  std::vector<ICD *> icds;
  for (auto comp : compartments_) {
    icds.push_back(comp.icd);
  }
  return icds;
}

std::vector<Compartment> Well::GetCompartments() const {
  assert(is_segmented_);
  return compartments_;
}

std::vector<Segment> Well::GetSegments() {
  assert(is_segmented_);

  std::vector<Segment> segments;
  auto root_segment = Segment(Segment::SegType::TUBING_SEGMENT,
                              1, 1, -1, -1.0, -1.0,
                              compartments_.front().start_packer->tvd(), -1.0,
                              -1.0, 0.0, 0.0, zero_pt_, zero_pt_,
                              1, 1, -1);

  root_segment.AddParentBlock(well_blocks_->at(0));
  segments.push_back(root_segment);

  std::vector<int> tubing_indexes = createTubingSegments(segments);
  std::vector<int> icd_indexes = createICDSegments(segments, tubing_indexes);
  createAnnulusSegments(segments, icd_indexes);
  return segments;
}

void Well::createAnnulusSegments(std::vector<Segment> &segments,
                                 const std::vector<int> &icd_indexes) {
  assert(is_segmented_);
  std::vector<int> annulus_indexes;

  int tot_comps = 0;
  int aux_idx;

  for (int i = 0; i < compartments_.size(); ++i) {

    // Get compartment blocks (already executed in initializeSegmentedWell)
    auto comp_blocks = trajectory()->GetWellBlocksByMdRange(
      compartments_[i].start_packer->md(trajectory()->GetLength()) + 0.1,
      compartments_[i].end_packer->md(trajectory()->GetLength() - 0.1));

    // dbg checks
    // for (int k=0; k < comp_blocks.size(); k++) {
    //   // assert(comp_blocks[i]->i() == compartments_[i].block_range_[k]->i());
    //   cout << "comp_blocks[k]->i(): " << comp_blocks[k]->i()
    //   << " -- compartments_[i].block_range_[k]->i()): " << compartments_[i].block_range_[k]->i() << endl;
    // }

    // tot_comps += comp_blocks.size();
    // cout << "# of blocks in compartment [i:" << i << "] = " << comp_blocks.size() << " total: " << tot_comps << endl;

    // Should be compartment.block range, but then we need to update well_blocks_
    // with the compartment.block range found in initializeSegmentedWell()

    double outlet_md = -1;

    for (int j = 0; j < comp_blocks.size() - 2; ++j) {
      // for (int j = 0; j < compartments_[i].block_range_.size() - 2; ++j) {

      int index = 2*compartments_.size() + annulus_indexes.size() + 2;
      int outlet_index;
      // int block_index = 2*compartments_.size() + annulus_indexes.size() + 2;

      if (j == 0) { // Outlet to ICD
        outlet_index = icd_indexes[i];
        for (int k = 0; k < segments.size(); ++k) {
          if (segments[k].Type() == Segment::ICD_SEGMENT
            && segments[k].Index() == icd_indexes[i]) {
            // outlet_md = segments[k].OutletMD();
            segments[k].AddInlet(index);
            break;
          }
        }

        outlet_md = compartments_[i].block_range_.at(0)->getEntryMd() + length_delta_; // 0 block (tubing)
        outlet_md += 0.08, // 1 block (ICD)
        outlet_md += comp_blocks[j]->getLength();

        assert(outlet_md != -1); // Ensure that the outlet md has been set
        aux_idx = outlet_index - 1;

      } else { // Outlet to previous annulus segment
        outlet_index = segments.back().Index();
        // outlet_md = segments.back().OutletMD();
        segments.back().AddInlet(index);
        aux_idx = -1;

        outlet_md += comp_blocks[j]->getLength();
        if (j == comp_blocks.size() - 3) {
          outlet_md += comp_blocks[j]->getLength(); // next-to last block
          outlet_md += comp_blocks[j]->getLength(); // last block
        }

      }

      // cout << "compartments_[i].block_idx_range_[j+2]: " << compartments_[i].block_idx_range_[j+2] << endl;

      auto ann_seg = Segment(
        Segment::ANNULUS_SEGMENT,
        index,                         // index
        compartments_.size() + 2 + i,  // branch
        outlet_index,                  // outlet
        //
        // compartments_[i].block_range_[j+2]->getLength(),
        // compartments_[i].block_range_[j+2]->getDepthChange(),
        // comp_blocks[j]->getLength(),
        // comp_blocks[j]->getDepthChange(),
        comp_blocks[j]->getLength(),
        outlet_md,
        comp_blocks[j]->getDepth(),
        // (comp_blocks[j]->getExitPoint() - comp_blocks[j]->getEntryPoint()).norm(), // length
        // comp_blocks[j]->getExitPoint().z() - comp_blocks[j]->getEntryPoint().z(),  // tvd delta
        //
        0.0889,
        // ann_diam_,
        ann_roughness_,
        //
        // compartments_[i].block_range_[j+2]->getLength()/2,
        // outlet_md,
        // comp_blocks[j]->getLength()/2,
        comp_blocks[j]->getEntryMd() + length_delta_,
        comp_blocks[j]->getExitMd() + length_delta_,
        comp_blocks[j]->getEntryPoint(),
        comp_blocks[j]->getExitPoint(),
        //
        compartments_[i].block_idx_range_[j+2],
        compartments_.size() + 2 + i,
        aux_idx
      );

      ann_seg.AddParentBlock(comp_blocks[j+2]);
      annulus_indexes.push_back(index);
      segments.push_back(ann_seg);
    }
  }
}

std::vector<int> Well::createICDSegments(std::vector<Segment> &segments,
                                         std::vector<int> &tubing_indexes) const {
  assert(is_segmented_);
  std::vector<int> icd_indexes;
  int outlet = 1;
  stringstream ss;

  double seg_length = 0.0;
  double md_length = 0.0;
  double depth = 0.0;

  for (int i = 0; i < compartments_.size(); ++i) {
    // int index = tubing_indexes.back() + 1 + i;
    int index = i + 2;
    // if (i > 0) {
      outlet = index + compartments_.size();
    // }

    if (vp_.vMOD >=5) {
      ss << "@[Well::createICDSegments] " + compartments_[i].block_range_[1]->getPropString() << "|";
    }

    seg_length = compartments_[i].block_range_.front()->getExitMd() - compartments_[i].block_range_.front()->getEntryMd();
    md_length = compartments_[i].block_range_.front()->getEntryMd() + length_delta_;
    depth = compartments_[i].block_range_.front()->getEntryPoint().z();

    auto icd_segment = Segment(
      Segment::ICD_SEGMENT,
      index, // index
      i + 2, // branch
      outlet, // outlet
      seg_length + 0.08,
      md_length + 0.08,
      depth,
      // compartments_[i].block_range_[1]->getLength(),
      // compartments_[i].block_range_[1]->getDepthChange(),
      // well_blocks_->at(index)->getLength(),
      // well_blocks_->at(index)->getDepthChange(),
      // 0.1,   // length
      // 0.0,   // tvd delta
      icd_diam_,
      icd_roughness_,
      //
      // segments[tubing_indexes[i]].OutletMD(), // outlet md
      // compartments_[i].block_range_[1]->getLength()/2,
      compartments_[i].block_range_[1]->getEntryMd() + length_delta_,
      compartments_[i].block_range_[1]->getExitMd() + length_delta_,
      compartments_[i].block_range_[1]->getEntryPoint(),
      compartments_[i].block_range_[1]->getExitPoint(),
      //
      compartments_[i].block_idx_range_[1],
      index + compartments_.size(),
      -1
    );
    icd_segment.AddParentBlock(compartments_[i].block_range_[1]);

    icd_segment.AddParentICD(compartments_[i].icd);
    icd_indexes.push_back(index);
    segments[i].AddInlet(index);

    segments.push_back(icd_segment);
  }

  if (vp_.vMOD >=5) {
    ext_info(ss.str(), md_, cl_);
  }
  return icd_indexes;
}

std::vector<int> Well::createTubingSegments(std::vector<Segment> &segments) {
  assert(is_segmented_);
  std::vector<int> tubing_indexes;
  int outlet = 1;
  int seg_index;

  // start-dbg
  // int test_idx = 0;
  // std::vector<int> test_idx_vec;
  //
  // for (int i = 0; i < compartments_.size(); ++i) {
  //   seg_index = i + 2 + compartments_.size();
  //
  //   cout << "seg_idx: " << seg_index << "; block_idx: " << i << " -> block.I: " << well_blocks_->at(i)->i();
  //   cout << "; blck_range.sz: " << compartments_[i].block_range_.size() << " -- ";
  //   for (auto block : compartments_[i].block_range_) {
  //     cout << block->i() << " ";
  //   }
  //   cout << " | ";
  //
  //   for (int num : compartments_[i].block_idx_range_) {
  //     cout << num << " ";
  //   }
  //   cout << " |" << endl;
  // }
  //
  // // Build seg_idx <-> block_idx table
  // for (int j = 0; j < well_blocks_->size() - 2; ++j) {
  //   test_idx = 2 * compartments_.size() + test_idx_vec.size() + 2;
  //   test_idx_vec.push_back(test_idx);
  //   cout << "test_idx: " << test_idx << "; test_idx_vec.sz: " << test_idx_vec.size() << endl;
  // }
  // end-dbg

  double seg_length = 0.0;
  double md_length = 0.0;
  double depth = 0.0;

  for (int i = 0; i < compartments_.size(); ++i) {
    seg_index = i + 2 + compartments_.size();

    if (i > 0) { outlet = seg_index - 1; }

    seg_length = compartments_[i].block_range_.front()->getExitMd() - compartments_[i].block_range_.front()->getEntryMd();
    md_length = compartments_[i].block_range_.front()->getEntryMd() + length_delta_;
    depth = compartments_[i].block_range_.front()->getEntryPoint().z();

    // if (i==0) {
    //
    // } else {
    //   length = compartments_[i].block_range_.front()->getEntryMd();
    // }

    auto tubing_segment = Segment(
      Segment::TUBING_SEGMENT, // type
      seg_index,  // seg_index
      1,  // branch
      outlet, // outlet
      //
      seg_length,
      md_length,
      depth,
    // compartments_[i].block_range_.back()->getExitMd() - compartments_[i].block_range_.front()->getEntryMd(),
    //   compartments_[i].block_range_.back()->getExitPoint().z() - compartments_[i].block_range_.front()->getEntryPoint().z(),
      //
      // compartments_[i].block_range_[0]->getLength(),
      // compartments_[i].block_range_[0]->getDepthChange(),
      //
      // well_blocks_->at(block_index)->getLength(),
      // well_blocks_->at(block_index)->getDepthChange(),
      // compartments_[i].GetLength(trajectory_->GetLength()), // length
      // compartments_[i].GetTVDDifference(),  // tvd delta

      tub_diam_,
      tub_roughness_,
      //
      // segments.back().OutletMD(),
      // compartments_[i].block_range_[0]->getLength()/2,
      compartments_[i].block_range_[0]->getEntryMd() + length_delta_,
      compartments_[i].block_range_[0]->getExitMd() + length_delta_,
      compartments_[i].block_range_[0]->getEntryPoint(),
      compartments_[i].block_range_[0]->getExitPoint(),
      //
      compartments_[i].block_idx_range_[0], // block_index,
      seg_index,
      -1
    );
    // tubing_segment.AddParentBlock(well_blocks_->at(block_index));
    tubing_segment.AddParentBlock(compartments_[i].block_range_[0]);

    segments.push_back(tubing_segment);
    tubing_indexes.push_back(seg_index);
    segments[segments.size() - 2].AddInlet(seg_index);
  }
  return tubing_indexes;
}

std::vector<Segment> Well::GetTubingSegments() {
  assert(is_segmented_);
  std::vector<Segment> segments;
  auto root_segment = Segment(Segment::SegType::TUBING_SEGMENT,
    1, 1, -1, -1.0, -1.0,
    compartments_.front().start_packer->tvd(),-1.0,
    -1.0, 0.0, 0.0, zero_pt_, zero_pt_,
    1, 1, -1);
  segments.push_back(root_segment);
  std::vector<int> tubing_indexes = createTubingSegments(segments);
  return segments;
}

std::vector<Segment> Well::GetICDSegments() {
  assert(is_segmented_);
  std::vector<Segment> segments;

  auto root_segment = Segment(Segment::SegType::TUBING_SEGMENT,
    1, 1,-1, -1.0, -1.0,
    compartments_.front().start_packer->tvd(),-1.0,
    -1.0, 0.0, 0.0, zero_pt_, zero_pt_,
    1, 1, -1);
  segments.push_back(root_segment);

  std::vector<int> tubing_indexes = createTubingSegments(segments);
  std::vector<int> icd_indexes = createICDSegments(segments, tubing_indexes);
  std::vector<Segment> icd_segments;
  for (auto seg : segments) {
    if (seg.Type() == Segment::SegType::ICD_SEGMENT) {
      icd_segments.push_back(seg);
    }
  }
  return icd_segments;
}

std::vector<Segment> Well::GetAnnulusSegments() {
  assert(is_segmented_);
  std::vector<Segment> segments;
  auto root_segment = Segment(Segment::SegType::TUBING_SEGMENT,
    1, 1, -1, -1.0, -1.0,
    compartments_.front().start_packer->tvd(),-1.0,
    -1.0, 0.0, 0.0, zero_pt_, zero_pt_,
    1, 1, -1);
  segments.push_back(root_segment);

  std::vector<int> tubing_indexes = createTubingSegments(segments);
  std::vector<int> icd_indexes = createICDSegments(segments, tubing_indexes);
  createAnnulusSegments(segments, icd_indexes);
  std::vector<Segment> ann_segments;
  for (auto seg : segments) {
    if (seg.Type() == Segment::SegType::ANNULUS_SEGMENT) {
      ann_segments.push_back(seg);
    }
  }
  return ann_segments;
}

std::vector<int> Well::GetICDSegmentIndices() {
  std::vector<int> indices;
  for (auto seg : GetICDSegments()) {
    indices.push_back(seg.Index());
  }
  return indices;
}

}
}
