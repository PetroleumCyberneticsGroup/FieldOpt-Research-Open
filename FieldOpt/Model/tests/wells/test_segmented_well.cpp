/***********************************************************
Created by einar on 6/6/18.
Copyright (C) 2018
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

#include <gtest/gtest.h>
#include <algorithm>
#include <wells/well.h>
#include "Settings/settings.h"
#include "Settings/tests/test_resource_multispline_wells_settings.hpp"
#include "Settings/paths.h"
#include "Settings/tests/test_resource_example_file_paths.hpp"

using namespace Model::Wells;

namespace {

class SegmentedWellTest : public ::testing::Test {
 protected:
  SegmentedWellTest() {
    auto partial_deck = model_adtl_pts;
    paths_.SetPath(Paths::SIM_DRIVER_FILE, TestResources::ExampleFilePaths::deck_flow_5spot_);
    paths_.SetPath(Paths::GRID_FILE, TestResources::ExampleFilePaths::grid_5spot_);
    // paths_.SetPath(Paths::SIM_SCH_FILE, TestResources::ExampleFilePaths::norne_sch_);
    mod_json_ = partial_deck["Model"].toObject();
    mod_settings_ = new Settings::Model(mod_json_, paths_);
    varcont_ = new Model::Properties::VariablePropertyContainer();
    grid_ = new Reservoir::Grid::ECLGrid(paths_.GetPath(Paths::GRID_FILE));
  }

  QJsonObject sim_json_;
  QJsonObject mod_json_;
  Settings::Model *mod_settings_;
  Model::Properties::VariablePropertyContainer *varcont_;
  Reservoir::Grid::Grid *grid_;
  Paths paths_;
  ::Model::Wells::Well *well;
};

TEST_F(SegmentedWellTest, Constructor ) {
  // Verifying settings
  EXPECT_TRUE(mod_settings_->wells()[0].use_segmented_model);
  EXPECT_EQ(3, mod_settings_->wells()[0].seg_n_compartments);

  // Construct well
  well = new Model::Wells::Well(*mod_settings_, 0, varcont_, grid_, nullptr);
  double length = well->trajectory()->GetLength();
}

TEST_F(SegmentedWellTest, Compartments ) {
  well = new Model::Wells::Well(*mod_settings_, 0, varcont_, grid_, nullptr);
  double length = well->trajectory()->GetLength();

  EXPECT_EQ(3, well->GetCompartments().size());
  EXPECT_EQ(4, well->GetPackers().size());
  EXPECT_EQ(3, well->GetICDs().size());

  EXPECT_FLOAT_EQ(0.0, well->GetCompartments()[0].start_packer->md(length));

  // some error is allowed here, because the packers are snapped to
  // spline-cell intersections.
  EXPECT_NEAR(length, well->GetCompartments()[2].end_packer->md(length), 6);

  EXPECT_TRUE(well->GetCompartments()[0].end_packer == well->GetCompartments()[1].start_packer);
  EXPECT_TRUE(well->GetCompartments()[1].end_packer == well->GetCompartments()[2].start_packer);
  EXPECT_TRUE(well->GetCompartments()[0].end_packer->md(length) > well->GetCompartments()[0].start_packer->md(length));
  EXPECT_TRUE(well->GetCompartments()[1].end_packer->md(length) > well->GetCompartments()[1].start_packer->md(length));
  EXPECT_TRUE(well->GetCompartments()[2].end_packer->md(length) > well->GetCompartments()[2].start_packer->md(length));

  EXPECT_NEAR(0.0, well->GetCompartments()[0].start_packer->md(length), 1.0);
  EXPECT_NEAR(length, well->GetCompartments()[2].end_packer->md(length), 1.0);

  EXPECT_FLOAT_EQ(7.85E-5, well->GetCompartments()[0].icd->valveSize());
  EXPECT_FLOAT_EQ(7.85E-5, well->GetCompartments()[1].icd->valveSize());
  EXPECT_FLOAT_EQ(7.85E-5, well->GetCompartments()[2].icd->valveSize());
  EXPECT_NEAR(0.0,          well->GetCompartments()[0].icd->md(length), 30);
  EXPECT_NEAR(well->GetCompartments()[2].start_packer->md(length), well->GetCompartments()[2].icd->md(length), 1);
}

TEST_F(SegmentedWellTest, SegmentTypes) {
  well = new Model::Wells::Well(*mod_settings_, 0, varcont_, grid_, nullptr);
  auto segs = well->GetSegments();
  EXPECT_EQ(Segment::SegType::TUBING_SEGMENT, segs[0].Type()); // Root segment
  for (int i = 1; i < 4; ++i) {
    EXPECT_EQ(Segment::SegType::TUBING_SEGMENT, segs[i].Type());
  }
  for (int i = 4; i < 7; ++i) {
    EXPECT_EQ(Segment::SegType::ICD_SEGMENT, segs[i].Type());
  }
  for (int i = 7; i < segs.size(); ++i) {
    EXPECT_EQ(Segment::SegType::ANNULUS_SEGMENT, segs[i].Type());
  }
}

TEST_F(SegmentedWellTest, SegmentConnections) {
  well = new Model::Wells::Well(*mod_settings_, 0, varcont_, grid_, nullptr);
  auto segs = well->GetSegments();
  auto tub_segs = well->GetTubingSegments();
  auto icd_segs = well->GetICDSegments();
  auto ann_segs = well->GetAnnulusSegments();
  EXPECT_EQ(4, tub_segs.size()); // Includes root segment
  EXPECT_EQ(3, icd_segs.size());
  EXPECT_EQ(segs.size() - 7, ann_segs.size());

  for (int i = 1; i < tub_segs.size(); ++i) {
    EXPECT_EQ(tub_segs[i].Outlet(), tub_segs[i-1].Index());

    // Check that the next tub index is in the list of inlets for the prev one.
    EXPECT_TRUE(std::find(tub_segs[i-1].GetInlets().begin(),
                          tub_segs[i-1].GetInlets().end(),
                          tub_segs[i].Index()) != tub_segs[i-1].GetInlets().end());
  }
  for (int i = 0; i < icd_segs.size(); ++i) {
    EXPECT_EQ(icd_segs[i].Outlet(), tub_segs[i].Index());
  }
  for (auto seg : segs) {
    std::cout << seg.ToString();
  }
}

}
