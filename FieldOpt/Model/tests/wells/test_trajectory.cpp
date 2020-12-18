/***********************************************************
Copyright (C) 2017
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
#include "Reservoir/tests/test_resource_grids.h"
#include "Model/tests/test_resource_model.h"
#include "Model/wells/well_exceptions.h"
#include "Settings/tests/test_resource_example_file_paths.hpp"

using namespace Model::Wells;

namespace {

class TrajectoryTest : public ::testing::Test,
                       public TestResources::TestResourceModel,
                       public TestResources::TestResourceGrids
{

 protected:
  TrajectoryTest() {
    prod_well_trajectory_ = model_->wells()->first()->trajectory();
  }
  virtual ~TrajectoryTest() {}
  virtual void SetUp() {}
  ::Model::Wells::Wellbore::Trajectory *prod_well_trajectory_;
  Settings::VerbParams vp_ = {};
};

TEST_F(TrajectoryTest, Constructor) {
  EXPECT_TRUE(true);
}

TEST_F(TrajectoryTest, GetWellBlock) {
  EXPECT_NO_THROW(prod_well_trajectory_->GetWellBlock(1, 4, 1));
  EXPECT_THROW(prod_well_trajectory_->GetWellBlock(9, 9, 9), WellBlockNotFoundException);

  EXPECT_EQ(1, prod_well_trajectory_->GetWellBlock(1,4,1)->i());
  EXPECT_EQ(4, prod_well_trajectory_->GetWellBlock(1,4,1)->j());
  EXPECT_EQ(1, prod_well_trajectory_->GetWellBlock(1,4,1)->k());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlock(1,4,1)->HasCompletion());

  EXPECT_EQ(2, prod_well_trajectory_->GetWellBlock(2,4,1)->i());
  EXPECT_EQ(4, prod_well_trajectory_->GetWellBlock(2,4,1)->j());
  EXPECT_EQ(1, prod_well_trajectory_->GetWellBlock(2,4,1)->k());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlock(2,4,1)->HasCompletion());

  EXPECT_EQ(3, prod_well_trajectory_->GetWellBlock(3,4,1)->i());
  EXPECT_EQ(4, prod_well_trajectory_->GetWellBlock(3,4,1)->j());
  EXPECT_EQ(1, prod_well_trajectory_->GetWellBlock(3,4,1)->k());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlock(3,4,1)->HasCompletion());

  EXPECT_EQ(4, prod_well_trajectory_->GetWellBlock(4,4,1)->i());
  EXPECT_EQ(4, prod_well_trajectory_->GetWellBlock(4,4,1)->j());
  EXPECT_EQ(1, prod_well_trajectory_->GetWellBlock(4,4,1)->k());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlock(4,4,1)->HasCompletion());
}

TEST_F(TrajectoryTest, AllWellBlocks) {
  EXPECT_EQ(4, prod_well_trajectory_->GetWellBlocks()->size());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlocks()->at(0)->HasCompletion());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlocks()->at(1)->HasCompletion());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlocks()->at(2)->HasCompletion());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlocks()->at(3)->HasCompletion());
}

TEST_F(TrajectoryTest, Completions) {
  EXPECT_EQ(::Model::Wells::Wellbore::Completions::Completion::CompletionType::Perforation,
            prod_well_trajectory_->GetWellBlock(1,4,1)->GetCompletion()->type());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlock(1,4,1)->HasPerforation());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlock(2,4,1)->HasPerforation());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlock(3,4,1)->HasPerforation());
  EXPECT_TRUE(prod_well_trajectory_->GetWellBlock(4,4,1)->HasPerforation());
  EXPECT_FLOAT_EQ(1.0, prod_well_trajectory_->GetWellBlock(1,4,1)->GetPerforation()->transmissibility_factor());
  EXPECT_FLOAT_EQ(1.0, prod_well_trajectory_->GetWellBlock(3,4,1)->GetPerforation()->transmissibility_factor());
}

TEST_F(TrajectoryTest, VariableContainerConsistencyAfterCreation) {
  // There should be three integer variables (i,j,k) for each of the four well block
  EXPECT_EQ(9, model_->variables()->DiscreteVariableSize());
}


TEST_F(TrajectoryTest, MultisplineWell) {
  Paths paths;
  paths.SetPath(Paths::GRID_FILE, TestResources::ExampleFilePaths::grid_5spot_);
  auto settings = Settings::Model(TestResources::TestResourceModelSettingSnippets::model_adtl_pts(), paths, vp_);
  auto wsettings = settings.wells()[0];
  auto varcont = new Model::Properties::VarPropContainer(vp_);
  auto well = Model::Wells::Wellbore::WellSpline(wsettings, varcont, TestResources::TestResourceGrids::grid_5spot_, nullptr);
  auto well_blocks = well.GetWellBlocks();

  EXPECT_EQ(well_blocks->size(), 26);
  EXPECT_EQ(well_blocks->first()->i(), 13);
  EXPECT_EQ(well_blocks->first()->j(), 38);
  EXPECT_EQ(well_blocks->last()->i(), 38);
  EXPECT_EQ(well_blocks->last()->j(), 38);
}


}

