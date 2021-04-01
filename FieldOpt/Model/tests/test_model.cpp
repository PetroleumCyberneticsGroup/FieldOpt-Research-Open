/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
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

#include <gtest/gtest.h>
#include "test_resource_model.h"

namespace {

class ModelTest : public ::testing::Test,
                  public TestResources::TestResourceModel {
 protected:
  ModelTest() { }

};

TEST_F(ModelTest, Constructor) {
  EXPECT_TRUE(true);
}

TEST_F(ModelTest, ChildObjects) {
  EXPECT_NO_THROW(model_->grid());
  EXPECT_NO_THROW(model_->wells());
  EXPECT_NO_THROW(model_->variables());
  EXPECT_GE(model_->wells()->size(), 2);
  EXPECT_GE(model_->variables()->BinaryVariableSize(), 0);
  EXPECT_GE(model_->variables()->DiscreteVariableSize(), 9);
  EXPECT_GE(model_->variables()->ContinuousVariableSize(), 5);
}

TEST_F(ModelTest, Variables) {
  // As of 2015.11.10, the variables are:
  // 3 Continuous variables (bhp at three time steps for the producer)
  QList<Model::Properties::ContinuousProperty *> prod_cont_variables = model_->variables()->GetWellControlVariables("PROD");
  EXPECT_EQ(prod_cont_variables.size(), 3);
  EXPECT_STREQ("BHP#PROD#0", prod_cont_variables[0]->name().toLatin1().constData());
  EXPECT_STREQ("BHP#PROD#50", prod_cont_variables[1]->name().toLatin1().constData());
  EXPECT_STREQ("BHP#PROD#365", prod_cont_variables[2]->name().toLatin1().constData());
  EXPECT_EQ(3, model_->wells()->at(0)->controls()->size());

  EXPECT_EQ(9, model_->variables()->GetWellBlockVariables("PROD").size());
  for (auto value : model_->variables()->GetDiscVarValues()) {
    EXPECT_GE(value.second, 0);
  }
}

TEST_F(ModelTest, ApplyCase) {
  Optimization::Case *c = new ::Optimization::Case(model_->variables()->GetBinVarValues(),
                                                   model_->variables()->GetDiscVarValues(),
                                                   model_->variables()->GetContVarValues());

  // Set all continous variables for the PROD well to 1. Should affect BHP and transmissibilty.
  auto producer_vars = model_->variables()->GetWellBHPVars("PROD");
  producer_vars.append(model_->variables()->GetTransmissibilityVariables("PROD"));
  for (auto var : producer_vars) {
    c->set_real_variable_value(var->id(), 1.0);
  }

  // Set all integer coordinates to 1 (should affect positions for all well blocks)
  auto producer_wb_vars = model_->variables()->GetWellBlockVariables("PROD");
  for (auto var : producer_wb_vars) {
    c->set_integer_variable_value(var->id(), 1);
  }

  model_->ApplyCase(c);

  for (Model::Wells::Control *control : *model_->wells()->first()->controls()) {
    EXPECT_FLOAT_EQ(1.0, control->bhp());
  }

  for (int i = 0; i < model_->wells()->first()->trajectory()->GetWellBlocks()->size(); ++i) {
    auto wb = model_->wells()->first()->trajectory()->GetWellBlocks()->at(i);
    if (i == 1) { //!< WB 1 is not variable, and should thus not have been changed
      EXPECT_EQ(2, wb->i());
      EXPECT_EQ(4, wb->j());
      EXPECT_EQ(1, wb->k());
    }
    else {
      EXPECT_EQ(1, wb->i());
      EXPECT_EQ(1, wb->j());
      EXPECT_EQ(1, wb->k());
    }
  }
}

TEST_F(ModelTest, Logging) {
  Optimization::Case *c = new ::Optimization::Case(model_->variables()->GetBinVarValues(),
                                                   model_->variables()->GetDiscVarValues(),
                                                   model_->variables()->GetContVarValues());

  // Set all continous variables for the PROD well to 1. Should affect BHP and transmissibilty.
  auto producer_vars = model_->variables()->GetWellBHPVars("PROD");
  producer_vars.append(model_->variables()->GetTransmissibilityVariables("PROD"));
  for (auto var : producer_vars) {
    c->set_real_variable_value(var->id(), 1.0);
  }

  // Set all integer coordinates to 1 (should affect positions for all well blocks)
  auto producer_wb_vars = model_->variables()->GetWellBlockVariables("PROD");
  for (auto var : producer_wb_vars) {
    c->set_integer_variable_value(var->id(), 1);
  }

  model_->ApplyCase(c);
  model_->SetResult("FieldOilProdTotal", vector<double>{1.0, 2.0, 3.0});
  model_->SetCompdatString("This is the compdat string.");
  model_->ApplyCase(c);
  model_->SetResult("FieldOilProdTotal", vector<double>{1.0, 2.0, 3.0});
  model_->SetResult("FieldWatProdTotal", vector<double>{4.0, 5.0, 6.0});
  model_->SetCompdatString("This is the compdat string.");
  model_->ApplyCase(c);
  model_->SetCompdatString("This is the compdat string.");
  model_->SetResult("FieldOilProdTotal", vector<double>{1.0, 2.0, 3.0});
  model_->SetResult("FieldWatProdTotal", vector<double>{4.0, 5.0, 6.0});
}

}
