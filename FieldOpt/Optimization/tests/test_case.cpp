/***********************************************************
Copyright (C) 2015-2016
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
#include "test_resource_cases.h"
#include <QList>

namespace {
class CaseTest : public ::testing::Test, public TestResources::TestResourceCases {
 public:
  CaseTest() : Test() { }
  virtual ~CaseTest() { }

 protected:
  virtual void TearDown() { }
  virtual void SetUp() { }


};

TEST_F(CaseTest, Constructor) {
  EXPECT_NO_THROW();
}

TEST_F(CaseTest, UUIDs) {
  EXPECT_FALSE(test_case_1_3i_->id().isNull());
  EXPECT_FALSE(test_case_2_3r_->id().isNull());
  EXPECT_FALSE(test_case_3_4b3i3r_->id().isNull());
  EXPECT_FALSE(test_case_4_4b3i3r->id().isNull());
}

TEST_F(CaseTest, ObjectiveValues) {
  EXPECT_THROW(test_case_1_3i_->objf_value(), Optimization::ObjectiveFunctionException);
  EXPECT_FLOAT_EQ(100.0, test_case_2_3r_->objf_value());
  EXPECT_FLOAT_EQ(-50.0, test_case_3_4b3i3r_->objf_value());
  EXPECT_FLOAT_EQ(-50.0, test_case_4_4b3i3r->objf_value());
}

TEST_F(CaseTest, Equals) {
  EXPECT_FALSE(test_case_1_3i_->Equals(test_case_2_3r_));
  EXPECT_TRUE(test_case_3_4b3i3r_->Equals(test_case_4_4b3i3r));
}

TEST_F(CaseTest, NumberOfVariables) {
  EXPECT_EQ(0, test_case_1_3i_->binary_variables().size());
  EXPECT_EQ(3, test_case_1_3i_->integer_variables().size());
  EXPECT_EQ(0, test_case_1_3i_->real_variables().size());

  EXPECT_EQ(0, test_case_2_3r_->binary_variables().size());
  EXPECT_EQ(0, test_case_2_3r_->integer_variables().size());
  EXPECT_EQ(3, test_case_2_3r_->real_variables().size());

  EXPECT_EQ(4, test_case_3_4b3i3r_->binary_variables().size());
  EXPECT_EQ(3, test_case_3_4b3i3r_->integer_variables().size());
  EXPECT_EQ(3, test_case_3_4b3i3r_->real_variables().size());

  EXPECT_EQ(4, test_case_4_4b3i3r->binary_variables().size());
  EXPECT_EQ(3, test_case_4_4b3i3r->integer_variables().size());
  EXPECT_EQ(3, test_case_4_4b3i3r->real_variables().size());
}

TEST_F(CaseTest, CopyConstructor) {
  Optimization::Case *copy = new Optimization::Case(test_case_1_3i_);
  EXPECT_FALSE(copy->id() == test_case_1_3i_->id());
  EXPECT_TRUE(copy->Equals(test_case_1_3i_));
}

TEST_F(CaseTest, PerturbMethod) {
  // Perturb the first case's first integer variable in positive direction
  QList<Optimization::Case *> perturbations_1 = test_case_1_3i_->Perturb(test_case_1_3i_->integer_variables().value(0).first,
                                                                         Optimization::Case::SIGN::PLUS, 7);
  EXPECT_EQ(1, perturbations_1.size());
  EXPECT_FALSE(perturbations_1.first()->id() == test_case_1_3i_->id());
  EXPECT_FALSE(perturbations_1.first()->Equals(test_case_1_3i_));

  EXPECT_TRUE(test_case_1_3i_->integer_variables().value(0).second + 7 == perturbations_1.first()->integer_variables().value(0).second);

  // Perturb the second case's first real variable in both directions
  QList<Optimization::Case *> perturbations_2 = test_case_2_3r_->Perturb(test_case_2_3r_->real_variables().value(0).first,
                                                                         Optimization::Case::SIGN::PLUSMINUS, 5);
  EXPECT_EQ(2, perturbations_2.size());
  EXPECT_FALSE(perturbations_2[0]->id() == test_case_2_3r_->id());
  EXPECT_FALSE(perturbations_2[1]->id() == test_case_2_3r_->id());
  EXPECT_FALSE(perturbations_2[0]->Equals(test_case_2_3r_));
  EXPECT_FALSE(perturbations_2[1]->Equals(test_case_2_3r_));
  EXPECT_FALSE(perturbations_2[0]->Equals(perturbations_2[1]));
  EXPECT_TRUE(test_case_2_3r_->real_variables().value(0).second + 5 == perturbations_2[0]->real_variables().value(0).second);
  EXPECT_TRUE(test_case_2_3r_->real_variables().value(0).second - 5 == perturbations_2[1]->real_variables().value(0).second);
}

TEST_F(CaseTest, VectorHelpers) {
  // Get the initial integer vector
  Eigen::VectorXi tc1_ivec_init = test_case_1_3i_->GetIntegerVarVector();

  // Create a delta vector
  Eigen::VectorXi delta_vec = Eigen::VectorXi(tc1_ivec_init.size());
  delta_vec[0] = 10;
  delta_vec[1] = 0;
  delta_vec[2] = 1;

  // Add the delta vector to the initial vector
  auto tc1_ivec_plus_delta = tc1_ivec_init + delta_vec;

  // Update the var values in the case to the new vector
  test_case_1_3i_->SetIntegerVarValues(tc1_ivec_plus_delta);

  // Get the updated vector
  Eigen::VectorXi tc1_updated = test_case_1_3i_->GetIntegerVarVector();

  // Check that the updated vector has the correct values
  EXPECT_EQ(tc1_updated[0], tc1_ivec_init[0] + delta_vec[0]);
  EXPECT_EQ(tc1_updated[1], tc1_ivec_init[1] + delta_vec[1]);
  EXPECT_EQ(tc1_updated[2], tc1_ivec_init[2] + delta_vec[2]);
}

}
