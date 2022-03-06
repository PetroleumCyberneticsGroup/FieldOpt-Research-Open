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
#include <Optimization/case_transfer_object.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/lexical_cast.hpp>
#include <sstream>

namespace {
using namespace Optimization;
using namespace boost::archive;

class CaseTransferObjectTest : public ::testing::Test, public TestResources::TestResourceCases {
 protected:
  CaseTransferObjectTest() {}
  virtual void SetUp() {}
  virtual void TearDown() {}

  // boost::uuids::uuid qUuidToBoostUuid(const QUuid quuid) const {
  //   boost::uuids::uuid uuid(string_generator()(quuid.toString().toStdString()));
  //   return uuid;
  // }
};

TEST_F(CaseTransferObjectTest, BaseConstructor) {
  CaseTransferObject cto();
}


TEST_F(CaseTransferObjectTest, CaseConstructor) {
  CaseTransferObject cto(test_case_3_4b3i3r_);
  EXPECT_STREQ(test_case_3_4b3i3r_->id_stdstr().c_str(), cto.id_stdstr().c_str());
  EXPECT_FLOAT_EQ(test_case_3_4b3i3r_->objf_value(),cto.objective_function_value());
  EXPECT_EQ(test_case_3_4b3i3r_->binary_variables().size(), cto.binary_variables().size());
  EXPECT_EQ(test_case_3_4b3i3r_->integer_variables().size(),cto.integer_variables().size());
  EXPECT_EQ(test_case_3_4b3i3r_->real_variables().size(),cto.real_variables().size());
  EXPECT_EQ(test_case_3_4b3i3r_->GetWICTime(), cto.wic_time_secs());
}

TEST_F(CaseTransferObjectTest, SerializationAndDeserializationAndGeneratedCase) {
  auto cto1 = CaseTransferObject(test_case_3_4b3i3r_); // create a populated cto

  // Serialize
  std::stringstream stream; // create a stream
  binary_oarchive oa(stream); // create an archive using the stream
  oa << cto1; // store the cto in the archive

  // Deserialize
  auto cto2 = CaseTransferObject(); // create an empty cto
  binary_iarchive ia(stream); // create an archive using the stream
  ia >> cto2; // deserialize the archive into cto2

  // Regenerate case object
  auto c = cto2.CreateCase();

  // Check that the generated case matches with the original test_case_3
  EXPECT_TRUE(test_case_3_4b3i3r_->Equals(c));
  EXPECT_STREQ(test_case_3_4b3i3r_->id_stdstr().c_str(), c->id_stdstr().c_str());
  EXPECT_FLOAT_EQ( test_case_3_4b3i3r_->objf_value(), c->objf_value());
  EXPECT_EQ(test_case_3_4b3i3r_->binary_variables().size(), c->binary_variables().size());
  EXPECT_EQ(test_case_3_4b3i3r_->integer_variables().size(), c->integer_variables().size());
  EXPECT_EQ(test_case_3_4b3i3r_->real_variables().size(), c->real_variables().size());
  EXPECT_EQ(test_case_3_4b3i3r_->GetWICTime(), c->GetWICTime());

  for (int ii=0; ii < test_case_3_4b3i3r_->integer_variables().size(); ii++) {
    EXPECT_EQ(test_case_3_4b3i3r_->integer_variables().at(ii).second,
              c->integer_variables().at(ii).second);
  }

  for (int ii=0; ii < test_case_3_4b3i3r_->binary_variables().size(); ii++) {
    EXPECT_EQ(test_case_3_4b3i3r_->binary_variables().at(ii).second,
              c->binary_variables().at(ii).second);
  }

  for (int ii=0; ii < test_case_3_4b3i3r_->real_variables().size(); ii++) {
    EXPECT_EQ(test_case_3_4b3i3r_->real_variables().at(ii).second,
              c->real_variables().at(ii).second);
  }

}
}
