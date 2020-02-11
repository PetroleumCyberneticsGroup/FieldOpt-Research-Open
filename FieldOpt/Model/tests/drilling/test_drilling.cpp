#include <gtest/gtest.h>
#include "Model/tests/test_resource_model.h"

namespace {

class DrillingTest : public ::testing::Test,  public TestResources::TestResourceModel {
 protected:
  DrillingTest() { }

};

TEST_F(DrillingTest, Constructor) {
    EXPECT_TRUE(true);
}

}
