#include <gtest/gtest.h>
#include "Model/properties/property.h"
#include "Model/properties/binary_property.h"
#include "Model/properties/discrete_property.h"
#include "Model/properties/continous_property.h"
#include "Model/properties/var_prop_container.h"
#include "Model/tests/test_resource_model.h"

using namespace Model::Properties;

namespace {

    class VariableContainerTest : public ::testing::Test, TestResources::TestResourceModel {
    protected:
        VariableContainerTest () {
            vc_ = model_->variables();
        }
        virtual ~VariableContainerTest () { }

        virtual void SetUp() { }

        VarPropContainer *vc_;
    };

    TEST_F(VariableContainerTest, Constructor) {
        // Check that all variables have been added
        // PROD well: 9 int pos vars, 3 float trans vars, 3 float bhp vars
        // INJ well: 6 float pos vars
        // TESTW well: 6 float pos vars, 1 float bhp var
        // Total: 12 int vars, 20 float vars, no binary vars
        EXPECT_EQ(9, vc_->GetDiscreteVariables()->size());
        EXPECT_EQ(19, vc_->GetContinuousVariables()->size());
        EXPECT_EQ(0, vc_->GetBinaryVariables()->size());
    }

    TEST_F(VariableContainerTest, PRODWellVariables) {
        EXPECT_EQ(9, vc_->GetWellBlockVariables("PROD").length());
        EXPECT_EQ(3, vc_->GetTransmissibilityVariables("PROD").length());
        EXPECT_EQ(3, vc_->GetWellBHPVars("PROD").length());
        EXPECT_EQ(0, vc_->GetWellRateVars("PROD").length());
        EXPECT_EQ(0, vc_->GetWSplineVars("PROD").length());
    }


    TEST_F(VariableContainerTest, INJWellVariables) {
        EXPECT_EQ(6, vc_->GetWSplineVars("INJ").length());
        EXPECT_EQ(0, vc_->GetWellRateVars("INJ").length());
        EXPECT_EQ(0, vc_->GetWellBlockVariables("INJ").length());
        EXPECT_EQ(0, vc_->GetTransmissibilityVariables("INJ").length());
        EXPECT_EQ(0, vc_->GetWellBHPVars("INJ").length());
    }

}
