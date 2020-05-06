#include <gtest/gtest.h>
#include "Model/tests/test_resource_model.h"
#include "Simulation/simulator_interfaces/adgprssimulator.h"

namespace {

class AdgprsDriverFileWriterTest : public ::testing::Test, TestResources::TestResourceModel {
protected:
    AdgprsDriverFileWriterTest() {
        settings_full_->paths().SetPath(Paths::BUILD_DIR, TestResources::ExampleFilePaths::bin_dir_);
        settings_full_->paths().SetPath(Paths::SIM_DRIVER_FILE, TestResources::ExampleFilePaths::gprs_drv_5spot_);
        settings_full_->paths().SetPath(Paths::SIM_DRIVER_DIR, GetParentDirectoryPath(settings_full_->paths().GetPath(Paths::SIM_DRIVER_FILE)));
        simulator_ = new Simulation::AdgprsSimulator(settings_full_, model_);
    }
    virtual ~AdgprsDriverFileWriterTest() {}
    virtual void SetUp() {}
    Simulation::AdgprsSimulator *simulator_;
};

TEST_F(AdgprsDriverFileWriterTest, Initialization) {
}

}


