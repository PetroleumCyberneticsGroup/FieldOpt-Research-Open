#ifndef FIELDOPT_TEST_RESOURCE_SETTINGS_H
#define FIELDOPT_TEST_RESOURCE_SETTINGS_H

#include "Utilities/settings/settings.h"
#include "Utilities/tests/test_resource_example_file_paths.h"
#include <iostream>

namespace TestResources {

    class TestResourceSettings {
    protected:
        TestResourceSettings() {

            settings_full_ = new Utilities::Settings::Settings(ExampleFilePaths::driver_example_, ExampleFilePaths::directory_output_);
            settings_optimizer_ = settings_full_->optimizer();
            settings_simulator_ = settings_full_->simulator();
            settings_model_ = settings_full_->model();

            // ESTABLISH OVERALL SETTINGS FOR TESTING
            if (settings_full_->verbosity()>2) std::cout << "... establishing overall settings for testing (test_resource_settings.h)" << std::endl;

            // ESTABLISH RESERVOIR BOUNDARY SETTINGS FOR TESTING
            constraint_settings_reservoir_boundary_.type = Utilities::Settings::Optimizer::ConstraintType::ReservoirBoundary;
            constraint_settings_reservoir_boundary_.box_imin = 0;
            constraint_settings_reservoir_boundary_.box_imax = 59;
            constraint_settings_reservoir_boundary_.box_jmin = 0;
            constraint_settings_reservoir_boundary_.box_jmax = 59;
            constraint_settings_reservoir_boundary_.box_kmin = 0;
            constraint_settings_reservoir_boundary_.box_kmax = 0;
            constraint_settings_reservoir_boundary_.well = "TESTW";
        }

        Utilities::Settings::Settings *settings_full_;
        Utilities::Settings::Optimizer *settings_optimizer_;
        Utilities::Settings::Simulator *settings_simulator_;
        Utilities::Settings::Model *settings_model_;

        Utilities::Settings::Optimizer::Constraint constraint_settings_reservoir_boundary_;
    };

}

#endif //FIELDOPT_TEST_RESOURCE_SETTINGS_H
