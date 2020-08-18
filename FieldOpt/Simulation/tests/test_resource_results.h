/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#ifndef FIELDOPT_TEST_RESOURCE_RESULTS_H
#define FIELDOPT_TEST_RESOURCE_RESULTS_H

#include "Settings/tests/test_resource_example_file_paths.hpp"
#include "Simulation/results/eclresults.h"
#include "Simulation/results/adgprsresults.h"

namespace TestResources {

class TestResourceResults {

 protected:
  TestResourceResults() {
    results_ecl_horzwell_ = new Simulation::Results::ECLResults();
    results_ecl_horzwell_->ReadResults(QString::fromStdString(ExampleFilePaths::ecl_base_horzwell));

    results_adgprs_5spot_ = new Simulation::Results::AdgprsResults();
    results_adgprs_5spot_->ReadResults(QString::fromStdString(ExampleFilePaths::gprs_smry_hdf5_5spot_));
  }

  virtual ~TestResourceResults() {
    results_ecl_horzwell_->DumpResults();
    results_adgprs_5spot_->DumpResults();
  }

 protected:
  Simulation::Results::ECLResults *results_ecl_horzwell_;
  Simulation::Results::AdgprsResults *results_adgprs_5spot_;
};
}

#endif //FIELDOPT_TEST_RESOURCE_RESULTS_H
