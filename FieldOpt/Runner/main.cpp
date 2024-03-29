/***********************************************************
Copyright (C) 2015-2016
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

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;
#include <iostream>

#include "runners/main_runner.h"

int main(int argc, const char *argv[])
{
  try {
    // Parse runtime settings
    auto *runtime_settings = new Runner::RuntimeSettings(argc, argv);

    // Initialize runner
    auto runner = Runner::MainRunner(runtime_settings);
    runner.Execute();

  } catch (std::exception &e) {
    std::cout << "error: " << e.what() << std::endl;
    return 1;

  } catch (...) {
    std::cout << "Exception of unknown type." << std::endl;
  }
  return 0;
}

