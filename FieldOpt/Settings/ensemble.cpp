/***********************************************************
Copyright (C) 2015-2018
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

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "ensemble.h"

#include "Utilities/filehandling.hpp"
#include "Utilities/verbosity.h"
#include "Utilities/printer.hpp"

using std::vector;
using std::string;
using std::pair;

namespace Settings {

using namespace Utilities::FileHandling;

Ensemble::Ensemble() {}

Ensemble::Ensemble(const std::string &ens_path) {

  ensemble_parent_dir_ = GetAbsoluteFilePath(GetParentDirectoryPath(ens_path));
  assert(FileExists(ens_path, false));
  assert(DirectoryExists(Ensemble::ensemble_parent_dir_, false));

  vector<string> file_contents = ReadFileToStdStringList(ens_path);
  for (auto line : file_contents) {

    if (line[0] == '#') {
      continue; // Skipping comment
    }

    vector<string> entries;
    boost::split(entries,line,boost::is_any_of(","));
    assert(entries.size() == 4);

    // Trim the whitespace
    boost::algorithm::trim(entries[0]);
    boost::algorithm::trim(entries[1]);
    boost::algorithm::trim(entries[2]);
    boost::algorithm::trim(entries[3]);

    string alias = entries[0];
    string data = ensemble_parent_dir_ + "/" + entries[1];
    string schedule = GetParentDirectoryPath(data) + "/" + entries[2];
    string grid = GetParentDirectoryPath(data) + "/" + entries[3];

    // Check that the alias has not already been used.
    assert(realizations_.count(entries[0]) == 0);
    assert(FileExists(data, false));
    assert(FileExists(schedule, false));
    assert(FileExists(grid, false));

    realizations_.insert(
        pair<string, Realization>
            ( alias,Realization(alias, data, schedule, grid) )
    );
  }
  n_select_ = realizations_.size();
}

Ensemble::Realization const &Ensemble::GetRealization(
    const string &alias) const {
  return realizations_.at(alias);
}

vector<string> Ensemble::GetAliases() const {
  vector<string> keys;
  for (auto realization : realizations_) {
    keys.push_back(realization.first);
  }
  return keys;
}

Ensemble::Realization::Realization(string alias,
                                   string data_rel_path,
                                   string schedule_rel_path,
                                   string grid_rel_path)
    : alias_(alias),
      data_rel_path_(data_rel_path),
      schedule_rel_path_(schedule_rel_path),
      grid_rel_path_(grid_rel_path) {
}

string Ensemble::Realization::alias() const { return alias_; }

string Ensemble::Realization::data() const { return data_rel_path_; }

string Ensemble::Realization::schedule() const { return schedule_rel_path_; }

string Ensemble::Realization::grid() const { return grid_rel_path_; }

int Ensemble::NSelect() const { return n_select_; }

void Ensemble::SetNSelect(const int n) { n_select_ = n; }

}
