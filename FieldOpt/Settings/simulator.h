/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#ifndef SETTINGS_SIMULATOR_H
#define SETTINGS_SIMULATOR_H

#include "settings.h"
#include "Settings/paths.h"
#include "Settings/ensemble.h"

#include <QStringList>

namespace Settings {

/*!
 * \brief The Simulator class contains simulator-specific
 * settings. Simulator settings objects may _only_ be
 * created by the Settings class. They are created when
 * reading a JSON-formatted "driver file".
 */
class Simulator
{
  friend class Settings;

 public:
  Simulator(const QJsonObject& json_simulator, Paths &paths, VerbParams vp);

  enum SimulatorType { ECLIPSE, ADGPRS, Flow, INTERSECT };
  enum SimulatorFluidModel { BlackOil, DeadOil };
  enum FileStructType { Flat, Branched };

  /*! @brief Get simulator type (e.g. ECLIPSE). */
  SimulatorType type() const { return type_; }

  /*! @brief Get simulator commands (commands used to execute a
   * simulation). Each list element is executed in sequence. */
  QStringList *sim_exec_cmds() const { return sim_exec_cmds_; }

  /*! @brief Get name of script used to execute simulations. */
  QString script_name() const { return script_name_; }

  bool is_ensemble() const { return is_ensemble_; }

  Ensemble get_ensemble() const { return ensemble_; }

  /*! @brief Get the fluid model. */
  SimulatorFluidModel fluid_model() const { return fluid_model_; }

  /*! @brief Get max # of minutes simulations are allowed to
   * run if no timeout value can be calculated. Returns -1
   * if field is not set. */
  int max_minutes() { return max_minutes_; }

  /*! @brief Check whether ACTIONX keywords should be used
   * when writing the driver file. Note that this disables
   * some functionality. */
  bool use_actionx() const { return ecl_use_actionx_; }

  /*! @brief Check whether or not to call a script named
   * FO_POSTSIM.sh in the same directory as the simulator
   * DATA file after simulation, if the script is found. */
  bool use_post_sim_script() const { return use_post_sim_script_; }
  bool use_pre_sim_script() const { return use_pre_sim_script_; }
  bool add_sim_scripts() const { return add_sim_scripts_; }

  QStringList* pre_sim_args() const { return pre_sim_args_; };
  QStringList* post_sim_args() const { return post_sim_args_; };

  /*! @brief Check whether or not to read external results
   * component from a file named FO_EXT_RESULTS.json in the
   * same directory as the simulator DATA file after simulation
   * (and potential call of PostSimScript). */
  bool read_external_json_results() const { return read_external_json_results_; }

  bool read_adj_grad_data() const { return read_adj_grad_data_; }

  struct FileStructure {
    FileStructType type_ = FileStructType::Flat;
    int levels_num_ = 0;
    string levels_str_ = "";
  };

  /*! @brief Get the fluid model. */
  FileStructure file_structure() { return file_structure_; }

  VerbParams verbParams() { return vp_; };

 private:
  SimulatorType type_;
  SimulatorFluidModel fluid_model_;

  string md_ = "Settings";
  string cl_ = "Simulator";
  VerbParams vp_;

  QString script_name_;
  QStringList* sim_exec_cmds_;
  QStringList* pre_sim_args_;
  QStringList* post_sim_args_;

  bool is_ensemble_ = false;
  bool ecl_use_actionx_ = false;

  bool use_post_sim_script_ = false;
  bool use_pre_sim_script_ = false;
  bool add_sim_scripts_ = false;
  QString sched_file_name_ = "";

  bool read_external_json_results_ = false;
  bool read_adj_grad_data_ = false;
  int max_minutes_ = -1;
  Ensemble ensemble_;

  FileStructure file_structure_;

  void setStructure(QJsonObject json_simulator);
  void setPaths(QJsonObject json_simulator, Paths &paths);
  void setType(QJsonObject json_simulator);
  void setParams(QJsonObject json_simulator);
  void setCommands(QJsonObject json_simulator);
  void setFluidModel(QJsonObject json_simulator);

};

}


#endif // SETTINGS_SIMULATOR_H
