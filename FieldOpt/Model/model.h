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

#ifndef MODEL_H
#define MODEL_H

#include <QString>
#include <QList>
#include "Reservoir/grid/eclgrid.h"
#include "properties/var_prop_container.h"
#include "wells/well.h"
#include "WellIndexCalculation/wicalc_rixx.h"
#include "Settings/model.h"
#include "Optimization/case.h"
#include "Model/wells/wellbore/wellblock.h"

#include "Runner/loggable.hpp"
#include "Runner/logger.h"

class Logger;

namespace Model {
class ModelSynchronizationObject;

using Reservoir::WellIndexCalculation::wicalc_rixx;
using Reservoir::Grid::ECLGrid;
using Optimization::Constraints::ConstraintHandler;

using WDefType = Settings::Model::WellDefinitionType;
using ES = Optimization::Case::CaseState::EvalStatus;

using Printer::info;
using Printer::idbg;
using Printer::ext_info;
using Printer::ext_warn;
using Printer::ext_error;
using Printer::num2str;

/*!
 * \brief The Model class represents the reservoir model
 * as a whole, including wells and any related variables,
 * and the reservoir grid.
 */
class Model : public Loggable
{
  friend class ModelSynchronizationObject;
 public:
  Model(::Settings::Settings settings, Logger *logger);

  LogTarget GetLogTarget() override;
  map<string, string> GetState() override;
  QUuid GetId() override;
  map<string, vector<double>> GetValues() override;

  /*!
   * \brief reservoir Get the reservoir (i.e. grid).
   */
  Reservoir::Grid::Grid *grid() const { return grid_; }

  void set_grid_path(const std::string &grid_path);

  /*!
   * \brief variables Get the set of variable properties of all types.
   */
  Properties::VarPropContainer *variables() const {
    return var_container_;
  }

  /*!
   * \brief wells Get a list of all the wells in the model.
   */
  QList<Wells::Well *> *wells() const { return wells_; }

  /*!
   * \brief ApplyCase Applies the variable values
   * from a case to the variables in the model.
   * \param c Case to apply the variable values of.
   */
  void ApplyCase(Optimization::Case *c);

  /*!
   * @brief Get the UUId of last case applied to the Model.
   * @return
   */
  QUuid GetCurrentCaseId() const { return current_case_id_; }

  void SetCompdatString(const QString compdat) { compdat_ = compdat; };

  void SetResult(const std::string key, std::vector<double> vec);

  /*!
   * @brief Should be called at the end of the optimization
   * run. Writes the last case to the extended log.
   */
  void Finalize();

  ConstraintHandler * constraintHandler() {
    return constraint_handler_;
  }

 private:
  Reservoir::Grid::Grid *grid_;
  Reservoir::WellIndexCalculation::wicalc_rixx *wic_;
  Properties::VarPropContainer *var_container_;

  QList<Wells::Well *> *wells_;

  //!< Verify the model. Throws an exception if it is not.
  void verify();

  void verifyWells();
  void verifyWellTrajectory(Wells::Well *w);
  void verifyWellBlock(Wells::Wellbore::WellBlock *wb);
  void verifyWellCompartments(Wells::Well *w);

  Logger *logger_;
  QUuid current_case_id_;

  //!< Pointer to current case. Kept for logging purposes.
  Optimization::Case *current_case_;

  //!< The compdat generated from the list of well blocks
  //!< corresponding to the current case. This is set by
  //!< the simulator library.
  QString compdat_;

  //!< The results of the last simulation (i.e.
  //!< the one performed with the current case).
  std::map<std::string, std::vector<double>> results_;

  QHash<QString, double> realization_ofv_map_;
  double ensemble_ofv_st_dev_;
  double ensemble_avg_ofv_;

  class Summary : public Loggable {
   public:
    Summary(Model *model) { model_  = model; }
    LogTarget GetLogTarget() override;
    map<string, string> GetState() override;
    map<string, WellDescription> GetWellDescriptions() override;
    QUuid GetId() override;
    map<string, vector<double>> GetValues() override;
   private:
    Model *model_;
  };

  /*! The Economy struct is intended for use in calculations
   * for drilling costs. The costs are represented in [$/m]
   * if used in tandem with NPV
   */
 public:
  struct Economy{
    map<string, double> well_xy;
    map<string, double> well_z;
    map<string, double> well_lengths;
    double costXY = 0;
    double costZ = 0;
    double cost = 0;
    bool separate = false;
    bool use_well_cost = false;
    QList<Wells::Well *> wells_pointer;
  };

  Economy well_economy_;
  void wellCost(Settings::Optimizer *settings);
  Economy* wellCostConstructor();

 private:
  vector<double> objf_scal_;
  ConstraintHandler *constraint_handler_;

  Settings::VerbParams vp_;

  string md_ = "Model";
  string cl_ = "Model";
  string im_ = "", wm_ = "", em_ = "";

 public:
  void setObjfScaler(vector<double> objf_scal) { objf_scal_ = objf_scal; }
  int getObjfScalSz() { return objf_scal_.size(); }

  double getObjfScalVal() {
    auto eps = std::numeric_limits<double>::epsilon();
    double scal = accumulate(objf_scal_.begin(), objf_scal_.end(), 0.0);

    if (scal >= -eps && scal <= eps) {
      return scal + 1.0;
    } else {
      return scal;
    }
  }

};

}

#endif // MODEL_H
