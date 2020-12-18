/***********************************************************
Copyright (C) 2015-2017
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

#ifndef CASE_H
#define CASE_H

#include <QMap>
#include <QHash>
#include <QUuid>
#include <Utilities/math.hpp>
#include <Eigen/Core>
#include <QtCore/QDateTime>
#include <Model/properties/var_prop_container.h>
#include <Utilities/verbosity.h>
#include "Runner/loggable.hpp"
#include "optimization_exceptions.h"

namespace Optimization {

using Printer::ext_warn;
using Printer::ext_info;
using Printer::info;

class CaseHandler;
class CaseTransferObject;

/*!
 * \brief The Case class represents a specific case for the
 * optimizer, i.e. a specific set of variable values and
 * the value of the objective function after evaluation.
 */
class Case : public Loggable
{
 public:
  friend class CaseHandler;
  friend class CaseTransferObject;

  Case();
  // Case(const QHash<QUuid, bool> &binary_variables,
  //      const QHash<QUuid, int> &integer_variables,
  //      const QHash<QUuid, double> &real_variables);

  // Case(const QMap<QUuid, bool> &binary_variables,
  //      const QMap<QUuid, int> &integer_variables,
  //      const QMap<QUuid, double> &real_variables);

  Case(const QList<QPair<QUuid, bool>> &binary_variables,
       const QList<QPair<QUuid, int>> &integer_variables,
       const QList<QPair<QUuid, double>> &real_variables);

  Case(const Case &c) = delete;
  Case(const Case *c);

  /*!
   * @brief The CaseState struct holds information about the
   * current status of the Case object, such as whether or
   * not it has been evaluated and whether or not it has been
   * modified by a constriant.
   */
  struct CaseState {
    enum EvalStatus : int {
      E_FAILED=-2, E_TIMEOUT=-1,
      E_PENDING=0,
      E_CURRENT=1, E_DONE=2,
      E_BOOKKEEPED=3
    };

    enum ConsStatus : int {
      C_PROJ_FAILED=-2, C_INFEASIBLE=-1, C_PENDING=0,
      C_FEASIBLE=1, C_PROJECTED=2, C_PENALIZED=3,
    };

    enum QueueStatus : int {
      Q_DISCARDED=-1,
      Q_QUEUED=0,
      Q_DEQUEUED=1
    };

    enum ErrorMessage : int {
      ERR_SIM=-4, ERR_WIC=-3, ERR_CONS=-2, ERR_UNKNOWN=-1,
      ERR_OK=0
    };

    CaseState() {
      eval = E_PENDING;
      cons = C_PENDING;
      queue = Q_QUEUED;
      err_msg = ERR_OK;
    }

    EvalStatus eval;
    ConsStatus cons;
    QueueStatus queue;
    ErrorMessage err_msg;
  };

  //!< The state of the Case, directly modifiable.
  CaseState state;

  /*!
   * \brief Equals Checks whether this case is equal
   * to another case within some tolerance.
   * \param other Case to compare with.
   * \param tolerance The allowed deviation between two cases.
   * \return True if the cases are equal within the tolerance,
   * otherwise false.
   */
  bool Equals(const Case *other, double tolerance=0.0) const;

  QUuid id() const { return id_; }

  //!< Get an std string representation of the case uuid.
  string id_stdstr() { return id_.toString().toStdString(); }

  /*!
   * @brief Get a string representation of this case,
   * suitable for console printing.
   * @param varcont Pointer to the variable container.
   * This is needed to get variable names.
   * @return An std string describing the case.
   */
  string StringRepresentation(Model::Properties::VarPropContainer *varcont);

  // QHash<QUuid, bool> binary_variables() const { return binary_variables_; }
  // QHash<QUuid, int> integer_variables() const { return integer_variables_; }
  // QHash<QUuid, double> real_variables() const { return real_variables_; }

  // QMap<QUuid, bool> binary_variables() const { return binary_variables_; }
  // QMap<QUuid, int> integer_variables() const { return integer_variables_; }
  // QMap<QUuid, double> real_variables() const { return real_variables_; }

  QList<QPair<QUuid, bool>> binary_variables() const {
    return binary_variables_;
  }
  QList<QPair<QUuid, int>> integer_variables() const {
    return integer_variables_;
  }
  QList<QPair<QUuid, double>> real_variables() const {
    return real_variables_;
  }

  // void set_binary_variables(const QHash<QUuid, bool> &binary_variables) { binary_variables_ = binary_variables; }
  // void set_integer_variables(const QHash<QUuid, int> &integer_variables) { integer_variables_ = integer_variables; }
  // void set_real_variables(const QHash<QUuid, double> &real_variables) { real_variables_ = real_variables; }

  // void set_binary_variables(const QMap<QUuid, bool> &binary_variables) { binary_variables_ = binary_variables; }
  // void set_integer_variables(const QMap<QUuid, int> &integer_variables) { integer_variables_ = integer_variables; }
  // void set_real_variables(const QMap<QUuid, double> &real_variables) { real_variables_ = real_variables; }

  void set_binary_variables(const QList<QPair<QUuid, bool>> &binary_variables) {
    binary_variables_ = binary_variables;
  }

  void set_integer_variables(const QList<QPair<QUuid, int>> &integer_variables) {
    integer_variables_ = integer_variables;
  }

  void set_real_variables(const QList<QPair<QUuid, double>> &real_variables) {
    real_variables_ = real_variables;
  }

  double objf_value() const; //!< Get the objective function value. Throws an exception if the value has not been defined.
  void set_objf_value(double objective_function_value);

  //!< Set the value of an integer variable in the case.
  void set_integer_variable_value(QUuid id, int val);

  //!< Set the value of a boolean variable in the case.
  void set_binary_variable_value(QUuid id, bool val);

  //!< Set the value of a real variable in the case.
  void set_real_variable_value(QUuid id, double val);

  int get_integer_variable_value(QUuid id) {
    for (const auto & integer_variable : integer_variables_) {
      if (integer_variable.first==id) {
        return integer_variable.second;
      }
    }
  };

  bool get_binary_variable_value(QUuid id) {
    for (const auto & binary_variable : binary_variables_) {
      if (binary_variable.first==id) {
        return binary_variable.second;
      }
    }
  };

  double get_real_variable_value(QUuid id) {
    for (const auto & real_variable : real_variables_) {
      if (real_variable.first==id) {
        return real_variable.second;
      }
    }
  };

  enum SIGN { PLUS, MINUS, PLUSMINUS};

  /*!
   * \brief Perturb Creates variations of this Case where the value of one variable has been changed.
   *
   * If PLUS or MINUS is selected as the sign, _one_ case is returned. If PLUSMINUS is selected, _two_
   * cases are returned.
   *
   * Note that this method only handles integer and real variables.
   * \param variabe_id The UUID of the variable to be perturbed.
   * \param sign The sign/direction of the perturbation.
   * \param magnitude The magnitude of the perturbaton.
   * \return One or two cases where one variable has been perturbed.
   */
  QList<Case *> Perturb(QUuid variabe_id, SIGN sign, double magnitude);

  /*!
   * Get the real variables of this case as a Vector.
   *
   * @note This function will not work with Case objects created from CaseTransferObject.
   * This implies that, when running in parallel, it will only work on the main process.
   *
   * This creates an ordering of the variables so that for future
   * use the i'th index in the vector will always correspond to the
   * same variable.
   * @return Values of the real variables in a vector
   */
  Eigen::VectorXd GetRealVarVector();

  /*!
   * Sets the real variable values of this case from a given vector.
   *
   * @note This function will not work with Case objects created from CaseTransferObject.
   * This implies that, when running in parallel, it will only work on the main process.
   *
   * The order of the variables as they appear in vector this case is preserved
   * given that they were taken from this same case from the function GetRealVector()
   * @param vec
   */
  void SetRealVarValues(Eigen::VectorXd vec);

  /*!
   * @brief Get a vector containing the variable UUIDs in the same order they appear
   * in in the vector from GetRealVarVector.
   */
  QList<QUuid> GetRealVarIdVector() { return real_id_index_map_; }

  /*!
   * Get the integer variables of this case as a Vector.
   *
   * @note This function will not work with Case objects created from CaseTransferObject.
   * This implies that, when running in parallel, it will only work on the main process.
   *
   * This creates an ordering of the variables so that for future
   * use the i'th index in the vector will always correspond to the
   * same variable.
   * @return Values of the integer variables in a vector
   */
  Eigen::VectorXi GetIntegerVarVector();

  /*!
   * Sets the integer variable values of this case from a given vector.
   *
   * @note This function will not work with Case objects created from CaseTransferObject.
   * This implies that, when running in parallel, it will only work on the main process.
   *
   * The order of the variables as they appear in vector this case is preserved
   * given that they were taken from this same case from the function GetIntegerVarVector()
   * @param vec
   */
  void SetIntegerVarValues(Eigen::VectorXi vec);

  /*!
   * @brief Set the origin info of this Case/trial point, i.e. which point it was generated
   * from, in which direction it was perturbed, and with what magnitude. This method is
   * needed by some optimization algorithms.
   * @param parent The Case/trial point this was generated from.
   * @param direction_index The direction index of the perturbation.
   * @param step_length The magnitude of the perturbation.
   */
  void set_origin_data(Case* parent, int direction_index, double step_length);

  Case* origin_case() const { return parent_; }
  int origin_direction_index() const { return direction_index_; }
  double origin_step_length() const { return step_length_; }

  void SetSimTime(const int sec) { sim_time_sec_ = sec; }
  int GetSimTime() const { return sim_time_sec_; }

  // Logger interface
  LogTarget GetLogTarget() override;
  map<string, string> GetState() override;
  QUuid GetId() override;
  map<string, vector<double>> GetValues() override;
  // End Logger interface

  /*!
   * @brief Set the time spent computing the well blocs for this case.
   * @param secs The number of seconds it took.
   */
  void SetWICTime(const int secs) { wic_time_sec_ = secs; }

  /*!
   * @brief Get the number of seconds spent computing the well blocks for this case.
   */
  int GetWICTime() const { return wic_time_sec_; }

  // Multiple realizations-support
  void SetEnsembleRealization(const QString &alias) { ensemble_realization_ = alias; }
  QString GetEnsembleRealization() const { return ensemble_realization_; }

  void SetRealizationOfv(const QString &alias, const double &ofv);
  bool HasRealizationOfv(const QString &alias);
  double GetRealizationOfv(const QString &alias);
  double GetEnsembleAverageOfv() const;
  /*!
   * Gets Ensemble Expected Objective Function Value (OFV). This includes the average
   * and standard deviation for all the OFVs in the ensemble that was previously run.
   */

  QPair<double, double> GetEnsembleExpectedOfv() const;
  QHash<QString, double> GetRealizationOFVMap() const { return ensemble_ofvs_; }

  bool integer_vars_contain(int &var_indx, QUuid var_id) {
    for(int ii=0; ii < integer_variables_.size(); ii++) {
      if (integer_variables_.at(ii).first==var_id) {
        var_indx = ii;
        return true;
      }
    }
  }

  bool real_vars_contain(int &var_indx, QUuid var_id) {
    for(int ii=0; ii < real_variables_.size(); ii++) {
      if (real_variables_.at(ii).first==var_id) {
        var_indx = ii;
        return true;
      }
    }
  }

  void SetVerbParams(Settings::VerbParams vp) { vp_ = vp; }

 private:
  QUuid id_; //!< Unique ID for the case.
  int sim_time_sec_;
  int wic_time_sec_; //!< The number of seconds spent computing the well index for this case.

  string im_ = "", wm_ = "", em_ = "";
  string md_ = "Optimization";
  string cl_ = "Case";
  Settings::VerbParams vp_;

  double objective_function_value_;
  // QHash<QUuid, bool> binary_variables_;
  // QHash<QUuid, int> integer_variables_;
  // QHash<QUuid, double> real_variables_;

  // QMap<QUuid, bool> binary_variables_;
  // QMap<QUuid, int> integer_variables_;
  // QMap<QUuid, double> real_variables_;

  QList<QPair<QUuid, bool>> binary_variables_;
  QList<QPair<QUuid, int>> integer_variables_;
  QList<QPair<QUuid, double>> real_variables_;

  QList<QUuid> real_id_index_map_;
  QList<QUuid> integer_id_index_map_;
  QList<QUuid> binary_id_index_map_;

  Case* parent_; //!< The parent of this trial point. Needed by the APPS algorithm.
  int direction_index_; //!< The direction index used to generate this trial point.
  double step_length_; //!< The step length used to generate this trial point.

  // Multiple realizations-support
  QString ensemble_realization_; //!< The realization to evaluate next. Used by workers when in parallel mode.
  QHash<QString, double> ensemble_ofvs_; //!< Map of objective function values from realization alias - value.
};

}

#endif // CASE_H
