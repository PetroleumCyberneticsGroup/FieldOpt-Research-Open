/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb.pcg@gmail.com>

Modified 2019-2020 Thiago Lima Silva
<thiagolims@gmail.com>

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

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Settings/optimizer.h"

#include "case.h"
#include "case_handler.h"
#include "normalizer.h"
#include "constraints/constraint_handler.h"
#include "optimization_exceptions.h"
#include "Optimization/objective/objective.h"

#include "Model/model.h"
#include "Simulation/simulator_interfaces/simulator.h"
#include "Model/properties/var_prop_container.h"

#include "Runner/loggable.hpp"
#include "Runner/logger.h"
#include "Utilities/verbosity.h"
#include "Utilities/printer.hpp"

using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace Model {
class Model;
}

namespace Simulation {
class Simulator;
}

class Logger;

namespace Optimization {

class HybridOptimizer;

// =========================================================
// The Optimizer class is the abstract parent class for all
// optimizers. It is primarily designed to support direct
// search optimization algorithms.

class Optimizer : public Loggable
{
  friend class HybridOptimizer;

 public:
  Optimizer() = delete;

  /*!
   * \brief GetCaseForEvaluation Get a new, unevaluated case
   * for evaluation.
   *
   * If no unevaluated cases are currently available in the
   * CaseHandler, the iterate() method is called to generate
   * new cases.
   * \return Pointer to a new, unevaluated case.
   */
  Case *GetCaseForEvaluation();

  /*!
   * \brief SubmitEvaluatedCase Submit an already evaluated
   case to the optimizer.
   *
   * Submitted case is marked "recently evaluated" in CaseHandler.
   * \param c Case to submit.
   */
  void SubmitEvaluatedCase(Case *c);

  // -------------------------------------------------------
  /*!
   * \brief GetTentativeBestCase Get best case found so far.
   * \return
   */
  Case *GetTentativeBestCase() const;

  /*!
   * \brief case_handler Get the case handler.
   Used by the bookkeeper in the runner lib.
   */
  CaseHandler *case_handler() const { return case_handler_; }

  // Status related methods
  int nr_evaluated_cases() const {
    return case_handler_->EvaluatedCases().size();
  }

  int nr_queued_cases() const {
    return case_handler_->QueuedCases().size();
  }

  int nr_recently_evaluated_cases() const {
    return case_handler_->RecentlyEvaluatedCases().size();
  }

  // -------------------------------------------------------
  // TerminationCondition enum enumerates the reasons for
  // ending the opt runs. Returned by IsFinished method.
  enum TerminationCondition : int { NOT_FINISHED=0,
    MAX_EVALS_REACHED=1, MIN_STEP_LENGTH_REACHED=2,
    MAX_ITERS_REACHED=3, OPT_CRITERIA_REACHED=4,
    DFTR_MIN_RADIUS_REACHED=5,
    DFTR_CRIT_NORM_1ST_ORD_LT_TOLF=6,
    DFTR_CRIT_NORM_2ND_ORD_LT_TOLF=7, // placeholder
    DFTR_MAX_NUM_RHO_INF_MET=8
  };

  /*!
   * \brief IsFinished Check whether the optimization is
   * finished, i.e. if the the optimizer has reached some
   * termination condition.
   *
   * This method should be called before attempting to get
   * a new case for evaluation.
   * \return NOT_FINISHED (0, =false) if the optimization
   * has not finished, otherwise the non-zero reason for
   * termination.
   */
  virtual TerminationCondition IsFinished() = 0;

  //!< Get the CSV header for the status string.
  virtual QString GetStatusStringHeader() const;

  //!< Get CSV string describing the current optimizer state.
  virtual QString GetStatusString() const;

  //!< Enable writing a text log for the constraint operations.
  void EnableConstraintLogging(const QString& output_dir_path);

  //!< Check if the optimizer is asynchronous.
  bool IsAsync() const { return is_async_; }

  /*!
   * @brief Get the simulation duration in seconds for a case.
   * @param c Case to get simulation duration for.
   * @return Simulation duration in seconds. -1 if the case
   * has not been successfully simulated.
   */
  int GetSimulationDuration(Case *c);

 protected:
  /*!
   * \brief Base constructor for optimizers.
   * Initializes constraints and sets some member values.
   * \param opt_settings Settings for the optimizer.
   * \param base_case Base case for optimizer. Must already
   * have been evaluated (i.e., have an objective function value).
   *
   * \param variables The variable property container from
   * the Model (needed for initialization of constraints).
   * \param grid The grid to be used (needed for initialization
   * of some constraints).
   *
   * \param logger Logger object passed from runner.
   * \param case_handler CaseHandler object. This is passed
   * from the HybridOptimizer; defaults to 0, in which case
   * a new one will be created.
   * \param constraint_handler ConstraintHandler object.
   * This is passed from the HybridOptimizer; defaults to 0,
   * in which case a new one will be created.
   */
  Optimizer(::Settings::Optimizer *opt_settings,
            Case *base_case,
            Model::Properties::VarPropContainer *variables,
            Reservoir::Grid::Grid *grid,
            Logger *logger,
            CaseHandler *case_handler = nullptr,
            Constraints::ConstraintHandler *constraint_handler = nullptr);

  /*!
   * @brief Handle an incomming evaluated case. This is
   * called at the end of the SubmitEvaluatedCase method.
   * @param c
   */
  virtual void handleEvaluatedCase(Case *c) = 0;

  /*!
   * @brief Check whether the Case c is an improvement on
   * the tentative best case.
   * @param c Case to be checked.
   * @return True if improvement; otherwise false.
   */
  bool isImprovement(const Case* c);

  /*!
   * @brief Check if Case c1 is better than Case c2, taking
   * into account if we're maximizing or minimizing.
   */
  bool isBetter(const Case* c1, const Case *c2) const;

  /*!
   * \brief iterate Performs an iteration, generating
   * new cases and adding them to the case_handler.
   */
  virtual void iterate() = 0;

  LogTarget GetLogTarget() override;
  map<string, string> GetState() override;
  QUuid GetId() override;
  map<string, vector<double>> GetValues() override;

 public:

  class EmbeddedProblem {
   public:

    // EmbeddedProblem& operator=(const EmbeddedProblem& prob);
    static EmbeddedProblem* pProb_;

    static EmbeddedProblem& getReference() {
      if (pProb_ == 0) {
        static EmbeddedProblem prob_;
        pProb_ = &prob_;
        return prob_;
      } else {
        return *pProb_;
      }
    }

    // In code:
    // Optimization::Optimizer::EmbeddedProblem& prob =
    // Optimization::Optimizer::EmbeddedProblem::getReference();

   public:
    EmbeddedProblem() {};

    // Set methods
    void setProbName(string pn) { prob_name_ = pn; }
    void setNunVars(int n) { n_vars_ = n; }
    void setNumLinConst(int n_lc) { n_lconst_ = n_lc; }
    void setNunNnlConst(int n_nlc) { n_nlconst_ = n_nlc; }

    void setXInit(VectorXd x_init) { x_init_ = x_init; }
    void setXUb(VectorXd x_ub) { x_ub_ = x_ub; }
    void setXLb(VectorXd x_lb) { x_lb_ = x_lb; }
    void setFUb(VectorXd F_ub) { F_ub_ = F_ub; }
    void setFLb(VectorXd F_lb) { F_lb_ = F_lb; }

    void setTrC(double tr_c) { tr_c_ = tr_c; }
    void setTrG(VectorXd tr_g) { tr_g_ = tr_g; }
    void setTrH(MatrixXd tr_H) { tr_H_ = tr_H; }

    void setXSol(vector<double> xsol) { xsol_ = xsol; }
    void setFSol(vector<double> fsol) { fsol_ = fsol; }

    void setSNOPTExitCode(int ec) { SNOPT_exit_code_ = ec; }
    void setSNOPTErrorMsg(string em) { SNOPT_error_msg_ = em; }

    // Get methods
    string getProbName() { return prob_name_; }
    int getNunVars() { return n_vars_; }
    int getNumLinConst() { return n_lconst_; }
    int getNunNnlConst() { return n_nlconst_; }

    VectorXd getXInit() { return x_init_; }
    VectorXd getXUb() { return x_ub_; }
    VectorXd getXLb() { return x_lb_; }
    VectorXd getFUb() { return F_ub_; }
    VectorXd getFLb() { return F_lb_; }

    double getTrC() { return tr_c_; }
    VectorXd getTrG() { return tr_g_; }
    MatrixXd getTrH() { return tr_H_; }



    // TEST
    double getSqpC() { return sqp_c_; }
    VectorXd getSqpG() { return sqp_g_; }
    MatrixXd getSqpH() { return sqp_H_; }

    Case *tcase_ = nullptr;
    Model::Model *model_ = nullptr;
    Simulation::Simulator *simulator_ = nullptr;
    // Runner::EmbeddedRunner *runner_ = nullptr;



    VectorXd getXSol() {
      return Eigen::Map<VectorXd>(xsol_.data(), xsol_.size());
    }

    VectorXd getFSol() {
      return Eigen::Map<VectorXd>(fsol_.data(), fsol_.size());
    }

    int getSNOPTExitCode() { return SNOPT_exit_code_ ; }
    string getSNOPTErrorMsg() { return SNOPT_error_msg_; }

   private:
    // -----------------------------------------------------
    string prob_name_ = "";
    int n_vars_;
    int n_lconst_;
    int n_nlconst_;

    VectorXd x_init_;
    VectorXd x_ub_;
    VectorXd x_lb_;
    VectorXd F_ub_;
    VectorXd F_lb_;

    double tr_c_;
    VectorXd tr_g_;
    MatrixXd tr_H_;

    double sqp_c_;
    VectorXd sqp_g_;
    MatrixXd sqp_H_;

    vector<double> xsol_;
    vector<double> fsol_;
    int tr_exit_flag_;

    int SNOPT_exit_code_;
    string SNOPT_error_msg_;
  };








 protected:
  void updateTentativeBestCase(Case *c);

  // All cases (base case, unevaluated cases and evaluated
  // cases) passed to or generated by the optimizer.
  CaseHandler *case_handler_ = nullptr;

  // All constraints defined for the optimization.
  Constraints::ConstraintHandler *constraint_handler_ = nullptr;

  int evaluated_cases_;  // Number of evaluated cases.

  // Maximum number of objective function evaluations
  // allowed before terminating.
  int max_evaluations_;

  int iteration_;  // Current iteration number

  // The optimization mode, i.e. whether the objective
  // function should be maximized or minimized.
  ::Settings::Optimizer::OptimizerMode mode_;
  ::Settings::Optimizer::OptimizerType type_;

  // Inidcates whether or not the optimizer is asynchronous;
  // defaults to false.
  bool is_async_;

  Logger *logger_ = nullptr;

  // Whether logging should be performed; this should be
  // set to false when the optimizer is a component in
  // HybridOptimizer.
  bool enable_logging_;

  // Disable logging for this optimizer;
  // this is called by HybridOptimizer.
  void DisableLogging();


  // Switch for whether or not to use penalty
  // function to account for constraints.
  bool penalize_;

  // The best case encountered thus far.
  Case *tentative_best_case_;

  // The iteration in which the current tentative best case was found.
  int tentative_best_case_iteration_;

  // Indicates this object is a hybrid optimization component.
  bool is_hybrid_component_ = false;

  // Normalizer for objective function values.
  Normalizer normalizer_ofv_;

  // Initialize all normalization parameters.
  void initializeNormalizers();

  string im_ = "", wm_ = "", em_ = "";
  ::Settings::VerbParams vp_;
  string md_ = "Optimization";
  string cl_ = "Optimizer";

  class Summary : public Loggable {
   public:
    Summary(Optimizer *opt, TerminationCondition cond,
            map<string, string> ext_state=map<string, string>()) {
      opt_ = opt;
      cond_ = cond;
      ext_state_ = ext_state;
    }

    LogTarget GetLogTarget() override;
    map<string, string> GetState() override;
    QUuid GetId() override;
    map<string, vector<double>> GetValues() override;

   private:
    Optimizer *opt_;
    Optimizer::TerminationCondition cond_;
    map<string, string> ext_state_;
  };

  /*!
   * @brief Calculate penalized objf function value for case.
   * @param c Case to calculate the penalized objective function value for.
   * @return The penalized objective function value.
   */
  double PenalizedOFV(Case *c);

 private:
  QDateTime start_time_;

  // The number of seconds spent in the iterate() method.
  int seconds_spent_in_iterate_;

  /*!
   * @brief Initialize the OFV normalizer, setting the
   * parameters for it from the cases that have been
   * evaluated so far.
   */
  void initializeOfvNormalizer();
};

}

#endif // OPTIMIZER_H
