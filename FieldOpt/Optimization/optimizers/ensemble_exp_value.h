/******************************************************************************
 * Created: 04.03.2020 by Amanda
 *
 * This file is part of the FieldOpt project.
 *
 * Copyright (C) 2020- Amanda Souza Machado <amanda.automacaoufsc@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 *****************************************************************************/

#ifndef FIELDOPT_ENSEMBLE_EXP_VALUE_H
#define FIELDOPT_ENSEMBLE_EXP_VALUE_H

#include <boost/random.hpp>
#include "Optimization/case.h"
#include <chrono>

namespace Optimization {
namespace Optimizers {

/*!
 * The EnsembleHelper class contains facilities that helps the
 * Runner classes deal with ensembles (multiple realizations).
 */
class EnsembleExpValue {


    public:
        EnsembleExpValue();
        //EnsembleExpValue(const Settings::Ensemble &ensemble, int rng_seed=0);

        /*!
         * Set a new active case and populate the realization queue
         * for it. You will not be able to set a new active case
         * until all the realizations selected for this one have been
         * evaluated.
         */
        void Evaluate(Eigen::VectorXd incubemt_solution_);

        void Evaluate(int current_case_, Eigen::VectorXd incubemt_solution_);
        /*!
         * Check whether all the realizations selected for the currently
         * active case have been evaluated.
         */
        double ExpValue();

        /*!
         * Chech whether there are any cases available for evaluation.
         */

        /*!
        * The case that is currently being considered.
        */
        bool IsCaseDone();

        void addFunction(double (*function)(Eigen::VectorXd));

        double initialValue(Eigen::VectorXd x0);

        void Print(Eigen::VectorXd incubemt_solution_);

        std::vector<std::string> getPrint();

        std::vector<double> function_values_;

        double current_case_;
        double number_of_functions_;

        /*!
         * Realization queue, containing the alias' of all realizatons
         * that should still be evaluated for the current case.
         */
        std::vector<double (*)(Eigen::VectorXd)> rzn_queue_;

        /*!
         * List of alias' for realizations currently being evaluated.
         */


        std::vector<double (*)(Eigen::VectorXd)> vector_of_functions_;

        std::string var;

        std::vector<std::string> points;


};

}
}

#endif //FIELDOPT_ENSEMBLE_EXP_VALUE_H

