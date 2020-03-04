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

#include <iostream>
#include "ensemble_exp_value.h"
#include "Utilities/math.hpp"
#include <numeric>
#include <string>
#include <sstream>

namespace Optimization {
namespace Optimizers {

        EnsembleExpValue::EnsembleExpValue() {
            current_case_ = 0.0;
            number_of_functions_ = 0.0;
            std::vector<double (*)(Eigen::VectorXd)> vector_of_functions_;
            std::vector<double> function_values_;
            std::string var = " ";
            std::vector<std::string> points;


        }

        void EnsembleExpValue::Evaluate(Eigen::VectorXd incubemt_solution_) {

            for (int i = 0; i < number_of_functions_; ++i){
                Evaluate(i, incubemt_solution_);
            }
            Print(incubemt_solution_);
        }

        void EnsembleExpValue::Print(Eigen::VectorXd incubemt_solution_){
            std::ostringstream oss;
            oss << "[ [" << incubemt_solution_(0) << ", " << incubemt_solution_(1) << "]; [";
            var = oss.str();

            for (int element = 0; element < number_of_functions_; ++element){
                var += std::to_string(function_values_[element]);
                if(element < (number_of_functions_-1)){
                    var += ", ";
                }
            }
            var += "]; ";

        }

        void EnsembleExpValue::Evaluate(int current_case_, Eigen::VectorXd incubemt_solution_){
            double value = (vector_of_functions_[current_case_])(incubemt_solution_);
            function_values_.insert(function_values_.begin() + current_case_, value);
        }

        double EnsembleExpValue::ExpValue() {
            // Debug
            //cout << "Function size: ";
            //cout << function_values_.size() << endl;
            //cout << "Accumulative :";
            //cout << accumulate(function_values_.begin(), function_values_.end(), 0.0) << endl;

            double value = (accumulate(function_values_.begin(), function_values_.end(), 0.0))/number_of_functions_;

            // Cleaning elements of Ensemble for the next case
            function_values_.clear();
            current_case_ = 0.0;


            var += std::to_string(value); // Debug
            var += "]";
            points.insert(points.end(), var);
            var = " ";


            return value;

        }

        bool EnsembleExpValue::IsCaseDone() {
            if (function_values_.size() == number_of_functions_) {
                return true;
            } else {
                return false;
            }
        }

        void EnsembleExpValue::addFunction(double (*function)(Eigen::VectorXd)) {

            vector_of_functions_.push_back(function);
            current_case_ = 0.0;
            ++number_of_functions_;

        }

        double EnsembleExpValue::initialValue(Eigen::VectorXd x0){
            for (int i = 0; i < number_of_functions_; ++i){
                double value = (vector_of_functions_[i])(x0);
                function_values_.insert(function_values_.begin() + current_case_, value);
            }

            // Just to visualize the intial value
            double value = ExpValue();
            //cout << "Initial value:";
            //cout << value << endl;

            return value;
        }

        std::vector<std::string> EnsembleExpValue::getPrint(){
            //for (int element; element <= points.size(); ++element){
            //    cout << points[element] << endl;
            //}
            return points;
        }


}
}
