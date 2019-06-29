/******************************************************************************
   Created by einar on 11/15/16.
   Copyright (C) 2016 Einar J.M. Baumann <einar.baumann@gmail.com>

   This file is part of the FieldOpt project.

   FieldOpt is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   FieldOpt is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FieldOpt.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/
#ifndef FIELDOPT_TEST_RESOURCE_TEST_FUNCTIONS_H
#define FIELDOPT_TEST_RESOURCE_TEST_FUNCTIONS_H

#include <Eigen/Core>

using namespace Eigen;

namespace TestResources {
    /*!
     * @brief This namespace contains synthetic functions meant to be used
     * for testQing optimzation algorithms.
     *
     * Note that these are all formulated for minimization.
     */
    namespace TestFunctions {

        /*!
         * @brief Sphere function.
         *
         * Formula: \f$ \sum_{i=1}^{n} x^2_i \f$
         * Minimum: \f$ f(x_1, ..., x_n) = f(0, ..., 0) = 0 \f$
         * Domain:  \f$ - \infty \leq x_i \leq \infty \f$
         * Dimensions: \f$ 1 \leq i \leq n \f$
         *
         * @param xs Vector of _continous_ variable values.
         * @return The function value at the given positon.
         */
        inline double Sphere(VectorXd xs) {
            return (xs.cwiseProduct(xs)).sum();
        }

        /*!
         * @brief Rosenbrock function.
         *
         * Formula: \f$ \sum_{i=1}^{n-1} \left[ 100 (x_{i+1} - x_i^2)^2 + (x_i -1)^2 \right] \f$
         * Minimum: \f$ f(1, ..., 1) = 0 \f$
         * Domain:  \f$ - \infty \leq x_i \leq \infty \f$
         * Dimensions: \f$ 1 \leq i \leq n \f$
         *
         * @param xs Vector of _continous_ variable values.
         * @return The function value at the given position.
         */
        inline double Rosenbrock(VectorXd xs) {
            VectorXd xhead = xs.head(xs.size() - 1);
            VectorXd xtail = xs.tail(xs.size() - 1);

            VectorXd p1 = xtail - xhead.cwiseProduct(xhead);
            VectorXd p2 = xhead - VectorXd::Ones(xhead.size());
            return (100 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
        }

        /*!
         * @brief CG.prob1 -> f = @(x) (1 - x(1))^2
         * Initial point: x0=[-1.2 2.0]
         */
        inline double tr_dfo_prob1(VectorXd xs) {
            double arg1 = 1 - xs(0);
            return pow(arg1,2);
        }

        /*!
         * @brief CG.prob2 -> f = @(x) log1p(x(1)^2) + x(2)^2;
         * Initial point: x0=[2.0 2.0]
         */
        inline double tr_dfo_prob2(VectorXd xs) {
            double arg1 = pow(xs(0),2);
            double arg2 = pow(xs(1),2);
            return log1p(arg1) + arg2;
        }

        /*!
         * @brief CG.prob3 -> f = @(x) sin(pi*x(1)/12) * cos(pi*x(2)/16);
         * Initial point: x0=[0.0 0.0]
         */
        inline double tr_dfo_prob3(VectorXd xs) {
            double arg1 = M_PI * xs(0)/12;
            double arg2 = M_PI * xs(1)/16;
            return sin(arg1) * cos(arg2);
        }

        /*!
         * @brief CG.prob4 -> f = @(x) 0.01*(x(1) - 1)^2 + (x(2) - x(1)^2)^2;
         * Initial point: x0=[2.0 2.0 2.0]
         */
        inline double tr_dfo_prob4(VectorXd xs) {
            double arg1 = 0.01 * pow((xs(0) - 1), 2);
            double arg2 = pow(xs(1) - pow(xs(0), 2), 2);
            return arg1 + arg2;
        }

        /*!
         * @brief CG.prob5 -> f = @(x) (x(1) - x(2))^2 + (x(2) - x(3))^4;
         * Initial point: x0=[-2.6 2.0 2.0]
         */
        inline double tr_dfo_prob5(VectorXd xs) {
            double arg1 = pow((xs(0) - xs(1)), 2);
            double arg2 = pow((xs(1) - xs(2)), 4);
            return arg1 + arg2;
        }

        /*!
         * @brief CG.prob6 -> f = @(x) (x(1) + x(2))^2 + (x(2) + x(3))^2;
         * Initial point: x0=[-4.0 1.0 1.0]
         */
        inline double tr_dfo_prob6(VectorXd xs) {
            double arg1 = pow((xs(0) - xs(1)), 2);
            double arg2 = pow((xs(1) - xs(2)), 2);
            return arg1 + arg2;
        }

        /*!
         * @brief CG.prob7 -> f = @(x) log1p(x(1)^2) + log1p((x(1)
         * - x(2))^2) + log1p((x(2) - x(3))^2) + log1p((x(3) - x(4))^2);
         * Initial point: x0=[2.0 2.0 2.0 2.0]
         */
        inline double tr_dfo_prob7(VectorXd xs) {
            double arg1 = log1p(pow(xs(0), 2));
            double arg2 = log1p(pow((xs(0) - xs(1)), 2));
            double arg3 = log1p(pow((xs(1) - xs(2)), 2));
            double arg4 = log1p(pow((xs(2) - xs(3)), 2));
            return arg1 + arg2 + arg3 + arg4;
        }

        /*!
         * @brief CG.prob8 -> f = @(x) (x(1)*x(2)*x(3)*x(4))^2;
         * Initial point: x0=[0.8 0.8 0.8 0.8]
         */
        inline double tr_dfo_prob8(VectorXd xs) {
            double arg1 = pow((xs(0) * xs(1) * xs(2) * xs(3)), 2);
            return arg1;
        }

        /*!
         * @brief CG.prob9 -> f = @(x) (x(1)-1)^2 + (x(2)-2)^2 + (x(3)-3)^2 + (x(4)-4)^2;
         * Initial point: x0=[1.0 1.0 1.0 1.0]
         */
        inline double tr_dfo_prob9(VectorXd xs) {
            double arg1 = pow((xs(0) - 1), 2);
            double arg2 = pow((xs(1) - 2), 2);
            double arg3 = pow((xs(2) - 3), 2);
            double arg4 = pow((xs(3) - 4), 2);
            return arg1 + arg2 + arg3 + arg4;
        }

        /*!
         * @brief CG.prob10 -> f = @(x)
         * (x(1) - x(2))^2 + (x(2) - x(3))^2 +
         * (x(3) - x(4))^4 + (x(4) - x(5))^4;
         * Initial point: x0=[2.0 sqrt(2) -1.0 2-sqrt(2) 0.5]
         */
        inline double tr_dfo_prob10(VectorXd xs) {
            double arg1 = pow((xs(0) - xs(1)), 2);
            double arg2 = pow((xs(1) - xs(2)), 2);
            double arg3 = pow((xs(2) - xs(3)), 4);
            double arg4 = pow((xs(3) - xs(4)), 4);
            return arg1 + arg2 + arg3 + arg4;
        }

        /*!
         * @brief CG.prob11 -> f = @(x) sum(2*x./(x.*x + 1));
         * Initial point: x0=[1.0 1.0 1.0 1.0]
         */
        inline double tr_dfo_prob11(VectorXd xs) {
            auto xc = xs; xc.fill(1);
            auto arg1 = 2*xs;
            auto arg2 = xs.cwiseProduct(xs) + xc;
            auto arg3 = arg1.cwiseQuotient(arg2);
            return arg3.sum();
        }

    }
}

#endif //FIELDOPT_TEST_RESOURCE_TEST_FUNCTIONS_H
