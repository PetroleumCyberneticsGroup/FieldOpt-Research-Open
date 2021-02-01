/***********************************************************
Created by einar on 11/15/16.
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

// =========================================================
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
  return (100.0 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
}

inline double Rosenbrock_variation1(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v1 = xhead - 0.8*VectorXd::Ones(xhead.size());

  VectorXd p1 = xtail + 4*VectorXd::Ones(xhead.size()) - v1.cwiseProduct(v1);

  return (95 * p1.cwiseProduct(p1) + v1.cwiseProduct(v1)).sum();
}

inline double Rosenbrock_variation2(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v1 = xhead - 0.4*VectorXd::Ones(xhead.size());

  VectorXd p1 = xtail + 0.3*VectorXd::Ones(xhead.size()) - v1.cwiseProduct(v1);
  VectorXd p2 = xhead - VectorXd::Ones(xhead.size());
  return (97 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
}


inline double Rosenbrock_variation3(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v1 = xhead + 0.4*VectorXd::Ones(xhead.size());

  VectorXd p1 = xtail + 0.3*VectorXd::Ones(xhead.size()) - v1.cwiseProduct(v1);
  VectorXd p2 = xhead - 1.2*VectorXd::Ones(xhead.size());
  return (103 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2) + 1*VectorXd::Ones(xhead.size())).sum();
}


inline double Rosenbrock_variation4(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v1 = xhead - 0.3*VectorXd::Ones(xhead.size());

  VectorXd p1 = xtail - 1.8*VectorXd::Ones(xhead.size()) - v1.cwiseProduct(v1);
  VectorXd p2 = xhead + 0.8*VectorXd::Ones(xhead.size());
  return (94 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
}


inline double Rosenbrock_variation5(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v1 = xhead + 0.7*VectorXd::Ones(xhead.size());

  VectorXd p1 = xtail - v1.cwiseProduct(v1);
  VectorXd p2 = xhead + 0.3*VectorXd::Ones(xhead.size());
  return (98 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
}

inline double Rosenbrock_variation6(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v1 = xhead - 0.5*VectorXd::Ones(xhead.size());

  VectorXd p1 = xtail - 1.8*VectorXd::Ones(xhead.size())- v1.cwiseProduct(v1);
  VectorXd p2 = xhead - VectorXd::Ones(xhead.size());
  return (95 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
}

inline double Rosenbrock_variation7(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v1 = xhead - 0.7*VectorXd::Ones(xhead.size());

  VectorXd p1 = xtail - v1.cwiseProduct(v1);
  VectorXd p2 = xhead - 0.2*VectorXd::Ones(xhead.size());
  return (106 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2) ).sum();
}

inline double Rosenbrock_variation8(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);


  VectorXd p1 = xtail + 4*VectorXd::Ones(xhead.size()) - xhead.cwiseProduct(xhead);
  VectorXd p2 = xhead - 1.3*VectorXd::Ones(xhead.size());
  return (96 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2) ).sum();
}


inline double Rosenbrock_variation9(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);


  VectorXd p1 = xtail - 2*VectorXd::Ones(xhead.size()) - xhead.cwiseProduct(xhead);
  VectorXd p2 = xhead + 0.7*VectorXd::Ones(xhead.size());
  return (105 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
}


inline double Rosenbrock_variation10(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v1 = xhead - 0.2*VectorXd::Ones(xhead.size());

  VectorXd p1 = xtail + 0.6*VectorXd::Ones(xhead.size()) - v1.cwiseProduct(v1);
  VectorXd p2 = xhead - VectorXd::Ones(xhead.size());
  return (90 * p1.cwiseProduct(p1) + p2.cwiseProduct(p2)).sum();
}

inline double Rastingi(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v2 = (2*M_PI*xhead).array().cos();
  VectorXd v3 = (2*M_PI*xtail).array().cos();

  VectorXd p1 = xhead.cwiseProduct(xhead) - 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = xtail.cwiseProduct(xtail) - 10*v3*VectorXd::Ones(xhead.size());
  return (10*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation1(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd v2 = (2*M_PI*xhead).array().cos();
  VectorXd v3 = (2*M_PI*xtail).array().cos();

  VectorXd p1 = xhead.cwiseProduct(xhead) - 15*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = xtail.cwiseProduct(xtail) - 15*v3*VectorXd::Ones(xhead.size());
  return (15*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation2(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead - 0.2*VectorXd::Ones(xhead.size());
  VectorXd y = xtail;

  VectorXd v2 = (2*M_PI*x).array().cos();
  VectorXd v3 = (2*M_PI*xtail).array().cos();


  VectorXd p1 = x.cwiseProduct(x) - 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y) - 10*v3*VectorXd::Ones(xhead.size());
  return (10*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation3(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead;
  VectorXd y = xtail + 0.2*VectorXd::Ones(xhead.size());

  VectorXd v2 = (2*M_PI*x+0.4*VectorXd::Ones(xhead.size())).array().cos();
  VectorXd v3 = (2*M_PI*y).array().cos();


  VectorXd p1 = x.cwiseProduct(x) - 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y) - 10*v3*VectorXd::Ones(xhead.size());
  return (10*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 +0.5*VectorXd::Ones(xhead.size())).sum();
}

inline double Rastingi_varation4(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead +0.6*VectorXd::Ones(xhead.size());
  VectorXd y = xtail + 0.15*VectorXd::Ones(xhead.size());

  VectorXd v2 = (2*M_PI*x).array().cos();
  VectorXd v3 = (2*M_PI*y).array().cos();


  VectorXd p1 = xhead.cwiseProduct(xhead) - 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y) - 10*v3*VectorXd::Ones(xhead.size());
  return (11*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation5(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead -0.2*VectorXd::Ones(xhead.size());
  VectorXd y = xtail;

  VectorXd v2 = (2*M_PI*x).array().cos();
  VectorXd v3 = (2*M_PI*y).array().cos();


  VectorXd p1 = x.cwiseProduct(x) - 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y) - 0.5*VectorXd::Ones(xhead.size()) - 10*v3*VectorXd::Ones(xhead.size());
  return (15*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation6(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead +0.3*VectorXd::Ones(xhead.size());
  VectorXd y = xtail -0.1*VectorXd::Ones(xhead.size());

  VectorXd v2 = (2*M_PI*x).array().cos();
  VectorXd v3 = (2*M_PI*xtail).array().cos();


  VectorXd p1 = x.cwiseProduct(x) - 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y)  - 10*v3*VectorXd::Ones(xhead.size());
  return (14*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation7(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead -0.4*VectorXd::Ones(xhead.size());
  VectorXd y = xtail -0.4*VectorXd::Ones(xhead.size());

  VectorXd v2 = (2*M_PI*x).array().cos();
  VectorXd v3 = (2*M_PI*y).array().cos();


  VectorXd p1 = x.cwiseProduct(x) + 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y)  - 10*v3*VectorXd::Ones(xhead.size());
  return (5*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation8(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead +0.2*VectorXd::Ones(xhead.size());
  VectorXd y = xtail -0.2*VectorXd::Ones(xhead.size());

  VectorXd v2 = (2*M_PI*xhead).array().cos();
  VectorXd v3 = (2*M_PI*xtail).array().cos();


  VectorXd p1 = x.cwiseProduct(x) - 15*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y)  -0.2*VectorXd::Ones(xhead.size()) - 15*v3*VectorXd::Ones(xhead.size());
  return (14*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation9(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead +0.6*VectorXd::Ones(xhead.size());
  VectorXd y = xtail +0.1*VectorXd::Ones(xhead.size());
  VectorXd y1 = xtail -0.2*VectorXd::Ones(xhead.size());

  VectorXd v2 = (2*M_PI*x).array().cos();
  VectorXd v3 = (2*M_PI*y1).array().cos();


  VectorXd p1 = xhead.cwiseProduct(xhead) - 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y)  -10*v3*VectorXd::Ones(xhead.size());
  return (13*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}

inline double Rastingi_varation10(VectorXd xs) {
  VectorXd xhead = xs.head(xs.size() - 1);
  VectorXd xtail = xs.tail(xs.size() - 1);

  VectorXd x = xhead -0.6*VectorXd::Ones(xhead.size());
  VectorXd x1 = xhead +1*VectorXd::Ones(xhead.size());
  VectorXd y = xtail -0.5*VectorXd::Ones(xhead.size());

  VectorXd v2 = (2*M_PI*x1).array().cos();
  VectorXd v3 = (2*M_PI*xtail).array().cos();


  VectorXd p1 = xhead.cwiseProduct(xhead) - 10*v2*VectorXd::Ones(xhead.size());
  VectorXd p2 = y.cwiseProduct(y)  -10*v3*VectorXd::Ones(xhead.size());
  return (12*xs.size()*VectorXd::Ones(xhead.size()) + p1 + p2 ).sum();
}
/*!
 * @brief CG.prob1 -> f = @(x) (1 - x(0))^2
 * Initial point: x0=[-1.2 2.0]
 */
inline double tr_dfo_prob1(VectorXd xs) {
  double arg1 = 1 - xs(0);
  return pow(arg1,2);
}

inline double tr_dfo_prob1_test(VectorXd xs) {
  double arg1 = 2 - xs(0);
  return pow(arg1, 2);
}

/*!
 * @brief CG.prob2 -> f = @(x) log1p(x(0)^2) + x(1)^2;
 * Initial point: x0=[2.0 2.0]
 */
inline double tr_dfo_prob2(VectorXd xs) {
  double arg1 = pow(xs(0), 2);
  double arg2 = pow(xs(1), 2);
  return log1p(arg1) + arg2;
}

inline double tr_dfo_prob2_variation1(VectorXd xs) {
  double arg1 = pow(xs(0),2);
  double arg2 = pow(xs(1),2);
  return 1.5*log1p(4*arg1) + 4*arg2;
}
inline double tr_dfo_prob2_variation2(VectorXd xs) {
  double arg1 = pow(xs(0),2);
  double arg2 = pow(xs(1),2);
  return 2*log1p(3*arg1) + arg2;
}
inline double tr_dfo_prob2_variation3(VectorXd xs) {
  double arg1 = pow(xs(0),2);
  double arg2 = pow(xs(1),2);
  return log1p(2*arg1) + 2*arg2 + 0.01;
}

/*!
 * @brief CG.prob3 -> f = @(x) sin(pi*x(0)/12) * cos(pi*x(1)/16);
 * Initial point: x0=[0.0 0.0]
 */
inline double tr_dfo_prob3(VectorXd xs) {
  double arg1 = M_PI * xs(0)/12;
  double arg2 = M_PI * xs(1)/16;
  return sin(arg1) * cos(arg2);
}

/*!
 * @brief CG.prob4 -> f = @(x) 0.01*(x(0) - 1)^2 + (x(1) - x(0)^2)^2;
 * Initial point: x0=[2.0 2.0 2.0]
 */
inline double tr_dfo_prob4(VectorXd xs) {
  double arg1 = 0.01 * pow((xs(0) - 1), 2);
  double arg2 = pow(xs(1) - pow(xs(0), 2), 2);
  return arg1 + arg2;
}

inline double tr_dfo_prob4_variation1(VectorXd xs) {
  double arg1 = 0.01 * pow((xs(0) + 10), 2);
  double arg2 = pow((xs(1)-10) - pow((xs(0)+10), 2), 2);
  return arg1 + arg2;
}

inline double tr_dfo_prob4_variation2(VectorXd xs) {
  double arg1 = 0.01 * pow((xs(0) - 10), 2);
  double arg2 = pow(10.0 * xs(1) - pow((xs(0)-10), 2), 2);
  return arg1 + arg2;
}

inline double tr_dfo_prob4_variation3(VectorXd xs) {
  double arg1 = 0.01 * pow((xs(0) - 1), 2);
  double arg2 = pow(xs(1) - pow((xs(0)+4), 2), 2);
  return arg1 + arg2;
}

/*!
 * @brief CG.prob5 -> f = @(x) (x(0) - x(1))^2 + (x(1) - x(2))^4;
 * Initial point: x0=[-2.6 2.0 2.0]
 */
inline double tr_dfo_prob5(VectorXd xs) {
  double arg1 = pow((xs(0) - xs(1)), 2);
  double arg2 = pow((xs(1) - xs(2)), 4);
  return arg1 + arg2;
}

/*!
 * @brief CG.prob6 -> f = @(x) (x(0) + x(1))^2 + (x(1) + x(2))^2;
 * Initial point: x0=[-4.0 1.0 1.0]
 */
inline double tr_dfo_prob6(VectorXd xs) {
  double arg1 = pow((xs(0) - xs(1)), 2);
  double arg2 = pow((xs(1) - xs(2)), 2);
  return arg1 + arg2;
}

/*!
   * @brief CG.prob7 -> f = @(x) log1p(x(0)^2) + log1p((x(0)
   * - x(1))^2) + log1p((x(1) - x(2))^2) + log1p((x(2) - x(3))^2);
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
 * @brief CG.prob8 -> f = @(x) (x(0)*x(1)*x(2)*x(3))^2;
 * Initial point: x0=[0.8 0.8 0.8 0.8]
 */
inline double tr_dfo_prob8(VectorXd xs) {
  double arg1 = pow((xs(0) * xs(1) * xs(2) * xs(3)), 2);
  return arg1;
}

/*!
 * @brief CG.prob9 -> f = @(x) (x(0)-1)^2 + (x(1)-2)^2 + (x(2)-3)^2 + (x(3)-4)^2;
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
 * (x(0) - x(1))^2 + (x(1) - x(2))^2 +
 * (x(2) - x(3))^4 + (x(3) - x(4))^4;
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

/*
 * Hock and Schittkowski problems from
 * https://apmonitor.com/wiki/index.php/Apps/HockSchittkowski
 */

/*
  HS1
  Initial point x0 = (-2.0, 1.0);
  Solution: 0.0
*/
inline double hs1(VectorXd x) {
  return 100*pow(x(1) - pow(x(0), 2), 2) + pow(1 - x(0), 2);
}

/*
  HS2
  Initial point x0 = (-2.0, 1.0);
  Bounds: ub = (inf, 1.5);
  Solution: 0.0
*/
inline double hs2(VectorXd x){
  return hs1(x);
}

/*
  HS3
  Initial point x0 = (10.0, 1.0);
  Bounds: lb = (-inf, 0.0);
  Solution: 0.0
*/
inline double hs3(VectorXd x){
  return x(1) + 1e-5*pow(x(1)-x(0), 2);
}


/*
  HS4
  Initial point x0 = (1.125, 0.125)
  Bounds: lb = (1.0, 0.0)
  Solution: 8.0/3.0
*/
inline double hs4(VectorXd x){
  return pow(x(0) + 1, 3)/3 + x(1);
}

/*
  HS5
  Initial point x0 = (0.0, 0.0)
  Bounds: lb = (-1.5, -3), ub = (4, 3)
  Solution: -(sqrt(3)/2 + 3.14159/3) = -1.91322207
*/
inline double hs5(VectorXd x){
  return sin(x(0) + x(1)) + pow(x(0) - x(1), 2) - 1.5*x(0) + 2.5*x(1) + 1;
}

/*
  HS25
  Initial point x0 = (100.0, 12.5, 3.0)
  Bounds: lb = (0.1, 0, 0), ub = (100.0, 25.6, 5.0);
  Solution: 0.0
*/
inline double hs25(VectorXd x){
  double result = 0;
  double u;

  for (double i = 1; i < 100; i++){
    u = 25.0 + pow(-50.0*log(i/100), (2/3));
    result = result + pow(exp(-u - pow(x(1), x(2))/x(0)) - i/100, 2);
  }
  return result;
}

/*
  HS38
  Initial point x0 = (-3.0, -1.0, -3.0, -1)
  Bounds: lb = (-10.0, -10.0, -10.0, -10.0), ub = (10.0, 10.0, 10.0, 10.0)
  Solution: 0
*/
inline double hs38(VectorXd x){
  return 100*pow(x(1) - pow(x(0), 2), 2)
      + pow(1 - x(0), 2)
      + 90*pow(x(3) - pow(x(2), 2), 2)
      + pow(1 - x(2), 2)
      + 10.1*(pow(x(1) - 1, 2)
          + pow(x(3) - 1, 2))
      + 19.8*(x(1) - 1)*(x(3) - 1);
}

/*
  HS45
  Initial point x0 = (0.0, 0.0, 0.0, 0.0, 0.0)
  Bounds: lb = (0.0, 0.0, 0.0, 0.0, 0.0), ub = (1.0, 2.0, 3.0, 4.0, 5.0);
  Solution: 1
*/
inline double hs45(VectorXd x){
  return 2 - x(0)*x(1)*x(2)*x(3)*x(4)/120;
}



}
}

#endif //FIELDOPT_TEST_RESOURCE_TEST_FUNCTIONS_H
