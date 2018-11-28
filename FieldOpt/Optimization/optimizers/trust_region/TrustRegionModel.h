/******************************************************************************
   Created by thiagols on 27/11/18.
   Copyright (C) 2018 Thiago Lima Silva <thiagolims@gmail.com>

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

#ifndef FIELDOPT_TRUSTREGIONMODEL_H
#define FIELDOPT_TRUSTREGIONMODEL_H

#include <Settings/optimizer.h>
#include <Eigen/Core>

typedef struct Polynomial {
    int dimension;
    Eigen::VectorXd coefficients;
};

namespace Optimization {
namespace Optimizers {

/*!
 * @brief This class is an implementation of a Trust Region model.
 *
 * It delimits a region within a radius, with a set of points
 * and the corresponding function values. A polynomial model
 * is built to represent the function within the delimited region.
 */
class TrustRegionModel {
 public:
  /*!
   * @brief Initialize the trust region model.
   */
    TrustRegionModel();
    TrustRegionModel(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> initial_points, Eigen::VectorXd initial_fvalues , Settings::Optimizer *settings);

    void moveToBestPoint();
    void criticalityStep(Settings::Optimizer::Parameters settings);
    double checkInterpolation();
    bool rebuildModel(Settings::Optimizer::Parameters settings);
    bool improveModelNfp(Settings::Optimizer::Parameters settings);
    bool ensureImprovement(Settings::Optimizer::Parameters settings);
    bool isLambdaPoised(Settings::Optimizer::Parameters settings);
    bool changeTrCenter(Eigen::VectorXd new_point, Eigen::VectorXd fvalues, Settings::Optimizer::Parameters settings);
    std::map<Eigen::VectorXd, Eigen::VectorXd> solveTrSubproblem(Settings::Optimizer::Parameters settings);
    Eigen::Matrix<Polynomial, Eigen::Dynamic, 1>  computePolynomialModels();


 private:
    Polynomial pivot_polynomials;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> points_abs;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> points_shitfted;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> cached_points;
    Eigen::Matrix<Polynomial, Eigen::Dynamic, 1> modeling_polynomials;

    Eigen::VectorXd fvalues;
    Eigen::VectorXd tr_center;
    Eigen::RowVector2d pivot_values;

    double radius;
    int cache_max;
};

}
}


#endif //FIELDOPT_TRUSTREGIONMODEL_H
