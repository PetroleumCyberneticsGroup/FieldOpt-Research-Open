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
#include "TrustRegionModel.h"
#include <Settings/optimizer.h>

namespace Optimization {
namespace Optimizers {

TrustRegionModel::TrustRegionModel() {
}

TrustRegionModel::TrustRegionModel(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> initial_points, Eigen::VectorXd initial_fvalues, Settings::Optimizer *settings) {
    //TODO: implement this method
}

void TrustRegionModel::moveToBestPoint() {
    //TODO: implement this method
}

void TrustRegionModel::criticalityStep(Settings::Optimizer::Parameters settings) {
    //TODO: implement this method
}

double TrustRegionModel::checkInterpolation() {
    //TODO: implement this method
}

bool TrustRegionModel::rebuildModel(Settings::Optimizer::Parameters settings) {
    //TODO: implement this method
}

bool TrustRegionModel::improveModelNfp(Settings::Optimizer::Parameters settings) {
    //TODO: implement this method
}

bool TrustRegionModel::ensureImprovement(Settings::Optimizer::Parameters settings) {
    //TODO: implement this method
}

bool TrustRegionModel::isLambdaPoised(Settings::Optimizer::Parameters settings) {
    //TODO: implement this method
}

bool TrustRegionModel::changeTrCenter(Eigen::VectorXd new_point, Eigen::VectorXd fvalues,
                                      Settings::Optimizer::Parameters settings) {
    //TODO: implement this method
}

std::map<Eigen::VectorXd, Eigen::VectorXd> TrustRegionModel::solveTrSubproblem(Settings::Optimizer::Parameters settings) {
    //TODO: implement this method
}

Eigen::Matrix<Polynomial, Eigen::Dynamic, 1> TrustRegionModel::computePolynomialModels() {
    //TODO: implement this method
}

}
}