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

TrustRegionModel::TrustRegionModel(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& initial_points, const Eigen::VectorXd& initial_fvalues, Settings::Optimizer *settings) {
    //TODO: implement this method
    settings_ = settings;

    points_abs_.setZero(initial_points.rows(), initial_points.cols());
    points_abs_ << initial_points;

    std::cout << points_abs_ << endl;

    fvalues_.setZero(initial_fvalues.size());
    fvalues_ << initial_fvalues;

    std::cout << fvalues_ << endl;


    radius_ = settings_->parameters().tr_initial_radius;
    tr_center_ = 0;
    dim_ = initial_points.cols();

    cache_max_ = 3*pow(dim_,2);

    rebuildModel();
    moveToBestPoint();
    computePolynomialModels();

    if (sizeof(points_abs_) < 2) {
        ensureImprovement();
    }
}

int TrustRegionModel::getDimension() {
    return dim_;
}

void TrustRegionModel::moveToBestPoint() {
    //TODO: implement this method
}

void TrustRegionModel::criticalityStep() {
    //TODO: implement this method
}

double TrustRegionModel::checkInterpolation() {
    //TODO: implement this method
}

void TrustRegionModel::rebuildModel() {
    //TODO: implement this method
}

void TrustRegionModel::improveModelNfp() {
    //TODO: implement this method
}

void TrustRegionModel::ensureImprovement() {
    //TODO: implement this method
}

bool TrustRegionModel::isLambdaPoised() {
    //TODO: implement this method
}

void TrustRegionModel::changeTrCenter(Eigen::VectorXd new_point, Eigen::VectorXd fvalues) {
    //TODO: implement this method
}

std::map<Eigen::VectorXd, Eigen::VectorXd> TrustRegionModel::solveTrSubproblem() {
    //TODO: implement this method
}

void TrustRegionModel::computePolynomialModels() {
    //TODO: implement this method
}

}
}