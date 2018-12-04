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
#include <numeric>
#include <functional>
#include <Settings/optimizer.h>
#include "Utilities/printer.hpp"
#include "Utilities/stringhelpers.hpp"
#include <Utilities/verbosity.h>

namespace Optimization {
namespace Optimizers {

bool compare(const int &lhs, const int &rhs, const double *distances) {
    return distances[lhs] <= distances[rhs];
}

TrustRegionModel::TrustRegionModel() {
}

TrustRegionModel::TrustRegionModel(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& initial_points, const Eigen::VectorXd& initial_fvalues, Settings::Optimizer *settings) {
    //TODO: implement this method
    settings_ = settings;

    points_abs_.setZero(initial_points.rows(), initial_points.cols());
    points_abs_ << initial_points;

    fvalues_.setZero(initial_fvalues.size());
    fvalues_ << initial_fvalues;

    radius_ =  settings_->parameters().tr_initial_radius;
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
    //Eigen::IOFormat frmt(3, 0, " ", "\n", "             [", "]");
    auto pivot_threshold = settings_->parameters().tr_pivot_threshold*std::fmin(1, radius_);
    Eigen::MatrixXd all_points(points_abs_.rows(), points_abs_.cols() + cached_points_.cols()); //!<All points we know>
    if (cached_points_.size() == 0) {
        all_points = points_abs_;
    } else {
        all_points << points_abs_;
        all_points << cached_points_;
    }
    Eigen::VectorXd all_fvalues(fvalues_.rows(), fvalues_.cols() + cached_fvalues_.cols()); //!<All function values we know>
    if (cached_fvalues_.size() == 0) {
        all_fvalues = fvalues_;
    } else {
        all_fvalues << fvalues_;
        all_fvalues << cached_fvalues_;
    }
    int dim = all_points.rows();
    int n_points = all_points.cols();
    if (tr_center_ != 0) {
        all_points.col(0).swap(all_points.col(tr_center_)); //!<center will be first>
        auto center_fvalue = all_points(tr_center_);
        all_fvalues(tr_center_) = all_fvalues(0);
        all_fvalues(0) = center_fvalue;
    }
    //!<calculate distances from tr center>
    points_shitfted_.resize(dim, n_points);
    distances_.resize(n_points);
    points_shitfted_.Zero(dim, n_points);
    distances_.Zero(n_points);
    for (int i=1; i<n_points; i++) { //!<Shift all points to TR center>
        points_shitfted_.col(i) << points_abs_.col(i) - points_abs_.col(0); //!<Compute distances>
        distances_(i) = points_shitfted_.col(i).lpNorm<Eigen::Infinity>(); //<!distances in infinity or 2-norm>
    }

    //!<dbg>
    //cout << "[          ] points_shifted_" << endl;
    //cout << points_shitfted_.format(frmt) << endl;
    //cout << "[          ] distances_" << endl;
    //cout << distances_.format(frmt) << endl;


    //!<Reorder points based on their distances to the tr center>
    index_vector_.setLinSpaced(distances_.size(),0,distances_.size()-1);
    std::sort(index_vector_.data(), index_vector_.data() + index_vector_.size(), std::bind(compare, std::placeholders::_1,  std::placeholders::_2, distances_.data()));
    sortVectorByIndex(distances_, index_vector_);
    sortVectorByIndex(fvalues_, index_vector_);

    //!<dbg>
    //cout << "[          ] distances_ord_" << endl;
    //cout << distances_.format(frmt) << endl;

    //cout << "[          ] fvalues_ord_" << endl;
    //cout << fvalues_.format(frmt) << endl;

    //cout << "[          ] index_vector_" << endl;
    //cout << index_vector_.format(frmt) << endl;

    sortMatrixByIndex(points_shitfted_, index_vector_);
    sortMatrixByIndex(points_abs_, index_vector_);

    //!<dbg>
    //cout << "[          ] points_shifted_ord" << endl;
    //cout << points_shitfted_.format(frmt) << endl;

    //cout << "[          ] points_abs_ord" << endl;
    //cout << points_abs_.format(frmt) << endl;

    //TODO: build polynomial model using the reordered points

}

void TrustRegionModel::sortVectorByIndex(Eigen::VectorXd &vec, const Eigen::VectorXd &ind) {
    Eigen::VectorXd vec_ord(vec.size());
    for (int i=0; i < vec.size(); i++) {
        int index = int(ind(i));
        vec_ord(i) = vec(index);
    }

    for(int i=0; i<vec.size(); i++) {
        vec(i) = vec_ord(i);
    }
    vec_ord.resize(0);
}

void TrustRegionModel::sortMatrixByIndex(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &points, const Eigen::VectorXd &ind) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> points_ord(points.rows(), points.cols());
    for (int i=0; i<points.cols(); i++) {
        int index = int(ind(i));
        points_ord.col(i) << points.col(index);
    }

    for (int i=0; i<points.cols(); i++) {
        points.col(i) << points_ord.col(i);
    }
    points_ord.resize(0,0);
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