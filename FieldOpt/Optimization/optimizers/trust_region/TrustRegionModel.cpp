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

template <typename T>;
void sort_indexes(const std::vector<T> &v, const std::vector<size_t> idx) {
//    std::vector<size_t> idx(v.size()); //!<initialize original index locations>
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];}); //!<sort indexes based on comparing values in v>
//    return idx;
}

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
    auto dim = all_points.rows();
    auto n_points = all_points.cols();
    if (tr_center_ != 0) {
        all_points.col(0).swap(all_points.col(tr_center_)); //!<center will be first>
        auto center_fvalue = all_points(tr_center_);
        all_fvalues(tr_center_) = all_fvalues(0);
        all_fvalues(0) = center_fvalue;
    }
    //!<calculate distances from tr center>
    points_shitfted_.Zero(dim, n_points);
    distances_.resize(n_points);
    distances_.Zero();
    for (int i=1; i<n_points; i++) { //!<Shift all points to TR center>
        points_shitfted_.col(i) << points_abs_.col(i) - points_abs_.col(0); //!<Compute distances>
        distances_(i) = points_shitfted_.col(i).lpNorm<Eigen::Infinity>(); //<!distances in infinity or 2-norm>
    }
    //!<Reorder points based on their distances to the tr center>
    distances_ord_.resize(distances_.cols());
    distances_ord_ << distances_;
    index_vector_.setLinSpaced(distances_ord_.size(),0,distances_ord_.size()-1);
    sort_indexes(distances_ord_.data(), index_vector_.data());

    //TODO: reorder points_shifted_, points_abs_ and fvalues_ based on order in index_vector.
    //TODO: build polynomial model using the reordered points

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