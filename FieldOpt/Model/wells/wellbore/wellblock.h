/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb.pcg@gmail.com>

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

#ifndef WELLBLOCK_H
#define WELLBLOCK_H

#include "Model/properties/discrete_property.h"
#include "Model/wells/wellbore/completions/completion.h"
#include "Model/wells/wellbore/completions/perforation.h"
#include <Eigen/Core>

namespace Model {
namespace Wells {
namespace Wellbore {

using Printer::num2str;

/*!
 * \brief The WellBlock class represents a single well block.
 * It contains references to any completion defined within it.
 *
 * Note that wellblocks should use 1-indexed IJK-indices,
 * instead of 0-indexing as is used in the rest of FieldOpt
 * and ERT. This is because the IJK-values will be used
 * directly by the simulator interface.
 */
class WellBlock
{
  friend class Trajectory;
 public:
  WellBlock(int i, int j, int k);

  virtual ~WellBlock();

/*!
 * \brief The DirOfPenetration enum indicates the
 * wells direction of penetration through block. The W
 * value should be used to indicate that the direction
 * could not be calculated.
 */
  enum DirOfPenetration : int { X=1, Y=2, Z=3, W=4, D=5 };

  //!< Add a completion to this well block.
  void AddCompletion(::Model::Wells::Wellbore::Completions::Completion *completion);

  bool HasCompletion(); //!< Check if this well block has a completion.
  Completions::Completion *GetCompletion(); //!< Get the completion defined for this block.

  bool HasPerforation(); //!< Check if this well block has a perforation-type completion.
  Completions::Perforation *GetPerforation(); //!< Get the perforation defined for this block.

  int i() const { return i_->value(); }
  int j() const { return j_->value(); }
  int k() const { return k_->value(); }

  DirOfPenetration directionOfPenetration() const {
    return dir_of_penetration_;
  }

  void setI(const int i) { i_->setValue(i); }
  void setJ(const int j) { j_->setValue(j); }
  void setK(const int k) { k_->setValue(k); }

  void setDirOfPenetration(const DirOfPenetration dop) {
    dir_of_penetration_ = dop;
  }

  string getDirPenetrationStr() {
    if (dir_of_penetration_ == X) return "X";
    if (dir_of_penetration_ == Y) return "Y";
    if (dir_of_penetration_ == Z) return "Z";
    if (dir_of_penetration_ == D) return "D";
    if (dir_of_penetration_ == W) return "W";
    return "N";
  }

  void setEntryPoint(const Eigen::Vector3d& entry_point) { entry_point_ = entry_point; }
  void setExitPoint(const Eigen::Vector3d& exit_point) { exit_point_ = exit_point; }

  Eigen::Vector3d getEntryPoint() const { return entry_point_; }
  Eigen::Vector3d getExitPoint() const { return exit_point_; }

  void setEntryMd(const double entry_md) { entry_md_ = entry_md; }
  void setExitMd(const double exit_md) { exit_md_ = exit_md; }
  void setLength(const double length) { length_ = length; }
  void setDxDyDz(const Eigen::VectorXd dxdydz) { dxdydz_ = dxdydz; }

  double getEntryMd() const { return entry_md_; }
  double getExitMd() const { return exit_md_; }
  double getLength() const { return length_; }

  double getDepthChange() const { return exit_point_[2] - entry_point_[2]; }
  double getDepth() const { return entry_point_[2] + (exit_point_[2] - entry_point_[2])/2; }
  Eigen::VectorXd getDxDyDz() { return dxdydz_; }

  string getPropString(double md=0) {
    stringstream ss;
    if (md > 0) {
      ss << "md:[" << num2str(md, 3, 1, 7)  << "] -- ";
    }
    ss << "block:[";
    ss << num2str(this->i(), 0, 0, 4) + ",";
    ss << num2str(this->j(), 0, 0, 4) + ",";
    ss << num2str(this->k(), 0, 0, 4) + "] ";
    ss << "entry_md:[" << num2str(this->entry_md_, 3, 1, 7)  << "] ";
    ss << "exit_md:[" << num2str(this->exit_md_, 3, 1, 7)  << "] ";
    return ss.str();
  }

 private:
  Model::Properties::DiscreteProperty *i_;
  Model::Properties::DiscreteProperty *j_;
  Model::Properties::DiscreteProperty *k_;

  Eigen::Vector3d entry_point_ = Eigen::Vector3d::Zero(); //!< Entry point for splines through this block.
  Eigen::Vector3d exit_point_ = Eigen::Vector3d::Zero();  //!< Exit point for splines through this block.

  double entry_md_ = 0.0;
  double exit_md_ = 0.0;
  double length_ = 0.0;
  Completions::Completion *completion_;

  Eigen::VectorXd dxdydz_;

  //!< Well's direction of penetration through block.
  DirOfPenetration dir_of_penetration_;
};

}
}
}
#endif // WELLBLOCK_H
