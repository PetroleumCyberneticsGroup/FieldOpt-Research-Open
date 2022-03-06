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

#ifndef WELL_H
#define WELL_H

#include "Settings/settings.h"
#include "Settings/model.h"
#include "Model/wells/wellbore/completions/completion.h"
#include "Model/wells/wellbore/completions/perforation.h"
#include "Model/properties/var_prop_container.h"
#include "Model/properties/discrete_property.h"
#include "Model/wells/control.h"
#include "Model/wells/wellbore/wellblock.h"
#include "Model/wells/wellbore/trajectory.h"
#include "Model/wells/wellbore/completions/icd.h"
#include "Reservoir/grid/eclgrid.h"
#include "WellIndexCalculation/wicalc_rixx.h"

#include "compartment.h"
#include "segment.h"

#include <QList>

namespace Model {
namespace Wells {

/*!
 * \brief The Well class represents any well in the model.
 */
class Well
{
 public:
  /*!
   * \brief Well Default constructor.
   * \param settings Settings object to create a well from.
   * \param well_number The index of the sepcific well
   * in the Model.Wells list to create a well from.
   * \param variables Variables object to add all new variables to.
   */
  Well(const Settings::Model& model_settings,
       int well_number,
       ::Model::Properties::VarPropContainer *variable_container,
       ::Reservoir::Grid::Grid *grid,
       ::Reservoir::WellIndexCalculation::wicalc_rixx *wic);

  struct Heel {
    int i = std::numeric_limits<int>::min();
    int j = std::numeric_limits<int>::min();
    int k = std::numeric_limits<int>::min();
  };

  template<typename T>
  std::vector<std::vector<T>> SplitVector(const std::vector<T>& vec, size_t n, bool dbg);

  enum PreferredPhase { Oil, Gas, Water, Liquid };

  QString name() const { return name_; }
  ::Settings::Model::WellType type() const { return type_; }
  QString group() const { return group_; }

  bool IsProducer();
  bool IsInjector();

  ::Settings::Model::PreferredPhase preferred_phase() const { return preferred_phase_; }
  double wellbore_radius() const { return wellbore_radius_->value(); }
  Wellbore::Trajectory *trajectory() { return trajectory_; }
  QList<Control *> *controls() { return controls_; }

  int heel_i() const { return heel_.i; }
  int heel_j() const { return heel_.j; }
  int heel_k() const { return heel_.k; }

  void Update();
  int GetTimeSpentInWIC() const { return trajectory_->GetTimeSpentInWic(); }

  bool HasSimpleICVs() const { return icds_.size() > 0; }
  vector<Wellbore::Completions::ICD> GetSimpleICDs() const { return icds_; }

  // Methods for segmented wells
  virtual bool IsSegmented() const { return is_segmented_; }
  std::vector<Compartment> GetCompartments() const;

  int GetNumCompartments() const { return compartments_.size(); };

  std::vector<Packer *> GetPackers() const;
  std::vector<ICD *> GetICDs() const;
  std::vector<Segment> GetSegments();
  std::vector<Segment> GetTubingSegments();
  std::vector<Segment> GetICDSegments();
  std::vector<Segment> GetAnnulusSegments();
  std::vector<int> GetICDSegmentIndices();

 protected:
  Settings::Model::Well well_settings_;
  QString name_;
  ::Settings::Model::WellType type_;
  QString group_;
  ::Settings::Model::PreferredPhase preferred_phase_;
  Properties::ContinuousProperty *wellbore_radius_;
  Wellbore::Trajectory *trajectory_;

  Eigen::Vector3d zero_pt_ = Eigen::Vector3d::Zero(3);

  QList<Wellbore::WellBlock *> *well_blocks_;

  string im_ = "", wm_ = "", em_ = "";
  string md_ = "Model";
  string cl_ = "Well";
  Settings::VerbParams vp_;


  //!< Whether the trajectory is defined. It does not
  //!< need to be for, e.g., control optimization.
  bool trajectory_defined_ = true;

  Heel heel_;
  QList<Control *> *controls_;

  // Fields for segmented wells
  bool is_segmented_ = false;
  void initializeSegmentedWell(Properties::VarPropContainer *variable_container);
  void initSegWellStruct(Properties::VarPropContainer *variable_container);

  double tub_diam_;            //!< Tubing (inner) diameter.
  double ann_diam_;            //!< Annular diameter.
  double icd_diam_;            //!< ICD diameter.

  double tub_cross_sect_area_; //!< Tubing cross section area.
  double ann_cross_sect_area_; //!< Annular cross section area.
  double icd_cross_sect_area_; //!< ICD cross section area.


  double tub_roughness_;       //!< Roughness for tubing segments.
  double ann_roughness_;       //!< Roughness for annulus segments.
  double icd_roughness_;       //!< Roughness for ICD.
  std::vector<Compartment> compartments_; //!< List of compartments.

  //!< List of icds for when we have neither a
  //!< defined compartmentalization or a trajectory.
  std::vector<Wellbore::Completions::ICD> icds_;

  // Methods for segmented wells
  std::vector<int> createTubingSegments(std::vector<Segment> &segments);
  std::vector<int> createICDSegments(std::vector<Segment> &segments, std::vector<int> &tubing_indexes) const;
  void createAnnulusSegments(std::vector<Segment> &segments, const std::vector<int> &icd_indexes);

  double length_delta_ = 2131.06288;
  double depth_delta_ = 2034.43714;

};

}
}

#endif // WELL_H
