/***********************************************************
Copyright (C) 2015-2017
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

#ifndef FIELDOPT_TEST_RESOURCE_VARIABLE_PROPERTY_LIST_H
#define FIELDOPT_TEST_RESOURCE_VARIABLE_PROPERTY_LIST_H

#include "Model/properties/var_prop_container.h"

using namespace Model::Properties;

namespace TestResources {
class TestResourceVariablePropertyContainer {

 public:
  TestResourceVariablePropertyContainer() {

    prod_heel_x_ = new ContinuousProperty(5.0);
    prod_heel_y_ = new ContinuousProperty(5.0);
    prod_heel_z_ = new ContinuousProperty(200.0);
    prod_toe_x_ = new ContinuousProperty(12);
    prod_toe_y_ = new ContinuousProperty(120.1);
    prod_toe_z_ = new ContinuousProperty(300.0);

    inje_heel_x_ = new ContinuousProperty(25.0);
    inje_heel_y_ = new ContinuousProperty(5.0);
    inje_heel_z_ = new ContinuousProperty(200.0);
    inje_toe_x_ = new ContinuousProperty(32);
    inje_toe_y_ = new ContinuousProperty(120.1);
    inje_toe_z_ = new ContinuousProperty(300.0);

    pseudocont_x_ = new ContinuousProperty(200);
    pseudocont_y_ = new ContinuousProperty(200);

    prod_heel_x_->setName("SplinePoint#TESTW#heel#x");
    prod_heel_y_->setName("SplinePoint#TESTW#heel#y");
    prod_heel_z_->setName("SplinePoint#TESTW#heel#z");
    prod_toe_x_->setName("SplinePoint#TESTW#toe#x");
    prod_toe_y_->setName("SplinePoint#TESTW#toe#y");
    prod_toe_z_->setName("SplinePoint#TESTW#toe#z");

    inje_heel_x_->setName("SplinePoint#INJE#heel#x");
    inje_heel_y_->setName("SplinePoint#INJE#heel#y");
    inje_heel_z_->setName("SplinePoint#INJE#heel#z");
    inje_toe_x_->setName("SplinePoint#INJE#toe#x");
    inje_toe_y_->setName("SplinePoint#INJE#toe#y");
    inje_toe_z_->setName("SplinePoint#INJE#toe#z");

    pseudocont_x_->setName("PseudoContVert#TESTW#x");
    pseudocont_y_->setName("PseudoContVert#TESTW#y");
    
    Settings::VerbParams vp_ = {};

    varcont_prod_spline_ = new VarPropContainer(vp_);
    varcont_prod_spline_->AddVariable(prod_heel_x_);
    varcont_prod_spline_->AddVariable(prod_heel_y_);
    varcont_prod_spline_->AddVariable(prod_heel_z_);
    varcont_prod_spline_->AddVariable(prod_toe_x_);
    varcont_prod_spline_->AddVariable(prod_toe_y_);
    varcont_prod_spline_->AddVariable(prod_toe_z_);

    varcont_two_spline_wells_ = new VarPropContainer(vp_);
    varcont_two_spline_wells_->AddVariable(prod_heel_x_);
    varcont_two_spline_wells_->AddVariable(prod_heel_y_);
    varcont_two_spline_wells_->AddVariable(prod_heel_z_);
    varcont_two_spline_wells_->AddVariable(prod_toe_x_);
    varcont_two_spline_wells_->AddVariable(prod_toe_y_);
    varcont_two_spline_wells_->AddVariable(prod_toe_z_);
    varcont_two_spline_wells_->AddVariable(inje_heel_x_);
    varcont_two_spline_wells_->AddVariable(inje_heel_y_);
    varcont_two_spline_wells_->AddVariable(inje_heel_z_);
    varcont_two_spline_wells_->AddVariable(inje_toe_x_);
    varcont_two_spline_wells_->AddVariable(inje_toe_y_);
    varcont_two_spline_wells_->AddVariable(inje_toe_z_);

    varcont_pseudocont_ = new VarPropContainer(vp_);
    varcont_pseudocont_->AddVariable(pseudocont_x_);
    varcont_pseudocont_->AddVariable(pseudocont_y_);

    prod_bhp_0_ = new ContinuousProperty(1.0);
    prod_bhp_10_ = new ContinuousProperty(1.0);
    prod_bhp_0_->setName("BHP#PRODUCER#0");
    prod_bhp_10_->setName("BHP#PRODUCER#10");

    varcont_prod_bhp_ = new VarPropContainer(vp_);
    varcont_prod_bhp_->AddVariable(prod_bhp_0_);
    varcont_prod_bhp_->AddVariable(prod_bhp_10_);
  }

  VarPropContainer *varcont_prod_spline_;
  VarPropContainer *varcont_pseudocont_;
  VarPropContainer *varcont_two_spline_wells_;
  VarPropContainer *varcont_prod_bhp_;
  VarPropContainer *varcont_6r_;

 private:
  ContinuousProperty *prod_heel_x_;
  ContinuousProperty *prod_heel_y_;
  ContinuousProperty *prod_heel_z_;
  ContinuousProperty *prod_toe_x_;
  ContinuousProperty *prod_toe_y_;
  ContinuousProperty *prod_toe_z_;
  ContinuousProperty *inje_heel_x_;
  ContinuousProperty *inje_heel_y_;
  ContinuousProperty *inje_heel_z_;
  ContinuousProperty *inje_toe_x_;
  ContinuousProperty *inje_toe_y_;
  ContinuousProperty *inje_toe_z_;

  ContinuousProperty *prod_bhp_0_;
  ContinuousProperty *prod_bhp_10_;

  ContinuousProperty *pseudocont_x_;
  ContinuousProperty *pseudocont_y_;

};
}

#endif //FIELDOPT_TEST_RESOURCE_VARIABLE_PROPERTY_LIST_H
