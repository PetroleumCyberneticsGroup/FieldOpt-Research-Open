/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2021 Mathias Bellout
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

#include "oneoff_runner.h"

namespace Runner {

OneOffRunner::OneOffRunner(RuntimeSettings *runtime_settings)
  : AbstractRunner(runtime_settings) {
  std::cout << "Initializing one-off runner." << std::endl;
  InitializeSettings();
  InitializeModel();
  InitializeSimulator();
  //EvaluateBaseModel();
  InitObjF();
  InitBaseCase();
  //InitializeOptimizer();
  //InitializeBookkeeper();
  InitializeLogger();
}

void OneOffRunner::Execute() {
  applyWellPositionFromArguments();
  EvaluateBaseModel();
  // Write objective function value to file
  model_->wellCost(settings_->optimizer());
  //Utilities::FileHandling::WriteLineToFile(QString("%1").arg(objf_->value()),
  //                                         QString::fromStdString(settings_->paths().GetPath(Paths::OUTPUT_DIR))
  //                                         + "/f.out");
//    std::cout << objf_->value() << std::endl;
}

void OneOffRunner::applyWellPositionFromArguments() {
  QList<Model::Properties::ContinuousProperty *> prod_vars = model_->variables()->GetWSplineVars("PRODUCER");
  QList<Model::Properties::ContinuousProperty *> inje_vars = model_->variables()->GetWSplineVars("INJECTOR");
  if (prod_vars.length() != 6 || inje_vars.length() != 6)
    throw std::runtime_error("Found an incorrect number of variables for the wells.");

  QUuid phx, phy, phz, ihx, ihy, ihz, ptx, pty, ptz, itx, ity, itz;
  for (int i = 0; i < prod_vars.length(); ++i) {
    if (prod_vars[i]->propertyInfo().spline_end == Model::Properties::Property::SplineEnd::Heel) {
      if (prod_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::x)
        phx = prod_vars[i]->id();
      else if (prod_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::y)
        phy = prod_vars[i]->id();
      else if (prod_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::z)
        phz = prod_vars[i]->id();

    } else if (prod_vars[i]->propertyInfo().spline_end == Model::Properties::Property::SplineEnd::Toe) {
      if (prod_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::x)
        ptx = prod_vars[i]->id();
      else if (prod_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::y)
        pty = prod_vars[i]->id();
      else if (prod_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::z)
        ptz = prod_vars[i]->id();
    }
  }

  for (int i = 0; i < inje_vars.length(); ++i) {
    if (inje_vars[i]->propertyInfo().spline_end == Model::Properties::Property::SplineEnd::Heel) {
      if (inje_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::x)
        ihx = inje_vars[i]->id();
      if (inje_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::y)
        ihy = inje_vars[i]->id();
      if (inje_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::z)
        ihz = inje_vars[i]->id();
    }

    if (inje_vars[i]->propertyInfo().spline_end == Model::Properties::Property::SplineEnd::Toe) {
      if (inje_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::x)
        itx = inje_vars[i]->id();
      if (inje_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::y)
        ity = inje_vars[i]->id();
      if (inje_vars[i]->propertyInfo().coord == Model::Properties::Property::Coordinate::z)
        itz = inje_vars[i]->id();
    }
  }

  base_case_->set_real_variable_value(phx, rts_->prod_coords().first[0]);
  base_case_->set_real_variable_value(phy, rts_->prod_coords().first[1]);
  base_case_->set_real_variable_value(phz, rts_->prod_coords().first[2]);
  base_case_->set_real_variable_value(ptx, rts_->prod_coords().second[0]);
  base_case_->set_real_variable_value(pty, rts_->prod_coords().second[1]);
  base_case_->set_real_variable_value(ptz, rts_->prod_coords().second[2]);
  base_case_->set_real_variable_value(ihx, rts_->inje_coords().first[0]);
  base_case_->set_real_variable_value(ihy, rts_->inje_coords().first[1]);
  base_case_->set_real_variable_value(ihz, rts_->inje_coords().first[2]);
  base_case_->set_real_variable_value(itx, rts_->inje_coords().second[0]);
  base_case_->set_real_variable_value(ity, rts_->inje_coords().second[1]);
  base_case_->set_real_variable_value(itz, rts_->inje_coords().second[2]);

  auto constraint_handler_ =
    new Optimization::Constraints::ConstraintHandler(settings_->optimizer(),
                                                     model_->variables(),
                                                     model_->grid());
  constraint_handler_->SnapCaseToConstraints(base_case_);
  model_->ApplyCase(base_case_);
}

}
