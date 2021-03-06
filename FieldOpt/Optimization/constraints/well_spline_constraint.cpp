/***********************************************************
Copyright (C) 2016-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
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

#include "well_spline_constraint.h"
#include "Utilities/verbosity.h"
#include "Utilities/printer.hpp"
#include <boost/lexical_cast.hpp>

namespace Optimization {
namespace Constraints {

using namespace Model::Properties;

WellSplineConstraint::Well WellSplineConstraint::
initWSplineConstraint(QList<Model::Properties::ContinuousProperty *> vars,
                      Settings::VerbParams vp) {
  xvp_ = vp;

  Well well;
  if ((vars.length() >= 6)
    && (vars.length() % 3) == 0
    && vars[0]->propertyInfo().prop_type == Property::PropertyType::SplinePoint) {

    if (xvp_.vOPT >= 3) {
      xim_ = "Using heel-toe parameterization for well spline constraint";
      ext_info(xim_, xmd_, xcl_, xvp_.lnw);
    }

    for (auto var : vars) {
      if (var->propertyInfo().spline_end == Property::SplineEnd::Heel) {

        if (var->propertyInfo().coord == Property::Coordinate::x) {
          well.heel.x = var->id();
        } else if (var->propertyInfo().coord == Property::Coordinate::y) {
          well.heel.y = var->id();
        } else if (var->propertyInfo().coord == Property::Coordinate::z) {
          well.heel.z = var->id();
        } else {
          xem_ = "Unable to parse variable " + var->name().toStdString();
          throw std::runtime_error(xem_);
        }

      } else if (var->propertyInfo().spline_end == Property::SplineEnd::Toe) {

        if (var->propertyInfo().coord == Property::Coordinate::x) {
          well.toe.x = var->id();
        } else if (var->propertyInfo().coord == Property::Coordinate::y) {
          well.toe.y = var->id();
        } else if (var->propertyInfo().coord == Property::Coordinate::z) {
          well.toe.z = var->id();
        } else {
          xem_ = "Unable to parse variable " + var->name().toStdString();
          throw std::runtime_error(xem_);
        }

      } else {
        xem_ = "Unable to parse variable " + var->name().toStdString();
        throw std::runtime_error(xem_);
      }
    }

    // Adding additional points
    if ((vars.length() > 6)
      && vars[0]->propertyInfo().prop_type == Property::PropertyType::SplinePoint) {
      if (xvp_.vOPT >= 3) {
        xim_ = "Using multi-point parameterization for well spline constraint";
        ext_info(xim_, xmd_, xcl_, xvp_.lnw);
      }

      std::map<int, Coord> addtl_points;
      for (auto var : vars) {
        // use hash/map when creating the ponts (it will autosort them alphabeticaly,
        // which is what we want)
        if (var->propertyInfo().spline_end == Property::SplineEnd::Middle) {
          if (var->propertyInfo().coord == Property::Coordinate::x) {
            addtl_points[var->propertyInfo().index].x = var->id();
          } else if (var->propertyInfo().coord == Property::Coordinate::y) {
            addtl_points[var->propertyInfo().index].y = var->id();
          } else if (var->propertyInfo().coord == Property::Coordinate::z) {
            addtl_points[var->propertyInfo().index].z = var->id();
          }
        }
      }

      for (int i = 0; i < addtl_points.size(); ++i) {
        well.additional_points.push_back(addtl_points[i+1]);
      }
    }

  } else if ((vars.length() == 3)
    && vars[0]->propertyInfo().prop_type == Property::PropertyType::SplinePoint) {
    for (auto var : vars) {
      if (var->propertyInfo().coord == Property::Coordinate::x) {
        well.toe.x = var->id();
      } else if (var->propertyInfo().coord == Property::Coordinate::y) {
        well.toe.y = var->id();
      } else if (var->propertyInfo().coord == Property::Coordinate::z) {
        well.toe.z = var->id();
      } else {
        xem_ = "Unable to parse variable " + var->name().toStdString();
        throw std::runtime_error(xem_);
      }
    }

  } else if ((vars.length() > 0)
    && vars[0]->propertyInfo().prop_type == Property::PropertyType::PolarSpline) {
    if (xvp_.vOPT >= 2) {
      xim_ = "Using PolarSpline parameterization for well spline constraint";
      ext_info(xim_, xmd_, xcl_, xvp_.lnw);
    }
    for (auto var : vars) {
      if (var->propertyInfo().polar_prop == Property::PolarProp::Midpoint) {
        if (var->propertyInfo().coord == Property::Coordinate::x) {
          well.midpoint.x = var->id();
        } else if (var->propertyInfo().coord == Property::Coordinate::y) {
          well.midpoint.y = var->id();
        } else if (var->propertyInfo().coord == Property::Coordinate::z) {
          well.midpoint.z = var->id();
        }
      }
    }

  } else {
    xem_  = "Incorrect number of variables (" + num2str(vars.length());
    xem_ += ") passed to the initialize well method.";
    throw std::runtime_error(xem_);
  }
  return well;
}

QPair<Eigen::Vector3d, Eigen::Vector3d>
WellSplineConstraint::GetEndpointValueVectors(Case *c,
                                              Well well) {

  double hx = c->get_real_variable_value(well.heel.x);
  double hy = c->get_real_variable_value(well.heel.y);
  double hz = c->get_real_variable_value(well.heel.z);

  double tx = c->get_real_variable_value(well.toe.x);
  double ty = c->get_real_variable_value(well.toe.y);
  double tz = c->get_real_variable_value(well.toe.z);

  Eigen::Vector3d heel(hx, hy, hz);
  Eigen::Vector3d toe(tx, ty, tz);
  return qMakePair(heel, toe);
}

std::vector<Eigen::Vector3d>
WellSplineConstraint::GetPointValueVectors(Case *c,
                                           WellSplineConstraint::Well well) {

  std::vector<Eigen::Vector3d> points;
  auto endpoints = GetEndpointValueVectors(c, well);
  points.push_back(endpoints.first);

  for (auto p : well.additional_points) {
    double x = c->get_real_variable_value(p.x);
    double y = c->get_real_variable_value(p.y);
    double z = c->get_real_variable_value(p.z);
    Eigen::Vector3d ep = Eigen::Vector3d(x, y, z);
    points.push_back(ep);
  }

  points.push_back(endpoints.second);
  return points;
}

double WellSplineConstraint::GetWellLength(Case *c,
                                           WellSplineConstraint::Well well) {
  double length = 0.0;
  auto points = GetPointValueVectors(c, well);
  for (int i = 0; i < points.size() - 1; ++i) {
    length += (points[i+1] - points[i]).norm();
  }
  return length;
}

}
}
