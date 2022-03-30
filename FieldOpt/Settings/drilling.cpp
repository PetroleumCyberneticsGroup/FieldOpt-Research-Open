/***********************************************************
Copyright (C) 2020-
Thiago L. Silva <thiagolims@gmail.com>
Mathias Bellout <chakibbb-pcg@gmail.com>

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

#include "drilling.h"
#include "settings_exceptions.h"

namespace Settings {

Drilling::Drilling(const QJsonObject& json_drilling, VerbParams vp) {
  vp_ = vp;
  readDrilling(json_drilling);
}

void Drilling::readDrilling(QJsonObject json_drilling) {

  if (json_drilling.contains("WellName")) {
    well_name = json_drilling["WellName"].toString();
  }

  if (json_drilling.contains("LocalOptimizer")) {
    QJsonObject json_local_optimizer = json_drilling["LocalOptimizer"].toObject();
    local_optimizer_settings = new Optimizer(json_local_optimizer, vp_);
  }

  if (json_drilling.contains("GlobalOptimizer")) {
    QJsonObject json_global_optimizer = json_drilling["GlobalOptimizer"].toObject();
    global_optimizer_settings = new Optimizer(json_global_optimizer, vp_);
  }

  if (!json_drilling.contains("DrillingSchedule"))
    throw UnableToParseDrillingSectionException(
      "The Drilling schedule must be defined at least one time in the Drilling/Model section.");

  if (json_drilling.contains("DrillingSchedule")) {
    QJsonArray json_d_sched = json_drilling["DrillingSchedule"].toArray();
    Drilling::DrillingSchedule d_sched;
    
    for (int i = 0; i < json_d_sched.size(); i++) {
      QJsonObject json_drilling_step = json_d_sched[i].toObject();
      d_sched.drilling_steps.append(i);
      if (json_drilling_step.contains("TimeStep")) {
        d_sched.time_steps.insert(i, json_drilling_step["TimeStep"].toDouble());
      }

      if (json_drilling_step.contains("Operation")) {
        if (QString::compare("StartDrilling", json_drilling_step["Operation"].toString()) == 0)
          d_sched.drilling_operations.insert(i, Drilling::DrillingSchedule::DrillingOperation::StartDrilling);

        if (QString::compare("Drilling", json_drilling_step["Operation"].toString()) == 0)
          d_sched.drilling_operations.insert(i, Drilling::DrillingSchedule::DrillingOperation::Drilling);

        if (QString::compare("PullingOutOfHole", json_drilling_step["Operation"].toString()) == 0)
          d_sched.drilling_operations.insert(i, Drilling::DrillingSchedule::DrillingOperation::PullingOutOfHole);
      }

      if (json_drilling_step.contains("Execution")) {
        if (QString::compare("Serial", json_drilling_step["Execution"].toString()) == 0)
          d_sched.execution_modes.insert(i, Drilling::DrillingSchedule::Execution::Serial);

        if (QString::compare("Parallel", json_drilling_step["Execution"].toString()) == 0)
          d_sched.execution_modes.insert(i, Drilling::DrillingSchedule::Execution ::Parallel);
      }

      if (json_drilling_step.contains("Optimizer")) {
        QJsonObject json_optimizer = json_drilling_step["Optimizer"].toObject();
        d_sched.optimizer_settings.insert(i, new Optimizer(json_optimizer, vp_));
      }

      if (json_drilling_step.contains("DrillingPoints")) {
        QJsonObject json_drilling_points = json_drilling_step["DrillingPoints"].toObject();

        if (!json_drilling_points.contains("x") || !json_drilling_points.contains("y")
          || !json_drilling_points.contains("z"))
          throw UnableToParseDrillingSectionException(
            "The Drilling points are not properly defined Drilling/Model section.");
        else {
          QJsonArray json_drilling_points_x = json_drilling_points["x"].toArray();
          QJsonArray json_drilling_points_y = json_drilling_points["y"].toArray();
          QJsonArray json_drilling_points_z = json_drilling_points["z"].toArray();

          if ((json_drilling_points_x.size() != json_drilling_points_y.size())
            || (json_drilling_points_y.size() != json_drilling_points_z.size()))
            throw UnableToParseDrillingSectionException(
              "The Drilling points have incompatible sizes in the different coordinates in the Drilling/Model section.");

          QList<Drilling::DrillingPoint> drilling_points;
          for (int j = 0; j < json_drilling_points_x.size(); j++) {
            Drilling::DrillingPoint point;
            point.x = json_drilling_points_x.at(j).toDouble();
            point.y = json_drilling_points_y.at(j).toDouble();
            point.z = json_drilling_points_z.at(j).toDouble();

            if (json_drilling_step.contains("OptimizeDrillingPoints")) {
              if (json_drilling_step["OptimizeDrillingPoints"].toBool())
                point.is_variable = true;
              else
                point.is_variable = false;
            }
            drilling_points.append(point);
          }
          d_sched.drilling_points.insert(i, drilling_points);
        }
      }

      if (json_drilling_step.contains("OptimizeDrillingPoints")) {
        if (json_drilling_step["OptimizeDrillingPoints"].toBool())
          d_sched.is_variable_drilling_points.insert(i, true);
        else
          d_sched.is_variable_drilling_points.insert(i, false);
      }

      if (json_drilling_step.contains("ModelUpdate")) {
        if (json_drilling_step["ModelUpdate"].toBool())
          d_sched.is_model_updates.insert(i, true);
        else
          d_sched.is_model_updates.insert(i, false);
      }

      if (json_drilling_step.contains("OptimizeCompletion")) {
        if (json_drilling_step["OptimizeCompletion"].toBool())
          d_sched.is_variable_completions.insert(i, true);
        else
          d_sched.is_variable_completions.insert(i, false);
      }

      if (json_drilling_step.contains("ModelType")) {
        if (json_drilling_step["ModelType"].toString().compare("TrueModel") == 0)
          d_sched.model_types.insert(i, Drilling::DrillingSchedule::ModelType::TrueModel);

        if (json_drilling_step["ModelType"].toString().compare("Surrogate") == 0)
          d_sched.model_types.insert(i, Drilling::DrillingSchedule::ModelType::Surrogate);
      }

      if (json_drilling_step.contains("WarmStart")) {
        d_sched.is_warm_start.insert(i, json_drilling_step["WarmStart"].toBool());
      }

      if (json_drilling_step.contains("Trigger")) {
        QJsonObject json_trigger = json_drilling_step["Trigger"].toObject();

        if (!json_trigger.contains("min_model_deviation") &&
          !json_trigger.contains("max_model_deviation")  &&
          !json_trigger.contains("max_obj_improvement"))
          throw UnableToParseDrillingSectionException(
            "The optimization triggers are not properly defined in the Drilling/Model section.");
        else {
          Drilling::OptimizationTrigger trigger;

          if (json_trigger.contains("max_model_deviation")) {
            trigger.max_model_deviation = json_trigger.value("max_model_deviation").toDouble();
          } else {
            trigger.max_model_deviation = -1;
          }

          if (json_trigger.contains("min_model_deviation")) {
            trigger.min_model_deviation = json_trigger.value("min_model_deviation").toDouble();
          } else {
            trigger.min_model_deviation = -1;
          }

          if (json_trigger.contains("max_obj_improvement")) {
            trigger.max_objective_improvement = json_trigger.value("max_obj_improvement").toDouble();
          } else {
            trigger.max_objective_improvement = 100;
          }

          d_sched.optimization_triggers.insert(i, trigger);
        }
      }
    }
    drilling_sched_ = d_sched;
  }
}

}
