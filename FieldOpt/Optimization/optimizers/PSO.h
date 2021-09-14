/***********************************************************
Copyright (C) 2018
Brage Strand Kristoffersen <brage_sk@hotmail.com>
Created by Brage on 08/11/18.

Modified 2021 Mathias Bellout
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

#include <boost/random.hpp>
#include "optimizer.h"

#ifndef FIELDOPT_PSO_H
#define FIELDOPT_PSO_H
namespace Optimization {
namespace Optimizers {

/*!
 * @brief Class implements Particle Swarm Optimization (PSO)
 *
 * Implementation based on the description found at:
 * http://www.cleveralgorithms.com/nature-inspired/swarm/pso.html
 * https://github.com/clever-algorithms/CleverAlgorithms
 */
class PSO : public Optimizer {
 public:
  PSO(Settings::Optimizer *settings,
      Case *base_case,
      Model::Properties::VarPropContainer *variables,
      Reservoir::Grid::Grid *grid,
      Logger *logger,
      CaseHandler *case_handler=nullptr,
      Constraints::ConstraintHandler *constraint_handler=nullptr);

 protected:
  void handleEvaluatedCase(Case *c) override;
  void iterate() override;
  TerminationCondition IsFinished() override;

 protected:
  //!< Random number generator, random functions in math.hpp
  boost::random::mt19937 gen_;

 public:
  struct Particle{
    Eigen::VectorXd rea_vars; //!< Real variables
    Case *case_pointer; //!< Pointer to case

    //!< Velocity of the real variables
    Eigen::VectorXd rea_vars_velocity;

    Particle(Optimization::Case *c,
             boost::random::mt19937 &gen,
             Eigen::VectorXd v_max, int n_vars);

    Particle(){}

    void ParticleAdapt(Eigen::VectorXd rea_vars_velocity_swap,
                       Eigen::VectorXd rea_vars);
    double ofv() { return case_pointer->objf_value(); }
  };

  /*!
   * @brief
   * Generates a random set of cases within given upper and
   * lower bounds. The function also generates an initial
   * velocity based on the vMax parameter given through the
   * .json file.
   * @return
   */
  Case *generateRandomCase();

  /*!
   * @brief Looks through the memory of the swarm
   * to find the best evaluated perturbation.
   * @param swarm
   * @param current_best_particle_global
   * @return
   */
  Particle get_global_best();

  /*!
   * @brief Updates the velocity based on learning_factor_1_ (c1),
   * learning_factor_2_ (c2), the best evaluated perturbation of
   * the swarm and the best evaluated perturbation of that particle.
   * @param swarm_memory
   * @return
   */
  vector<PSO::Particle> update_velocity();

  /*!
   * @brief Updates the position based on the updated
   * velocities of the particles in the swarm.
   * @return
   */
  vector<PSO::Particle> update_position();

  /*!
   * @brief Prints the swarm and its current values
   * in a readable format, calls print particle
   * @param swarm
   */
  void printSwarm(vector<Particle> swarm = vector<Particle>()) const;

  /*!
   * @brief Prints the individual particle in a readable format
   * @param partic
   */

  void printParticle(Particle &partic) const;
  /*!
   * @brief Finds the best perturbation in the individual particle's memory
   * @param swarm_memory
   * @param particle_num
   * @return
   */
  Particle find_best_in_particle_memory(int particle_num);

  /*!
   * @brief Performs a check on the swarm, to figure out whether
   * it is stuck with particles that are too close to one another.
   * @return
   */
  bool is_stagnant();

  //!< Stagnation criterion:
  //!< standard deviation of all particle positions.
  double stagnation_limit_;

  //!< The memory of the swarm at previous time steps.
  vector<vector<Particle>> swarm_memory_;

  //!< # of particles in the swarm
  int number_of_particles_;

  double learning_factor_1_; //!< Learning factor 1 (c1)
  double learning_factor_2_; //!< Learning factor 2 (c2)

  //!< Max velocity of particle
  Eigen::VectorXd v_max_;

  int max_iterations_; //!< Max iterations

  //!< Current swarm of particles
  vector<Particle> swarm_;

  //!< Global best particle position
  Particle current_best_particle_global_;

  //!< Lower bounds for the variables (used for
  //!< randomly generating populations and mutation)
  Eigen::VectorXd lower_bound_;

  //!< Upper bounds for the variables (used for
  //!< randomly generating populations and mutation)
  Eigen::VectorXd upper_bound_;

  int n_vars_; //!< Number of variables in the problem.

  string im_ = "", wm_ = "", em_ = "";
  string cl_ = "PSO";
  string md_ = "Optimization::Optimizers";

};
}
}

#endif //FIELDOPT_PSO_H
