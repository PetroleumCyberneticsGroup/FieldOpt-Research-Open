/***********************************************************
Copyright (C) 2017
Einar J.M. Baumann <einar.baumann@gmail.com>
Created by einar on 11/23/18.

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

#include <gtest/gtest.h>
#include <Runner/tests/test_resource_runner.hpp>
#include "Optimization/optimizers/PSO.h"
#include "Optimization/tests/test_resource_optimizer.h"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"

#include "Model/tests/test_resource_model.h"

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;

namespace {

class PSOTest : public ::testing::Test,
                public TestResources::TestResourceOptimizer,
                public TestResources::TestResourceGrids
{
 protected:
  PSOTest() {
    base_ = base_case_;
    vars_ = model_rstrt_->variables();
    base_case_rstrt_ = new Optimization::Case(vars_->GetBinVarValues(),
                                              vars_->GetDiscVarValues(),
                                              vars_->GetContVarValues());

    base_case_rstrt_->set_objf_value(0.0);
    base_case_rstrt_->StringRepresentation(vars_);

    optmzr_ = new PSO(settings_rstrt_opt_,
                      base_case_rstrt_,
                      vars_,
                      grid_olympr37_,
                      logger_);
  }
  virtual ~PSOTest() {}
  virtual void SetUp() {}

  Optimization::Case *base_;
  Optimization::Case* base_case_rstrt_;
  VarPropContainer* vars_;
  PSO *optmzr_;
};

TEST_F(PSOTest, Constructor) {
}

TEST_F(PSOTest, TestFunctionSpherical) {
//    test_case_2r_->set_objective_function_value(abs(Sphere(test_case_2r_->GetRealVarVector())));
//    test_case_ga_spherical_30r_->set_objective_function_value(abs(Sphere(test_case_ga_spherical_30r_->GetRealVarVector())));
  test_case_ga_spherical_6r_->set_objf_value(abs(Sphere(test_case_ga_spherical_6r_->GetRealVarVector())));
  Optimization::Optimizer *minimizer = new PSO(settings_pso_min_,
                                               test_case_ga_spherical_6r_,
                                               varcont_6r_,
                                               grid_5spot_,
                                               logger_);

  while (!minimizer->IsFinished()) {
    auto next_case = minimizer->GetCaseForEvaluation();
    next_case->set_objf_value(abs(Sphere(next_case->GetRealVarVector())));
    minimizer->SubmitEvaluatedCase(next_case);
  }
  auto best_case = minimizer->GetTentBestCase();
  EXPECT_NEAR(0.0, best_case->objf_value(), 0.12);
  EXPECT_NEAR(0.0, best_case->GetRealVarVector()[0], 0.12);
  EXPECT_NEAR(0.0, best_case->GetRealVarVector()[1], 0.12);
}

TEST_F(PSOTest, TestFunctionRosenbrock) {
  test_case_ga_spherical_6r_->set_objf_value(abs(Rosenbrock(test_case_ga_spherical_6r_->GetRealVarVector())));
  settings_pso_min_->SetRngSeed(5);
  Optimization::Optimizer *minimizer = new PSO(settings_pso_min_,
                                               test_case_ga_spherical_6r_,
                                               varcont_6r_,
                                               grid_5spot_,
                                               logger_);

  while (!minimizer->IsFinished()) {
    auto next_case = minimizer->GetCaseForEvaluation();
    next_case->set_objf_value(Rosenbrock(next_case->GetRealVarVector()));
    minimizer->SubmitEvaluatedCase(next_case);
  }
  auto best_case = minimizer->GetTentBestCase();

  EXPECT_NEAR(0.0, best_case->objf_value(), 1);
  EXPECT_NEAR(1.0, best_case->GetRealVarVector()[0], 0.5);
  EXPECT_NEAR(1.0, best_case->GetRealVarVector()[1], 0.5);
}

TEST_F(PSOTest, TestRestart) {
  auto start = std::chrono::steady_clock::now();
  // std::cout << "[PSOTest]-t1 (ms)=" << Printer::since(start).count() << std::endl;
  // std::cout << "[PSOTest]-t2 (ms)=" << Printer::since(start).count() << std::endl;

  // testing against base_case value now returns error since base_case
  // is now loaded with restart value in the abstract_runner and not in
  // the optimizer constructor (as originally)

  map<string, map<string, double>> vmap;
  vmap.insert(
      {"bc", { // f: 0.7571595999999999
          { "Var#BHP#INJD-15#0", 199.96961185494197 },
          { "Var#BHP#INJD-15#1094", 209.54137953870247 },
          { "Var#BHP#INJD-15#1641", 250.0531920881195 },
          { "Var#BHP#INJD-15#2188", 197.0584663055752 },
          { "Var#BHP#INJD-15#2735", 200.6145349062493 },
          { "Var#BHP#INJD-15#3282", 195.43382947558968 },
          { "Var#BHP#INJD-15#3829", 253.89998863950598 },
          { "Var#BHP#INJD-15#4376", 218.0791077686635 },
          { "Var#BHP#INJD-15#4923", 250.38103437442064 },
          { "Var#BHP#INJD-15#547", 244.80148597015165 },
          { "Var#BHP#INJD-16#0", 236.31054225077148 },
          { "Var#BHP#INJD-16#1094", 238.83416476794562 },
          { "Var#BHP#INJD-16#1641", 201.13418297135198 },
          { "Var#BHP#INJD-16#2188", 239.78391825977764 },
          { "Var#BHP#INJD-16#2735", 254.14934644255555 },
          { "Var#BHP#INJD-16#3282", 228.72663503476676 },
          { "Var#BHP#INJD-16#3829", 195.01219735780882 },
          { "Var#BHP#INJD-16#4376", 247.82345714823884 },
          { "Var#BHP#INJD-16#4923", 202.85386597435658 },
          { "Var#BHP#INJD-16#547", 253.05227067702344 },
          { "Var#BHP#PRODX2#0", 167.36580743563388 },
          { "Var#BHP#PRODX2#1094", 105.43188190973554 },
          { "Var#BHP#PRODX2#1641", 117.34227644118684 },
          { "Var#BHP#PRODX2#2188", 157.77387174276436 },
          { "Var#BHP#PRODX2#2735", 173.02961678205642 },
          { "Var#BHP#PRODX2#3282", 106.125474949253 },
          { "Var#BHP#PRODX2#3829", 109.99697898592774 },
          { "Var#BHP#PRODX2#4376", 105.4197999414107 },
          { "Var#BHP#PRODX2#4923", 106.73833844615724 },
          { "Var#BHP#PRODX2#547", 105.00074442714711 },
          { "Var#ICD#PRODX2#0", 0.00653218586329762 },
          { "Var#ICD#PRODX2#1", 0.0001577555527071811 },
          { "Var#ICD#PRODX2#2", 0.00012941330135794304 },
          { "Var#ICD#PRODX2#3", 3.9575624016865985e-06 },
          { "Var#ICD#PRODX2#4", 1.2178067233547564e-05 },
          { "Var#ICD#PRODX2#5", 6.558205015438912e-07 },
          { "Var#ICD#PRODX2#6", 5.045821649893485e-06 },
          { "Var#ICD#PRODX2#7", 1.127772782414261e-05 }, }
      });

  vmap.insert(
      {"0", { // f: ??
          {"Var#BHP#INJD-15#0", 196.24936506946796},
          {"Var#BHP#INJD-15#1094", 200.0833560398039},
          {"Var#BHP#INJD-15#1641", 251.13343459311912},
          {"Var#BHP#INJD-15#2188", 200.1893227376486},
          {"Var#BHP#INJD-15#2735", 205.0350523328153},
          {"Var#BHP#INJD-15#3282", 216.1881559085568},
          {"Var#BHP#INJD-15#3829", 251.3247638358726},
          {"Var#BHP#INJD-15#4376", 253.4133455343952},
          {"Var#BHP#INJD-15#4923", 253.93635171506784},
          {"Var#BHP#INJD-15#547", 237.49513533302155},
          {"Var#BHP#INJD-16#0", 223.87405556061668},
          {"Var#BHP#INJD-16#1094", 215.13974789943242},
          {"Var#BHP#INJD-16#1641", 220.16614213459175},
          {"Var#BHP#INJD-16#2188", 252.87298163667134},
          {"Var#BHP#INJD-16#2735", 239.30586929969076},
          {"Var#BHP#INJD-16#3282", 253.67998487183758},
          {"Var#BHP#INJD-16#3829", 202.232645963778},
          {"Var#BHP#INJD-16#4376", 195.5902180542252},
          {"Var#BHP#INJD-16#4923", 195.0026225553358},
          {"Var#BHP#INJD-16#547", 254.15899803601545},
          {"Var#BHP#PRODX2#0", 165.9690380333624},
          {"Var#BHP#PRODX2#1094", 172.9728946737652},
          {"Var#BHP#PRODX2#1641", 105.23387404971061},
          {"Var#BHP#PRODX2#2188", 174.99681128750234},
          {"Var#BHP#PRODX2#2735", 161.93552396712982},
          {"Var#BHP#PRODX2#3282", 108.70769930787112},
          {"Var#BHP#PRODX2#3829", 131.12548574013024},
          {"Var#BHP#PRODX2#4376", 105.25991728061871},
          {"Var#BHP#PRODX2#4923", 105.27057804066898},
          {"Var#BHP#PRODX2#547", 127.19350439499866},
          {"Var#ICD#PRODX2#0", 0.005854802425182972},
          {"Var#ICD#PRODX2#1", 0.0006240847592807594},
          {"Var#ICD#PRODX2#2", 0.0003596574746618313},
          {"Var#ICD#PRODX2#3", 0.00044827061606302586},
          {"Var#ICD#PRODX2#4", 7.571113932905264e-06},
          {"Var#ICD#PRODX2#5", 2.7966798940608044e-05},
          {"Var#ICD#PRODX2#6", 3.4027344981138116e-06},
          {"Var#ICD#PRODX2#7", 2.5090771620907976e-05},}
      });

  vmap.insert(
      {"1", {
          { "Var#BHP#INJD-15#0", 230.9636994129001 },
          { "Var#BHP#INJD-15#1094", 237.4936347787147 },
          { "Var#BHP#INJD-15#1641", 253.72703095907366 },
          { "Var#BHP#INJD-15#2188", 195.31526491773323 },
          { "Var#BHP#INJD-15#2735", 212.6019643368932 },
          { "Var#BHP#INJD-15#3282", 195.08710708881216 },
          { "Var#BHP#INJD-15#3829", 253.92041949591834 },
          { "Var#BHP#INJD-15#4376", 210.3501986048761 },
          { "Var#BHP#INJD-15#4923", 254.1057123117979 },
          { "Var#BHP#INJD-15#547", 251.25862967494257 },
          { "Var#BHP#INJD-16#0", 211.11666552314384 },
          { "Var#BHP#INJD-16#1094", 247.10244486807932 },
          { "Var#BHP#INJD-16#1641", 197.24159581503685 },
          { "Var#BHP#INJD-16#2188", 242.64303515654933 },
          { "Var#BHP#INJD-16#2735", 253.65398580847489 },
          { "Var#BHP#INJD-16#3282", 236.2701161358217 },
          { "Var#BHP#INJD-16#3829", 195.43058665627728 },
          { "Var#BHP#INJD-16#4376", 254.64390054074892 },
          { "Var#BHP#INJD-16#4923", 220.03608143927158 },
          { "Var#BHP#INJD-16#547", 253.06691720141805 },
          { "Var#BHP#PRODX2#0", 105.03982648029881 },
          { "Var#BHP#PRODX2#1094", 106.3463855116187 },
          { "Var#BHP#PRODX2#1641", 127.33193085580048 },
          { "Var#BHP#PRODX2#2188", 172.40214022026504 },
          { "Var#BHP#PRODX2#2735", 155.05196665420254 },
          { "Var#BHP#PRODX2#3282", 105.03074031334575 },
          { "Var#BHP#PRODX2#3829", 108.97319784322359 },
          { "Var#BHP#PRODX2#4376", 120.71746318439087 },
          { "Var#BHP#PRODX2#4923", 106.18899656294829 },
          { "Var#BHP#PRODX2#547", 105.00473009290992 },
          { "Var#ICD#PRODX2#0", 0.004468277517232465 },
          { "Var#ICD#PRODX2#1", 0.007847513587302162 },
          { "Var#ICD#PRODX2#2", 5.841697280341111e-05 },
          { "Var#ICD#PRODX2#3", 1.7071353802747095e-05 },
          { "Var#ICD#PRODX2#4", 1.2007103018616717e-05 },
          { "Var#ICD#PRODX2#5", 1.4281191273194868e-06 },
          { "Var#ICD#PRODX2#6", 4.405052875854597e-06 },
          { "Var#ICD#PRODX2#7", 1.5999284711187833e-05 }, }
      });

  vmap.insert(
      { "2", {
          { "Var#BHP#INJD-15#0", 209.1366588747733 },
          { "Var#BHP#INJD-15#1094", 228.5130962351769 },
          { "Var#BHP#INJD-15#1641", 248.88729334105304 },
          { "Var#BHP#INJD-15#2188", 207.92926643878437 },
          { "Var#BHP#INJD-15#2735", 196.9768069479697 },
          { "Var#BHP#INJD-15#3282", 204.29494344219037 },
          { "Var#BHP#INJD-15#3829", 254.3262191015214 },
          { "Var#BHP#INJD-15#4376", 234.12597143479493 },
          { "Var#BHP#INJD-15#4923", 244.41777160293265 },
          { "Var#BHP#INJD-15#547", 236.57934566633367 },
          { "Var#BHP#INJD-16#0", 218.73327063581303 },
          { "Var#BHP#INJD-16#1094", 239.8087162468726 },
          { "Var#BHP#INJD-16#1641", 206.43274514501172 },
          { "Var#BHP#INJD-16#2188", 237.27564284908743 },
          { "Var#BHP#INJD-16#2735", 253.29989218404359 },
          { "Var#BHP#INJD-16#3282", 195.2395880432801 },
          { "Var#BHP#INJD-16#3829", 211.43397419698985 },
          { "Var#BHP#INJD-16#4376", 254.10855103049045 },
          { "Var#BHP#INJD-16#4923", 219.55725648463678 },
          { "Var#BHP#INJD-16#547", 252.14408626262463 },
          { "Var#BHP#PRODX2#0", 167.20927017455207 },
          { "Var#BHP#PRODX2#1094", 133.46360645315963 },
          { "Var#BHP#PRODX2#1641", 125.59475785181948 },
          { "Var#BHP#PRODX2#2188", 174.22878998426043 },
          { "Var#BHP#PRODX2#2735", 172.51124975890733 },
          { "Var#BHP#PRODX2#3282", 109.90517650991568 },
          { "Var#BHP#PRODX2#3829", 110.94303111873506 },
          { "Var#BHP#PRODX2#4376", 107.79839972431944 },
          { "Var#BHP#PRODX2#4923", 107.89898148659594 },
          { "Var#BHP#PRODX2#547", 125.18465586858451 },
          { "Var#ICD#PRODX2#0", 0.005981304960689753 },
          { "Var#ICD#PRODX2#1", 0.0003006601202991568 },
          { "Var#ICD#PRODX2#2", 0.0005311153419789694 },
          { "Var#ICD#PRODX2#3", 3.024245665072308e-05 },
          { "Var#ICD#PRODX2#4", 8.187637950250693e-05 },
          { "Var#ICD#PRODX2#5", 1.4763030178836095e-06 },
          { "Var#ICD#PRODX2#6", 3.077253825368477e-06 },
          { "Var#ICD#PRODX2#7", 9.408154603518168e-06 }, }
      });

  vmap.insert(
      { "3", {
          { "Var#BHP#INJD-15#0", 201.26959535017778 },
          { "Var#BHP#INJD-15#1094", 228.42623259123948 },
          { "Var#BHP#INJD-15#1641", 250.86595138911682 },
          { "Var#BHP#INJD-15#2188", 195.00000762918356 },
          { "Var#BHP#INJD-15#2735", 197.1003480434137 },
          { "Var#BHP#INJD-15#3282", 254.5689857887454 },
          { "Var#BHP#INJD-15#3829", 253.48445866817974 },
          { "Var#BHP#INJD-15#4376", 212.26307208748955 },
          { "Var#BHP#INJD-15#4923", 195.21909719784713 },
          { "Var#BHP#INJD-15#547", 223.46583705975038 },
          { "Var#BHP#INJD-16#0", 254.5573001089755 },
          { "Var#BHP#INJD-16#1094", 236.89869457819546 },
          { "Var#BHP#INJD-16#1641", 197.9911762733692 },
          { "Var#BHP#INJD-16#2188", 198.86458569992243 },
          { "Var#BHP#INJD-16#2735", 243.41477992389605 },
          { "Var#BHP#INJD-16#3282", 206.86137275868924 },
          { "Var#BHP#INJD-16#3829", 195.12490460718016 },
          { "Var#BHP#INJD-16#4376", 227.32892560126248 },
          { "Var#BHP#INJD-16#4923", 195.87578460893056 },
          { "Var#BHP#INJD-16#547", 237.8074799300802 },
          { "Var#BHP#PRODX2#0", 122.12636519950598 },
          { "Var#BHP#PRODX2#1094", 115.06377366511958 },
          { "Var#BHP#PRODX2#1641", 155.4211326018974 },
          { "Var#BHP#PRODX2#2188", 171.78281554588634 },
          { "Var#BHP#PRODX2#2735", 170.79330043469733 },
          { "Var#BHP#PRODX2#3282", 113.54577353651757 },
          { "Var#BHP#PRODX2#3829", 108.56294785278985 },
          { "Var#BHP#PRODX2#4376", 110.03756136483824 },
          { "Var#BHP#PRODX2#4923", 108.85666867543799 },
          { "Var#BHP#PRODX2#547", 109.72614483224864 },
          { "Var#ICD#PRODX2#0", 0.007773448308697457 },
          { "Var#ICD#PRODX2#1", 6.9278098882574746e-12 },
          { "Var#ICD#PRODX2#2", 0.00012952460210536015 },
          { "Var#ICD#PRODX2#3", 3.4893675225370225e-05 },
          { "Var#ICD#PRODX2#4", 3.9549037222983644e-05 },
          { "Var#ICD#PRODX2#5", 1.2757495913367188e-05 },
          { "Var#ICD#PRODX2#6", 1.0417444207641569e-05 },
          { "Var#ICD#PRODX2#7", 3.7220256079980706e-05 }, }
      });

  // EXPECT_NEAR(0.7571595999999999, optmzr_->GetTentBestCase()->objf_value(), 1e-12); // bc
  // temp commented out until finding the fval of the particle:
  // EXPECT_NEAR(??????, optmzr_->getSwarmMem()[0][0].case_pointer->objf_value(), 1e-12); //

  // check base case
  for ( auto const &v : vmap.at("0")) {
    auto uuid = vars_->GetContVarQUuid(QString::fromStdString(v.first).replace("Var#",""));
    auto val = optmzr_->getSwarmMem()[0][0].case_pointer->get_real_variable_value(uuid);
    vars_->GetContinuousVariable(uuid)->scaleBackOtherValue(val);
    // auto val = optmzr_->GetTentBestCase()->get_real_variable_value(uuid);
    EXPECT_NEAR(val, v.second, 1e-12);
  }
  // std::cout << "[PSOTest]-t3 (ms)=" << Printer::since(start).count() << std::endl;

  // Figure out QHash numbering to select corresponding case,
  // i.e., which particle ended in which rank-dir since
  // QHash<QUuid, Case *> cases_; and
  // QQueue<QUuid> evaluation_queue_;
  // in case_handler

  // for (int ii=0; ii < vmap.size(); ii++) {
  //   auto c = optmzr->case_handler()->AllCases()[ii];
  //   for ( auto v : vmap.at(to_string(ii))) {
  //     auto uuid = vars->GetContVarQUuid(QString::fromStdString(v.first).replace("Var#",""));
  //     auto val = c->get_real_variable_value(uuid);
  //     EXPECT_NEAR(val, v.second, 1e-12);
  //   }
  // }
  // std::cout << "[PSOTest]-t4 (ms)=" << Printer::since(start).count() << std::endl;
}

}
