/***********************************************************
 Created by Amanda Souza Machado on 12.02.20
 Copyright (C) 2020-2021 Amanda Machado
 <amanda.automacaoufsc@gmail.com>

 Modified 2019-2020 Thiago Lima Silva
 <thiagolims@gmail.com>

 Modified 2019-2020 Mathias Bellout
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

#include <gtest/gtest.h>
#include <iostream>
#include <fstream>

// Define settings for tr-dfo_exp-value
#include "Optimization/tests/test_resource_optimizer.h"
#include "Optimization/tests/test_resource_cases.h"

#include "Runner/tests/test_resource_runner.hpp"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"

#include "Optimization/optimizers/trust_region/TrustRegionOptimization.h"
#include "Optimization/optimizers/ensemble_exp_value.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include "Utilities/math.hpp"
#include "Utilities/colors.hpp"
#include "Utilities/stringhelpers.hpp"

#include "test_tr-model-data.hpp"
// #include "test_tr-support.hpp"


namespace fs = boost::filesystem;

//TODO make header file for test_tr-support.hpp
namespace TestResources {
void FindVarSequence(TestResources::TrustRegionModelData::prob &prob,
                     Optimization::Case &base_case_tr_dfo_probs);

void OverrideSecondPointEn(TestResources::TrustRegionModelData::prob &prob,
                           Optimization::Case &c,
                           Optimization::Optimizers::EnsembleExpValue *tr_en_);
}

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;
using namespace std;

namespace {

class EnTrTest : public ::testing::Test,
                 public TestResources::TestResourceOptimizer,
                 public TestResources::TestResourceGrids {

 protected:

  EnTrTest() {}

  virtual ~EnTrTest() {}
  virtual void SetUp() {}

  TrustRegionOptimization *tr_dfo_;
  EnsembleExpValue *tr_en_;
  Optimization::Case *test_case_tr_dfo_probs_;
  VarPropContainer *varcont_tr_dfo_probs_;
  TestResources::TrustRegionModelData tr_mdata;


  void SetUpOptimizer(TestResources::TrustRegionModelData::prob &prob){

    VectorXd x0 = prob.xm.col(0);

    // Dummy var container based on initial point
    varcont_tr_dfo_probs_ = new VarPropContainer();
    QString base_varname = "BHP#PRODUCER#"; // dummy var name

    for (int i = 0; i < x0.rows(); ++i) {
      // Use initial point values to construct container
      auto *prop = new ContinuousProperty(x0(i));
      prop->setName(base_varname + QString::number(i));
      varcont_tr_dfo_probs_->AddVariable(prop);
    }

    // Set up base case using dummy var containter
    test_case_tr_dfo_probs_ = new Optimization::Case(
        QList<QPair<QUuid, bool>>(),
        QList<QPair<QUuid, int>>(),
        varcont_tr_dfo_probs_->GetContVarValues());

    TestResources::FindVarSequence(prob,
                                   *test_case_tr_dfo_probs_);

    VectorXd ordered_vec(test_case_tr_dfo_probs_->GetRealVarVector().size());
    auto vec =test_case_tr_dfo_probs_->GetRealVarVector();
    for (int i=0; i <prob.idx.size(); i++) {
      ordered_vec(i) = vec(prob.idx[i]);
    }
    test_case_tr_dfo_probs_->SetRealVarValues(ordered_vec);

    // Use initial point from Matlab data
    test_case_tr_dfo_probs_->set_objf_value(tr_en_->initialValue(x0));

    tr_dfo_ = new TrustRegionOptimization(
        settings_tr_opt_max_,
        test_case_tr_dfo_probs_,
        varcont_tr_dfo_probs_,
        grid_5spot_,
        logger_);
  }

  bool RunnerSubs(TestResources::TrustRegionModelData::prob prob,
                  double (*tr_dfo_prob)(VectorXd xs)){

    stringstream ss; ss << "[          ] " << FMAGENTA;
    double tol = 1e-06;
    int p_count = 0;
    //bool is_ensemble_run_ = true;
    while (tr_dfo_->IsFinished()
        == Optimization::Optimizer::TerminationCondition::NOT_FINISHED) {

      auto new_case = tr_dfo_->GetCaseForEvaluation();

      while (new_case == nullptr) {
        if (tr_dfo_->IsFinished()) {
          break;
        } else {
          new_case = tr_dfo_->GetCaseForEvaluation();
        }
      }


      if (tr_dfo_->IsFinished()) {
        break;
      }

      // Evaluate for one ensemble element
      tr_en_->Evaluate(new_case->GetRealVarVector());


      // Check if all ensembles were evaluated
      if (tr_en_->IsCaseDone()) {
        // Calculate ExpValue and submit the case
        new_case->set_objf_value(tr_en_->ExpValue());
        if (p_count == 1 && prob.xm.cols() > 1) {
          TestResources::OverrideSecondPointEn(prob, *new_case, tr_en_);
        }
        tr_dfo_->SubmitEvaluatedCase(new_case);
        p_count++;

      }


    }

    // Printing the points calculated by EnsembleExpValue --------------------------
    Print_Log();

    return true;
  }
  void TearDown() override {
    /*
        TrustRegionOptimization *tr_dfo_;
        Optimization::Case *test_case_tr_dfo_probs_;
        VarPropContainer *varcont_tr_dfo_probs_;
        TestResources::TrustRegionModelData tr_mdata;
     */
    delete(tr_dfo_);
    delete(test_case_tr_dfo_probs_);
    delete (varcont_tr_dfo_probs_);
  }
  void Print_Log(){
    std::vector<std::string> points = tr_en_ -> getPrint();
    fs::path path =  fs::current_path();
    string path_string = path.string();
    path_string = path_string.substr(0,path_string.length() - 30);
    string subdir;
    subdir = path_string + "tools/python_scripts/ensemble_postprocessing/data.txt";
    //std::cout << "Current path is (data): " << subdir << '\n';
    std::ofstream myfile (subdir);
    if(myfile.is_open()){
      for (int element = 1; element < points.size(); ++element){
        if ( !(element == 3) ){
          myfile << points[element] << "\n";
        }
      }
      myfile.close();
    }
    else {cout << "Not possible to open file.";}
  }

};
// TEST_F(TrustRegionTest, trHS1) {
//     cout << endl << FMAGENTA << "[          ] ========p======"
//          << "=========================================== " << endl
//          << "[ HS1 ] "
//          << " 100*pow(x(1) - pow(x(0), 2), 2) + pow(1 - x(0), 2)" << END << endl;

//     SetUpOptimizer(tr_mdata.prob_hs1, hs1);
//     EXPECT_TRUE(RunnerSubs(tr_mdata.prob_hs1, hs1));
// }

TEST_F(EnTrTest, EntrDfoProb1) {
  cout << endl << FMAGENTA << "[          ] =============="
       << "=========================================== " << endl
       << "[ CG.prob1 ] "
       << " f = @(x) (1 - x(1))^2; x0=[-1.2 2.0]" << END << endl;

  tr_en_ = new EnsembleExpValue();
  tr_en_-> addFunction(Rosenbrock);
  //tr_en_-> addFunction(Rastingi_varation1);
//    tr_en_-> addFunction(Rastingi_varation2);
//    tr_en_-> addFunction(Rastingi_varation3);
//    tr_en_-> addFunction(Rastingi_varation4);
//    tr_en_-> addFunction(Rastingi_varation5);
//    tr_en_-> addFunction(Rastingi_varation6);
//    tr_en_-> addFunction(Rastingi_varation7);
//    tr_en_-> addFunction(Rastingi_varation8);
//    tr_en_-> addFunction(Rastingi_varation9);
//    tr_en_-> addFunction(Rastingi_varation10);

  // Clear file
  fs::path path =  fs::current_path();
  string path_string = path.string();
  path_string = path_string.substr(0,path_string.length() - 30);
  string subdir;
  subdir = path_string + "tools/python_scripts/ensemble_postprocessing/data_iterations.txt";
  //std::cout << "Current path is (data_iterations):" << subdir << '\n';
  std::ofstream ofs;
  ofs.open(subdir, std::ofstream::out | std::ofstream::trunc);
  ofs.close();

  SetUpOptimizer(tr_mdata.prob1);
  //cout << "SetUp completed" << endl;
  EXPECT_TRUE(RunnerSubs(tr_mdata.prob1, tr_dfo_prob1));
}

}
