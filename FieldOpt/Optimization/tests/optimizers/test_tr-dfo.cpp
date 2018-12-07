//
// Created by bellout on 11/27/18.
//

#include <gtest/gtest.h>

// Define settings for tr-dfo
#include "Optimization/tests/test_resource_optimizer.h"
#include "Optimization/tests/test_resource_cases.h"

#include "Runner/tests/test_resource_runner.hpp"
#include "Reservoir/tests/test_resource_grids.h"
#include "Optimization/tests/test_resource_test_functions.h"

#include "Optimization/optimizers/trust_region/TrustRegionOptimization.h"

#include <Eigen/Core>
#include "Utilities/math.hpp"
#include "Utilities/colors.hpp"
#include "Utilities/stringhelpers.hpp"

#include "test_tr-model-data.hpp"

using namespace TestResources::TestFunctions;
using namespace Optimization::Optimizers;
using namespace std;

namespace {

    class TrustRegionTest : public ::testing::Test,
                            public TestResources::TestResourceOptimizer,
                            public TestResources::TestResourceGrids {

    protected:
        TrustRegionTest() {}

        virtual ~TrustRegionTest() {}
        virtual void SetUp() {}

        TrustRegionOptimization *tr_dfo_;
        Optimization::Case *test_case_tr_dfo_probs_;
        VariablePropertyContainer *varcont_tr_dfo_probs_;
        TestResources::TrustRegionModelData tr_mdata;

        void SetUpOptimizer(TestResources::TrustRegionModelData::prob prob,
                            double (*tr_dfo_prob)(VectorXd xs)){

            VectorXd x0 = prob.xm.col(0);

            // Dummy var container based on initial point
            varcont_tr_dfo_probs_ = new VariablePropertyContainer();
            QString base_varname = "BHP#PRODUCER#"; // dummy var name

            for (int i = 0; i < x0.rows(); ++i) {
                // Use initial point values to construct container
                auto *prop = new ContinousProperty(x0(i));
                prop->setName(base_varname + QString::number(i));
                varcont_tr_dfo_probs_->AddVariable(prop);
            }

            FindVarSequence(prob);

            // Set up base case using dummy var containter
            test_case_tr_dfo_probs_ = new Optimization::Case(
                    QHash<QUuid, bool>(), QHash<QUuid, int>(),
                    varcont_tr_dfo_probs_->GetContinousVariableValues());

            // Use initial point from Matlab data
            test_case_tr_dfo_probs_->set_objective_function_value(tr_dfo_prob(x0));

            // cout << "SetUpOptimizer" << endl;
            tr_dfo_ = new TrustRegionOptimization(
                    settings_tr_opt_max_,
                    test_case_tr_dfo_probs_,
                    varcont_tr_dfo_probs_,
                    grid_5spot_,
                    logger_);
        }

        void FindVarSequence(TestResources::TrustRegionModelData::prob prob) {

          auto varcont = varcont_tr_dfo_probs_;
//
//          for (int ii=0; ii < prob.xm)
//
//          varcont.
        }

        void RunnerSubs(TestResources::TrustRegionModelData::prob prob,
                        double (*tr_dfo_prob)(VectorXd xs)){

            stringstream ss; ss << "[          ] " << FMAGENTA;
            cout << ss.str() << "-------------------------------------------------------" << END << endl;

            // Tentative best case should be equal to base case at this point
            CheckSameX(tr_dfo_->GetTentativeBestCase()->GetRealVarVector(),
                       prob.xm.col(0), 1e-6, "Check init point the same");

            while (!tr_dfo_->IsFinished()) {

                // RUNNER CALL (START)
                auto next_case = tr_dfo_->GetCaseForEvaluation();

                // COMPUTE OBJ.FUNCTION VALUE FOR CASE
                next_case->set_objective_function_value(
                        tr_dfo_prob(next_case->GetRealVarVector()));

                // PRINT CASE DATA (ID, X, F)
                cout << ss.str() << "-------------------------------------------------------" << END << endl;
                cout << ss.str() << "Case id:" << next_case->GetId().toString().toStdString() << END << endl;
                cout << ss.str() << "# initial points: " << tr_dfo_->GetNumInitPoints() << END << endl;
                cout << ss.str() << "x:" << next_case->GetRealVarVector().transpose() << END << endl;
                cout << ss.str() << "f:" << next_case->objective_function_value() << END << endl;


                // RUNNER CALL (FINISH)
                tr_dfo_->SubmitEvaluatedCase(next_case);

            }

            cout << ss.str() << "----------------------" << END << endl;
            auto best_case = tr_dfo_->GetTentativeBestCase();
            cout << best_case->objective_function_value() << endl;
        }

        void CheckSameX(VectorXd xa, VectorXd xb, double tol, std::string msg) {

          // cout << "xa: " << xa << endl;
          // cout << "xb: " << xb << endl;
          ASSERT_TRUE(xa.isApprox(xb,tol));
          stringstream ss; ss << "[          ] " << FMAGENTA;
          cout << ss.str() << msg.c_str() << " -> ok" << END << endl;
        }

    };

    //  auto values = variables->GetContinousVariableValues();
    //  VectorXd test_vec = VectorXd::Zero(values.size());
    //  Case * init_case = new Case(base_case);
    //  init_case->SetRealVarValues(test_vec);
    //  case_handler_->AddNewCase(init_case);

    TEST_F(TrustRegionTest, tr_dfo_prob1) {
        cout << FMAGENTA
             << "[CG.prob1  ] f = @(x) (1 - x(1))^2; x0=[-1.2 2.0]"
             << END << endl;
        SetUpOptimizer(tr_mdata.prob1, tr_dfo_prob1);
        RunnerSubs(tr_mdata.prob1, tr_dfo_prob1);
    }

    TEST_F(TrustRegionTest, tr_dfo_prob2) {
        cout << FMAGENTA
             << "[CG.prob2  ] f = @(x) log1p(x(1)^2) + x(2)^2; x0=[2.0 2.0]"
             << END << endl;
        SetUpOptimizer(tr_mdata.prob1, tr_dfo_prob2);
        RunnerSubs(tr_mdata.prob1, tr_dfo_prob2);
    }

    TEST_F(TrustRegionTest, tr_dfo_prob3) {
        cout << FMAGENTA
             << "[CG.prob3  ] f = @(x) sin(pi*x(1)/12) * cos(pi*x(2)/16); x0=[0.0 0.0]"
             << END << endl;
        SetUpOptimizer(tr_mdata.prob1, tr_dfo_prob3);
        RunnerSubs(tr_mdata.prob1, tr_dfo_prob3);
    }

    TEST_F(TrustRegionTest, tr_dfo_initial_model_created) {
        bool is_model_present = false;
        SetUpOptimizer(tr_mdata.prob1, tr_dfo_prob1);

        // RUNNER CALL (START)
        auto next_case = tr_dfo_->GetCaseForEvaluation();

        next_case->set_objective_function_value(
                tr_dfo_prob1(next_case->GetRealVarVector()));

        cout << "----------------------" << END << endl;
        cout << "Case id:" << next_case << END << endl;
        cout << "x:" << next_case->GetRealVarVector().transpose() << END << endl;
        cout << "f:" << next_case->objective_function_value() << END << endl;

        // RUNNER CALL (FINISH)
        tr_dfo_->SubmitEvaluatedCase(next_case);

        auto tr_model = tr_dfo_->getTrustRegionModel();

        if (tr_model) {
            if (tr_model->getDimension() >=2) { //<!at least 2 points are needed in the
                is_model_present = true;
            } else {
                cout << "Not enough points in the model." << endl;
            }
        } else {
            cout << "model not created" << endl;
        }
        EXPECT_TRUE(is_model_present);
    }
}
