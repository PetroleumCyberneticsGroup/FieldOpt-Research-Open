#ifndef RUNNER_H
#define RUNNER_H

#include <stdexcept>
#include "runtime_settings.h"
#include "abstract_runner.h"

namespace Runner {

/*!
 * \brief The MainRunner class takes care of initializing and starting the actual runner indicated in the runtime settings.
 */
    class MainRunner
    {
    public:
        MainRunner(RuntimeSettings *rts);
        MainRunner(RuntimeSettings *rts, Optimization::Case* base_case, Model::ModelSynchronizationObject* mso);

        /*!
         * \brief Execute Start the optimization run by calling the Execute function in the simulator.
         */
        void Execute();

        Model::Model* getModel() { return runner_->getModel(); }
        Settings::Settings* getSettings() { return runner_->getSettings(); }
        Optimization::Optimizer* getOptimizer() { return runner_->getOptimizer(); }
        Simulation::Simulator* getSimulator() { return runner_->getSimulator(); }

        /*!
         * @brief Initializes runner modules.
        */
        void InitializeModules() { runner_->InitializeModules(); }

        /*!
         * @brief Replace the optimizer.
         * @param opt New optimizer setttings.
        */
        void ReplaceOptimizer(Settings::Optimizer *opt) { runner_->ReplaceOptimizer(opt); }


        /*!
         * @brief Evaluates the base case.
         * @return Base case with objective value.
        */
        Optimization::Case* evaluateBaseCase();

    private:
        RuntimeSettings *runtime_settings_;
        AbstractRunner *runner_;
    };

}

#endif // RUNNER_H
