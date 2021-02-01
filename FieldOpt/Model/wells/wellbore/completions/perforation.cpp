#include "perforation.h"

namespace Model {
namespace Wells {
namespace Wellbore {
namespace Completions {

Perforation::Perforation(Settings::Model::Well::Completion completion_settings,
                         Properties::VarPropContainer *variable_container)
    : Completion(completion_settings)
{
    transmissibility_factor_ = new Properties::ContinuousProperty(completion_settings.transmissibility_factor);
    if (completion_settings.is_variable) {
        transmissibility_factor_->setName(completion_settings.name);
        variable_container->AddVariable(transmissibility_factor_);
    }
}

Perforation::Perforation()
    : Completion(CompletionType::Perforation)
{
    transmissibility_factor_ = new Properties::ContinuousProperty(0.0);
}
}
}
}
}
