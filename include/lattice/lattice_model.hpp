//
// Created by lukas on 09.10.19.
//

#ifndef MAIN_MODEL_HPP
#define MAIN_MODEL_HPP


#include "param_helper/params.hpp"
#include "mcmc_simulation/measure_policy.hpp"


class LatticeModelParameters : public Parameters {
public:
    using Parameters::Parameters;

    const void write_to_file(const std::string& root_dir) {
        Parameters::write_to_file(root_dir, "model_params");
    }
};


template <typename ModelCL>
class LatticeModel
{
public:
    template<typename SB, typename SBP>
    MeasurePolicy<SB>* model_measure_factory(const std::string& measure, const SBP& system_base_parameters) {
        return nullptr;
    }

private:
    ModelCL& lattice_model() {
        return *static_cast<ModelCL*>(this);
    }
};

#endif //MAIN_MODEL_HPP
