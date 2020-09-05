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
        Parameters::write_to_file(root_dir, param_file_name());
    }

    static const std::string param_file_name()
    {
        return "model_params";
    }
};


template <typename ModelCL>
class LatticeModel
{
public:
    template<typename T>
    T normalize(T state)
    {
        return state;
    }

    template<typename SB>
    std::vector< std::unique_ptr<common_measures::MeasurePolicy<SB>> > generate_model_measures(const json& measure_names)
    {
        return std::vector< std::unique_ptr<common_measures::MeasurePolicy<SB>> > {};
    }

private:
    ModelCL& lattice_model() {
        return *static_cast<ModelCL*>(this);
    }
};

#endif //MAIN_MODEL_HPP
