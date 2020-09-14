//
// Created by lukas on 11.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_MCMC_UPDATE_BASE_HPP
#define LATTICEMODELIMPLEMENTATIONS_MCMC_UPDATE_BASE_HPP


#include "param_helper/params.hpp"
#include "mcmc_simulation/util/random.hpp"


class MCMCUpdateBaseParameters : public Parameters {
public:
    using Parameters::Parameters;

    void write_to_file(const std::string& root_dir) {
        Parameters::write_to_file(root_dir, param_file_name());
    }

    static const std::string param_file_name()
    {
        return "mcmc_update_params";
    }
};

template <typename MCMCUpdate>
class MCMCUpdateBase {
public:
    template<typename Site>
    void initialize (const Site &site)
    {
        return mcmc_update_update().initialize_mcmc_update(site);
    }

    template<typename Site>
    void initialize_mcmc_update(const Site &site)
    {}

private:
    MCMCUpdate& mcmc_update_update() {
        return *static_cast<MCMCUpdate*>(this);
    }

    const MCMCUpdate& mcmc_update_update() const {
        return *static_cast<const MCMCUpdate*>(this);
    }
};

#endif //LATTICEMODELIMPLEMENTATIONS_MCMC_UPDATE_BASE_HPP
