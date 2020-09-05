//
// Created by lukas on 09.01.20.
//

#ifndef MAIN_LOCAL_UPDATE_FORMALISM_HPP
#define MAIN_LOCAL_UPDATE_FORMALISM_HPP

#include "param_helper/params.hpp"
#include "mcmc_simulation/util/random.hpp"


class LocalUpdateFormalismParameters : public Parameters {
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

template <typename LocalUpdate>
class LocalUpdateFormalism {
public:
    template<typename Lattice>
    void initialize (const Lattice &lattice)
    {
        return lattice_update().initialize_local_update(lattice);
    }

private:
    LocalUpdate& lattice_update() {
        return *static_cast<LocalUpdate*>(this);
    }

    const LocalUpdate& lattice_update() const {
        return *static_cast<const LocalUpdate*>(this);
    }
};

#endif //MAIN_LOCAL_UPDATE_FORMALISM_HPP
