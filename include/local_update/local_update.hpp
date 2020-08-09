//
// Created by lukas on 09.01.20.
//

#ifndef MAIN_LOCAL_UPDATE_FORMALISM_HPP
#define MAIN_LOCAL_UPDATE_FORMALISM_HPP

#include "param_helper/params.hpp"


class LocalUpdateFormalismParameters : public Parameters {
public:
    using Parameters::Parameters;

    void write_to_file(const std::string& root_dir) {
        Parameters::write_to_file(root_dir, "local_updateformalism_params");
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
