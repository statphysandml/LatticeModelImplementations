//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_GLOBAL_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_GLOBAL_UPDATE_HPP

#include "param_helper/params.hpp"


class GlobalUpdateFormalismParameters : public Parameters {
public:
    using Parameters::Parameters;

    void write_to_file(const std::string& root_dir) {
        Parameters::write_to_file(root_dir, "global_updateformalism_params");
    }
};

template <typename GlobalUpdate>
class GlobalUpdateFormalism {
public:
    template<typename Lattice>
    void initialize (const Lattice &lattice)
    {
        return lattice_update().initialize_local_update(lattice);
    }

private:
    GlobalUpdate& lattice_update() {
        return *static_cast<GlobalUpdate*>(this);
    }

    const GlobalUpdate& lattice_update() const {
        return *static_cast<const GlobalUpdate*>(this);
    }
};

#endif //LATTICEMODELIMPLEMENTATIONS_GLOBAL_UPDATE_HPP
