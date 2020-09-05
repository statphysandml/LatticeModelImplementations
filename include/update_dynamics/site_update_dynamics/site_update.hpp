//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SITE_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SITE_UPDATE_HPP

#include "param_helper/params.hpp"


class SiteUpdateParameters : public Parameters {
public:
    using Parameters::Parameters;

    const void write_to_file(const std::string& root_dir) {
        Parameters::write_to_file(root_dir, param_file_name());
    }

    static const std::string param_file_name()
    {
        return "update_dynamics_params";
    }
};


template<typename Derived>
class SiteUpdate
{
public:
    template<typename Site>
    void initialize (const Site &site)
    {
        return site_update().initialize_update(site);
    }

    template<typename Site>
    void operator() (Site &site, uint measure_interval=1)
    {
        return site_update().update(site, measure_interval);
    }
private:
    Derived& site_update() {
        return *static_cast<Derived*>(this);
    }

    const Derived& site_update() const {
        return *static_cast<const Derived*>(this);
    }
};

#endif //LATTICEMODELIMPLEMENTATIONS_SITE_UPDATE_HPP
