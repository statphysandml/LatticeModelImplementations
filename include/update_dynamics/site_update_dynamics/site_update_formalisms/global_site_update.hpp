//
// Created by lukas on 02.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_GLOBAL_SITE_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_GLOBAL_SITE_UPDATE_HPP

#include "../site_update.hpp"


struct GlobalSiteUpdate;

struct GlobalSiteUpdateParameters : SiteUpdateParameters
{
    explicit GlobalSiteUpdateParameters(const json params_) : SiteUpdateParameters(params_)
    {}

    explicit GlobalSiteUpdateParameters() : GlobalSiteUpdateParameters(json {})
    {}

    static std::string name() {
        return "GlobalSiteUpdate";
    }

    typedef GlobalSiteUpdate SiteUpdate;
};

struct GlobalSiteUpdate : public SiteUpdate<GlobalSiteUpdate>
{
    explicit GlobalSiteUpdate(const GlobalSiteUpdateParameters &lp_) : lp(lp_)
    {}

    template<typename Site>
    void initialize_update(const Site& site)
    {}

    template<typename Site>
    void update (Site& site, uint measure_interval=1)
    {
        for(auto k = 0; k < measure_interval; k++)
        {
            global_lattice_update(site.get_update_formalism(), site);
        }
    }

    const GlobalSiteUpdateParameters &lp;
};

#endif //LATTICEMODELIMPLEMENTATIONS_GLOBAL_SITE_UPDATE_HPP
