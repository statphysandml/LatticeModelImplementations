//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP

#include "../site_update.hpp"


struct SiteSimpleUpdate;

struct SiteSimpleUpdateParameters : SiteUpdateParameters
{
    explicit SiteSimpleUpdateParameters(const json params_) : SiteUpdateParameters(params_)
    {}

    static std::string name() {
        return "SiteSimpleUpdate";
    }

    typedef SiteSimpleUpdate SiteUpdate;
};

struct SiteSimpleUpdate : public SiteUpdate<SiteSimpleUpdate>
{
    explicit SiteSimpleUpdate(const SiteSimpleUpdateParameters &lp_) : lp(lp_)
    {}

    template<typename Site>
    void initialize_update(const Site& site)
    {}

    template<typename Site>
    void update (Site& site, uint measure_interval=1)
    {
        for(auto k = 0; k < measure_interval; k++)
        {
            site.get_site() = update_lattice_site(site.get_update_formalism(), site.get_site());
        }
    }

    const SiteSimpleUpdateParameters &lp;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SIMPLE_UDATE_HPP
