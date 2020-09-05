//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP

#include "../lattice_update.hpp"


struct GlobalLatticeUpdate;

struct GlobalLatticeUpdateParameters : LatticeUpdateParameters
{
    explicit GlobalLatticeUpdateParameters(const json params_) : LatticeUpdateParameters(params_)
    {}

    explicit GlobalLatticeUpdateParameters() : GlobalLatticeUpdateParameters(json {})
    {}

    static std::string name() {
        return "GlobalLatticeUpdate";
    }

    typedef GlobalLatticeUpdate LatticeUpdate;
};

struct GlobalLatticeUpdate : public LatticeUpdate<GlobalLatticeUpdate>
{
    explicit GlobalLatticeUpdate(const GlobalLatticeUpdateParameters &lp_) : lp(lp_)
    {}

    template<typename Lattice>
    void initialize_update(const Lattice& lattice)
    {}

    template<typename Lattice>
    void update (Lattice& lattice, uint measure_interval=1)
    {
        for(auto k = 0; k < measure_interval; k++)
        {
            global_lattice_update(lattice.get_update_formalism(), lattice);
        }
    }

    const GlobalLatticeUpdateParameters &lp;
};

#endif //LATTICEMODELIMPLEMENTATIONS_SIMPLE_UPDATE_HPP
