//
// Created by lukas on 05.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP

#include "../lattice_update.hpp"


struct ParallelUpdate;

struct ParallelUpdateParameters : LatticeUpdateParameters
{
    explicit ParallelUpdateParameters(const json params_) : LatticeUpdateParameters(params_)
    {}

    static std::string name() {
        return "ParallelUpdate";
    }

    typedef ParallelUpdate LatticeUpdate;
};

struct ParallelUpdate : public LatticeUpdate<ParallelUpdate>
{
    explicit ParallelUpdate(const ParallelUpdateParameters &lp_) : lp(lp_)
    {}

    template<typename Lattice>
    void initialize_update(const Lattice& lattice)
    {}

    template<typename Lattice>
    void update (Lattice& lattice, uint measure_interval=1)
    {
        // ToDo: Introduce boost!
        for(auto j = 0; j < measure_interval; j++)
        {
            std::vector<typename Lattice::SiteType> lattice_grid_new(lattice.get_size(), typename Lattice::SiteType(0));

            // #pragma omp parallel for
            for(uint i = 0; i < lattice.get_size(); i++)
            {
                lattice_grid_new[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i], lattice.neighbours_at(i));
            }

            // ToDo: Rewrite?
            auto& lattice_grid = lattice.get_lattice_grid();
            lattice_grid = lattice_grid_new;
        }
    }

    const ParallelUpdateParameters &lp;
};

#endif //LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_HPP
