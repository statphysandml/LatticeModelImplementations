//
// Created by lukas on 05.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP
#define LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP

#include "../lattice_update.hpp"


struct ParallelUpdateWithAdpativeStepsize;

struct ParallelUpdateWithAdpativeStepsizeParameters : LatticeUpdateParameters
{
    explicit ParallelUpdateWithAdpativeStepsizeParameters(const json params_) : LatticeUpdateParameters(params_)
    {
        thermalization_steps = get_value_by_key<int>("thermalization_steps", 2000);
    }

    static std::string name() {
        return "ParallelUpdateWithAdpativeStepsize";
    }

    typedef ParallelUpdateWithAdpativeStepsize LatticeUpdate;

    int thermalization_steps;
};

struct ParallelUpdateWithAdpativeStepsize : public LatticeUpdate<ParallelUpdateWithAdpativeStepsize>
{
    explicit ParallelUpdateWithAdpativeStepsize(const ParallelUpdateWithAdpativeStepsizeParameters &lp_) : lp(lp_)
    {}

    template<typename Lattice>
    void initialize_update(const Lattice& lattice)
    {
        thermalization_counter = 0;
        /* function_ptr_on_update_func = (void(
                Lattice&)) &ParallelUpdateWithAdpativeStepsize::thermalization_phase_with_adpative_stepsize; */
    }

    template<typename Lattice>
    void update (Lattice& lattice, uint measure_interval=1)
    {
        // (this->*function_ptr_on_update_func) (lattice);
        if(thermalization_counter > lp.thermalization_steps)
            parallel_update_with_adpative_stepsize(lattice, measure_interval);
        else
            thermalization_phase_with_adpative_stepsize(lattice, measure_interval);
    }

    template<typename Lattice>
    void thermalization_phase_with_adpative_stepsize(Lattice& lattice, uint measure_interval=1)
    {
        // static_assert(detail::is_updateable<T, typename UpdateFormalismParameters::UpdateFormalism>::value, "is not estimate_drift_term");

        // ToDo Adapt measure interval
        if(measure_interval > 1)
            std::cout << "Thermalization with adaptive stepsize needs to be completed" << std::endl;

        double KMax = 0;
        for(uint i = 0; i < lattice.get_size(); i++)
        {
            const double K = std::fabs(lattice.get_update_formalism().estimate_drift_term(lattice[i], lattice.neighbours_at(i)));
            if(K > KMax)
                KMax = K;
        }

        if(thermalization_counter < lp.thermalization_steps) {
            KExpectation += KMax;
            thermalization_counter++;
        }
        else if(thermalization_counter == lp.thermalization_steps)
        {
            KExpectation /= lp.thermalization_steps;
            thermalization_counter++;
            // thermalization_counter = 0;
            /* function_ptr_on_update_func = (void(
                    Lattice)) &ParallelUpdateWithAdpativeStepsize::parallel_update_with_adpative_stepsize; */
        }

        std::vector<typename Lattice::SiteType> lattice_grid_new(lattice.get_size(), typename Lattice::SiteType(0));

        // #pragma omp parallel for
        for(uint i = 0; i < lattice.get_size(); i++)
        {
            lattice_grid_new[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i], lattice.neighbours_at(i));
        }
    }


    template<typename Lattice>
    void parallel_update_with_adpative_stepsize(Lattice& lattice, uint measure_interval=1)
    {
        // ToDo Adapt measure interval
        if(measure_interval > 1)
            std::cout << "Parallel update with adaptive stepsize needs to be completed" << std::endl;

        double KMax = 0;
        for(uint i = 0; i < lattice.get_size(); i++)
        {
            const double K = std::fabs(lattice.get_update_formalism().estimate_drift_term(lattice[i], lattice.neighbours_at(i)));
            if(K > KMax)
                KMax = K;
        }

        std::vector<typename Lattice::SiteType> lattice_grid_new(lattice.get_size(), typename Lattice::SiteType(0));

        // #pragma omp parallel for
        for(uint i = 0; i < lattice.get_size(); i++)
        {
            lattice_grid_new[i] = update_lattice_site(lattice.get_update_formalism(), lattice[i], lattice.neighbours_at(i), KMax, KExpectation);
        }

        // ToDo: Rewrite?
        auto& lattice_grid = lattice.get_lattice_grid();
        lattice_grid = lattice_grid_new;
    }

    const ParallelUpdateWithAdpativeStepsizeParameters &lp;

    // For adaptive step size
    double KExpectation = 0;
    int thermalization_counter = 0;

    // void (ParallelUpdateWithAdpativeStepsize::*function_ptr_on_update_func) (void);
};

#endif //LATTICEMODELIMPLEMENTATIONS_PARALLEL_UPDATE_WIH_ADAPATIVE_STEPSIZE_HPP
