cat >"${include_path}/simulation_header.hpp" <<EOL
#ifndef MAIN_SIMULATION_HEADER_HPP
#define MAIN_SIMULATION_HEADER_HPP

#include <Python.h>

#include "mcmc_simulation/util/random.hpp"

#include "lattice_model_impl/site/site_header.hpp"
#include "lattice_model_impl/lattice/lattice_header.hpp"
#include "lattice_model_impl/link_lattice/link_lattice_header.hpp"

#include "lattice_model_impl/update_dynamics/update_dynamics_header.hpp"
#include "lattice_model_impl/mcmc_update/mcmc_update_header.hpp"

#include "execution/executer.hpp"

namespace from_file_simulation {
    /* typedef double BasicType;
    typedef XYModelParameters ModelParams;
    typedef HybridMonteCarloUpdateParameters<BasicType, ModelParams, GaussianSampler> UpdateParams;
    typedef LatticeParameters< BasicType, ModelParams, UpdateParams, GlobalLatticeUpdateParameters> SystemBaseParams; */

    template<typename ModelParameters>
    struct ComplexLangevinDynamics
    {
        typedef std::complex<double> BasicType;
        typedef ComplexLangevinUpdateParameters<ModelParameters, GaussianSampler> MCMCUpdateParams;
        typedef MemorySiteSimpleUpdateParameters<BasicType> UpdateDynamicsParams;
        typedef SiteParameters< BasicType, ModelParameters, MCMCUpdateParams, UpdateDynamicsParams > SystemBaseParams;
    };

    template<typename ModelParameters>
    struct MetropolisDynamics
    {
        typedef double BasicType;
        typedef MetropolisUpdateParameters<ModelParameters, GaussianSampler> MCMCUpdateParams;
        typedef SequentialUpdateParameters UpdateDynamicsParams;
        typedef LatticeParameters< BasicType, ModelParameters, MCMCUpdateParams, UpdateDynamicsParams> SystemBaseParams;

    };

    template<typename ModelParameters>
    void run_based_on_algorithm(int argc, char **argv)
    {
        std::string algorithm = std::string(argv[5]);

        if(algorithm == "ComplexLangevinDynamics")
            run_from_file<typename ComplexLangevinDynamics<ModelParameters>::SystemBaseParams>(argc, argv);
    }

    void run_based_on_model_and_algorithm(int argc, char **argv)
    {
        std::string model = std::string(argv[6]);

        if(model == "ComplexPolynomialModel")
            run_based_on_algorithm<ComplexPolynomialModelParameters>(argc, argv);
    }
}

#endif //MAIN_SIMULATION_HEADER_HPP

EOL