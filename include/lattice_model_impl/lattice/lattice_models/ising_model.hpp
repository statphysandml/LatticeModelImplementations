//
// Created by lukas on 17.10.19.
//

#ifndef MAIN_ISING_MODEL_HPP
#define MAIN_ISING_MODEL_HPP

#include "../lattice_model.hpp"
#include "mcmc_simulation/util/random.hpp"
#include "param_helper/json.hpp"

namespace lm_impl
{
    namespace lattice_system
    {

        class IsingModel;

        class IsingModelParameters : public LatticeModelParameters {
        public:
            explicit IsingModelParameters(const json params_) : LatticeModelParameters(params_),
                                                                beta(get_entry<double>("beta", 0.4)),
                                                                J(get_entry<double>("J", 1.0)),
                                                                h(get_entry<double>("h", 0.0))
            {

            }

            explicit IsingModelParameters(double beta_, double J_, double h_) : IsingModelParameters(json{
                    {"beta", beta_},
                    {"J", J_},
                    {"h", h_}
            })
            {}

            const static std::string name() {
                return "IsingModel";
            }

            typedef IsingModel Model;

        private:
            friend class IsingModel;

            const double beta;
            const double J;
            const double h;
        };


        class IsingModel : public LatticeModel< IsingModel >
        {
        public:
            explicit IsingModel(const IsingModelParameters &mp_) : mp(mp_) {}

            template<typename T, typename T2=double_t>
            T2 get_potential(const T site, const std::vector<T*> neighbours)
            {
                double coupling = 0;
                for(size_t i = 0; i < neighbours.size(); i++) {
                    coupling += *neighbours[i];
                }
                return  -1.0 * mp.beta * site * (mp.J * coupling + mp.h);
            }

            template<typename T, typename T2=double_t>
            T2 get_energy_per_lattice_elem(const T site, const std::vector<T*> neighbours)
            {
                double coupling = 0;
                // Only neighbours in positive direction
                for(size_t i = 0; i < neighbours.size(); i += 2) {
                    coupling += *neighbours[i];
                }
                return  -1.0 * mp.beta * site * (mp.J * coupling + mp.h);
            }

        private:
            const IsingModelParameters &mp;
        };

        struct IsingModelSampler //  : public Sampler
        {
            IsingModelSampler(const double eps_)
            {
                uniint = std::uniform_int_distribution<int>(0, 1);
            }

            template<typename T>
            T random_state()
            {
                return 2 * uniint(mcmc::util::gen) - 1;
            }

            template<typename T>
            T propose_state(T site)
            {
                return -1 * site;
            }

            double get_eps() const
            {
                return 0.0;
            }

            const static std::string name() {
                return "IsingModelSampler";
            }

            std::uniform_int_distribution<int> uniint;
        };

    }
}



#endif //MAIN_ISING_MODEL_HPP
