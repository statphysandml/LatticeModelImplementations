//
// Created by lukas on 11.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_DUMMY_MCMC_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_DUMMY_MCMC_UPDATE_HPP


#include "mcmc_update_base.hpp"


namespace lm_impl {
    namespace mcmc_update {

        class DummyMCMCUpdate;


        class DummyMCMCUpdateParameters : public MCMCUpdateBaseParameters {
        public:
            explicit DummyMCMCUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_) {}

            explicit DummyMCMCUpdateParameters() : DummyMCMCUpdateParameters(json{}) {}

            static std::string name() {
                return "DummyMCMCUpdate";
            }

            typedef DummyMCMCUpdate MCMCUpdate;

        protected:
            friend class DummyMCMCUpdate;
        };


        class DummyMCMCUpdate : public MCMCUpdateBase<DummyMCMCUpdate, mcmc::sampler::DummySampler> {
        public:
            template<typename Model>
            explicit DummyMCMCUpdate(const DummyMCMCUpdateParameters &up_, Model &model_)
                    : MCMCUpdateBase<DummyMCMCUpdate, mcmc::sampler::DummySampler>(0.0), up(up_) {}

            template<typename T>
            T operator()(const T site) {
                return site;
            }

            template<typename T>
            T operator()(const T site, const std::vector<T *> neighbours) {
                return site;
            }

        protected:
            const DummyMCMCUpdateParameters &up;
        };

    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_DUMMY_MCMC_UPDATE_HPP
