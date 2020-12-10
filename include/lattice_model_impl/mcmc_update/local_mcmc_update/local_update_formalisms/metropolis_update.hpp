//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_METROPOLIS_UPDATE_HPP
#define MAIN_METROPOLIS_UPDATE_HPP


#include "../../mcmc_update_base.hpp"


namespace lm_impl {
    namespace mcmc_update {

        namespace detail {
            // T operator() (const T site)
            template<typename T, typename UpdateFormalismCl, typename = void>
            struct is_updateable : std::false_type {
            };

            template<typename T, typename UpdateFormalismCl>
            struct is_updateable<T, UpdateFormalismCl, std::void_t<decltype(std::declval<UpdateFormalismCl>().template operator()<T>(
                    std::declval<const T>())),
                    decltype(std::declval<UpdateFormalismCl>().template operator()<T>(std::declval<const T>(),
                                                                                      std::declval<const double>(),
                                                                                      std::declval<double>()))> >
                    : std::true_type {
            };


            template<typename T, typename UpdateFormalismCl, typename = void>
            struct is_estimate_drift_term : std::false_type {
            };
            template<typename T, typename UpdateFormalismCl>
            // T operator() (const T site)
            struct is_estimate_drift_term<T, UpdateFormalismCl, std::void_t<decltype(std::declval<UpdateFormalismCl>().template estimate_drift_term<T>(
                    std::declval<const T>())),
                    decltype(std::declval<UpdateFormalismCl>().template estimate_drift_term<T>(std::declval<const T>(),
                                                                                               std::declval<const double>(),
                                                                                               std::declval<double>()))> >
                    : std::true_type {
            };
        }

        template<typename ModelParameters, typename SamplerCl>
        class MetropolisUpdate;


        template<typename ModelParameters, typename SamplerCl>
        class MetropolisUpdateParameters : public MCMCUpdateBaseParameters {
        public:
            explicit MetropolisUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_) {}

            explicit MetropolisUpdateParameters(const double eps)
                    : MetropolisUpdateParameters(json{{"eps", eps}}) {}

            static std::string name() {
                return "MetropolisUpdate";
            }

            typedef MetropolisUpdate<ModelParameters, SamplerCl> MCMCUpdate;

        protected:
            friend class MetropolisUpdate<ModelParameters, SamplerCl>;
        };


        template<typename ModelParameters, typename SamplerCl>
        class MetropolisUpdate : public MCMCUpdateBase<MetropolisUpdate<ModelParameters, SamplerCl>, SamplerCl> {
        public:
            explicit MetropolisUpdate(const MetropolisUpdateParameters<ModelParameters, SamplerCl> &up_,
                                      typename ModelParameters::Model &model_) :
                    MCMCUpdateBase<MetropolisUpdate<ModelParameters, SamplerCl>, SamplerCl>(up_.eps), up(up_),
                    model(model_) {
                rand = std::uniform_real_distribution<double>(0, 1);
            }

            template<typename T>
            T operator()(const T site) {
                T proposed_site = this->sampler.propose_state(site);
                if (rand(mcmc::util::gen) <
                    std::min(1.0, std::exp(-1.0 * (model.get_potential(proposed_site) - model.get_potential(site)))))
                    return proposed_site;
                else
                    return site;
            }

            template<typename T>
            T operator()(const T site, const std::vector<T *> neighbours) {
                T proposed_site = this->sampler.propose_state(site);
                if (rand(mcmc::util::gen) < std::min(1.0, std::exp(-1.0 *
                                                                   (model.get_potential(proposed_site, neighbours) -
                                                                    model.get_potential(site, neighbours)))))
                    return proposed_site;
                else
                    return site;
            }

        protected:
            const MetropolisUpdateParameters<ModelParameters, SamplerCl> &up;
            typename ModelParameters::Model &model;

            std::uniform_real_distribution<double> rand;
        };

    }
}

#endif //MAIN_METROPOLIS_UPDATE_HPP
