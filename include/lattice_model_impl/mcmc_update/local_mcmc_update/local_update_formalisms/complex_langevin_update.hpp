//
// Created by lukas on 09.10.19.
//

#ifndef MAIN_COMPLEX_LANGEVIN_UPDATE_HPP
#define MAIN_COMPLEX_LANGEVIN_UPDATE_HPP


#include "../../mcmc_update_base.hpp"


namespace lm_impl {
    namespace mcmc_update {

        template<typename ModelParameters, typename SamplerCl>
        class ComplexLangevinUpdate;


        template<typename ModelParameters, typename SamplerCl>
        class ComplexLangevinUpdateParameters : public MCMCUpdateBaseParameters {
        public:
            explicit ComplexLangevinUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                           epsilon(get_entry<double>("epsilon", eps)),
                                                                           sqrt2epsilon(sqrt(2 * get_entry<double>(
                                                                                   "epsilon", eps))) {}

            explicit ComplexLangevinUpdateParameters(
                    const double epsilon_
            ) : ComplexLangevinUpdateParameters(json{
                    {"epsilon", epsilon_},
                    {"eps",     epsilon_}
            }) {}

            static std::string name() {
                return "ComplexLangevinUpdate";
            }

            typedef ComplexLangevinUpdate<ModelParameters, SamplerCl> MCMCUpdate;

        private:
            friend class ComplexLangevinUpdate<ModelParameters, SamplerCl>;

            const double epsilon;
            const double sqrt2epsilon;
        };


        template<typename ModelParameters, typename SamplerCl>
        class ComplexLangevinUpdate
                : public MCMCUpdateBase<ComplexLangevinUpdate<ModelParameters, SamplerCl>, SamplerCl> {
        public:
            explicit ComplexLangevinUpdate(const ComplexLangevinUpdateParameters<ModelParameters, SamplerCl> &up_,
                                           typename ModelParameters::Model &model_)
                    : MCMCUpdateBase<ComplexLangevinUpdate<ModelParameters, SamplerCl>, SamplerCl>(up_.eps), up(up_),
                      model(model_) {
                normal = std::normal_distribution<double>(0, 1);
            }

            template<typename T>
            T estimate_drift_term(const T site) {
                return model.get_drift_term(site);
            }

            template<typename T>
            T estimate_drift_term(const T site, const std::vector<T *> neighbours) {
                return model.get_drift_term(site, neighbours);
            }

            template<typename T>
            T operator()(const T site) {
                return update(site, model.get_drift_term(site), up.epsilon, up.sqrt2epsilon);
            }

            template<typename T>
            T operator()(const T site, const std::vector<T *> neighbours) {
                return update(site, model.get_drift_term(site, neighbours), up.epsilon, up.sqrt2epsilon);
            }

            template<typename T>
            T operator()(const T site, const double KMax, const double KExpectation) {
                T eps_drift_term = model.get_drift_term(site);
                double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);
                return update(site, eps_drift_term, epsilon, std::sqrt(2 * epsilon));
            }

            template<typename T>
            T
            operator()(const T site, const std::vector<T *> neighbours, const double KMax, const double KExpectation) {
                T eps_drift_term = model.get_drift_term(site, neighbours);
                double epsilon = std::min(up.epsilon, up.epsilon * KExpectation / KMax);
                return update(site, eps_drift_term, epsilon, std::sqrt(2 * epsilon));
            }

        private:
            const ComplexLangevinUpdateParameters<ModelParameters, SamplerCl> &up;
            typename ModelParameters::Model &model;

            std::normal_distribution<double> normal;

            template<typename T>
            T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon) {
                T new_site = {site.real(), site.imag()};
                new_site.real(site.real() - epsilon * drift_term.real() + sqrt2epsilon * normal(mcmc::util::gen));
                new_site.imag(site.imag() - epsilon * drift_term.imag());

                return model.normalize(new_site);
            }
        };

    }
}

#endif //MAIN_COMPLEX_LANGEVIN_UPDATE_HPP
