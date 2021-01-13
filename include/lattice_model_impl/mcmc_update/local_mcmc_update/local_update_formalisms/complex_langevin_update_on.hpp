//
// Created by lukas on 12.01.21.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_UPDATE_ON_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_UPDATE_ON_HPP


#include "../../mcmc_update_base.hpp"


namespace lm_impl {
    namespace mcmc_update {

        struct ComplexONModelSampler;

        template<typename ModelParameters>
        class ComplexLangevinUpdateON;


        template<typename ModelParameters>
        class ComplexLangevinUpdateONParameters : public MCMCUpdateBaseParameters {
        public:
            explicit ComplexLangevinUpdateONParameters(const json params_) : MCMCUpdateBaseParameters(params_),
                                                                           epsilon(get_entry<double>("epsilon", eps)),
                                                                           sqrt2epsilon(sqrt(2 * get_entry<double>(
                                                                                   "epsilon", eps))) {}

            explicit ComplexLangevinUpdateONParameters(
                    const double epsilon_
            ) : ComplexLangevinUpdateONParameters(json{
                    {"epsilon", epsilon_},
                    {"eps",     epsilon_}
            }) {}

            static std::string name() {
                return "ComplexLangevinUpdateON";
            }

            typedef ComplexLangevinUpdateON<ModelParameters> MCMCUpdate;

        private:
            friend class ComplexLangevinUpdateON<ModelParameters>;

            const double epsilon;
            const double sqrt2epsilon;
        };


        template<typename ModelParameters>
        class ComplexLangevinUpdateON
                : public MCMCUpdateBase<ComplexLangevinUpdateON<ModelParameters>, ComplexONModelSampler> {
        public:
            explicit ComplexLangevinUpdateON(const ComplexLangevinUpdateONParameters<ModelParameters> &up_,
                                           typename ModelParameters::Model &model_)
                    : MCMCUpdateBase<ComplexLangevinUpdateON<ModelParameters>, ComplexONModelSampler>(up_.eps), up(up_),
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
            const ComplexLangevinUpdateONParameters<ModelParameters> &up;
            typename ModelParameters::Model &model;
            std::vector<double> epsilon;

            std::normal_distribution<double> normal;

            template<typename T>
            T update(const T site, const T drift_term, const double &epsilon, const double &sqrt2epsilon) {
                T new_site(0);
                for(uint i = 0; i < new_site.dim(); i++) {
                    new_site(i).real(site(i).real() - epsilon * drift_term(i).real() + sqrt2epsilon * normal(mcmc::util::gen));
                    new_site(i).imag(site(i).imag() - epsilon * drift_term(i).imag());
                }
                return model.normalize(new_site);
            }
        };

        struct ComplexONModelSampler
        {
            ComplexONModelSampler(const double eps_) : eps(eps_)
            {
                normal = std::normal_distribution<double>(0, 1);
            }

            template<typename T>
            T random_state() {
                T new_site(0);
                for(uint i = 0; i < new_site.dim(); i++)
                    new_site(i) += std::sqrt(2 * eps) * normal(mcmc::util::gen);
                return new_site;
            }

            template<typename T>
            T propose_state(T site) {
                return site + random_state<T>();
            }

            double get_eps() const
            {
                return eps;
            }

            const static std::string name() {
                return "ComplexONModelSampler";
            }

            const double eps;
            std::normal_distribution<double> normal;
        };
    }
}

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_LANGEVIN_UPDATE_ON_HPP
