//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_SITE_HPP
#define MAIN_SITE_HPP

#include "mcmc_simulation/header.hpp"

#include "../util/measures/lattice_measures.hpp"
#include "../update_dynamics/site_update_dynamics/site_update_formalisms/simple_update.hpp"


namespace lm_impl {
    namespace site_system {

        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters=update_dynamics::SiteSimpleUpdateParameters>
        class SiteSystem;


        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters=update_dynamics::SiteSimpleUpdateParameters>
        class SiteParameters : public mcmc::simulation::SystemBaseParameters {
        public:
            explicit SiteParameters(const json params_) : SystemBaseParameters(params_) {
                model_parameters = std::make_unique<ModelParameters>(
                        mcmc::util::generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, ModelParameters>(
                                *this, model_parameters->param_file_name()));

                update_parameters = std::make_unique<UpdateFormalismParameters>(
                        mcmc::util::generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, UpdateFormalismParameters>(
                                *this, update_parameters->param_file_name()));

                site_update_parameters = std::make_unique<SiteUpdateFormalismParameters>(
                        mcmc::util::generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, SiteUpdateFormalismParameters>(
                                *this, site_update_parameters->param_file_name()));
            }

            // Move constructor
            SiteParameters(SiteParameters &&site_parameters) :
                    SystemBaseParameters({"measures", site_parameters.measures}),
                    model_parameters(std::move(site_parameters.model_parameters)),
                    update_parameters(std::move(site_parameters.update_parameters)),
                    site_update_parameters(std::move(
                            site_parameters.site_update_parameters)) { /* std::cout << "Move cosntructor is called" << std::endl;*/ }

            // Move assignment
            SiteParameters &operator=(
                    SiteParameters &site_parameters) // Changed on my own from no & to && (from DevDat other to &&other)
            {

                // std::cout << "Move assignment operator is called" << std::endl;
                using std::swap;
                swap(measures, site_parameters.measures);
                model_parameters = std::move(site_parameters.model_parameters);
                update_parameters = std::move(site_parameters.update_parameters);
                site_update_parameters = std::move(site_parameters.site_update_parameters);
                return *this;
            }

            typedef SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> System;

            std::unique_ptr<System> generate() {
                return std::make_unique<System>(*this);
            }

            void write_to_file(const std::string &rel_config_path) override {
                std::string model_params_path = get_entry<std::string>(model_parameters->param_file_name() + "_path",
                                                                       rel_config_path);
                model_parameters->write_to_file(model_params_path);

                std::string update_formalism_params_path = get_entry<std::string>(
                        update_parameters->param_file_name() + "_path", rel_config_path);
                update_parameters->write_to_file(update_formalism_params_path);

                std::string site_update_params_path = get_entry<std::string>(
                        site_update_parameters->param_file_name() + "_path", rel_config_path);
                site_update_parameters->write_to_file(site_update_params_path);

                json model_parameters_ = model_parameters->get_json();
                delete_entry(model_parameters->param_file_name());

                json update_parameters_ = update_parameters->get_json();
                delete_entry(update_parameters->param_file_name());

                json site_update_parameters_ = site_update_parameters->get_json();
                delete_entry(site_update_parameters->param_file_name());

                Parameters::write_to_file(rel_config_path, name());

                add_entry(model_parameters->param_file_name(), model_parameters_);
                add_entry(update_parameters->param_file_name(), update_parameters_);
                add_entry(site_update_parameters->param_file_name(), site_update_parameters_);
            }

            Parameters build_expanded_raw_parameters() const override {
                Parameters parameters(params);
                parameters.add_entry(model_parameters->param_file_name(), model_parameters->get_json());
                parameters.add_entry(update_parameters->param_file_name(), update_parameters->get_json());
                parameters.add_entry(site_update_parameters->param_file_name(), site_update_parameters->get_json());
                return parameters;
            }

            // Only used for testing in simulation.hpp
            typedef ModelParameters MP_;
            typedef UpdateFormalismParameters UP_;

        private:
            template<typename, typename, typename, typename>
            friend
            class SiteSystem;

            std::unique_ptr<ModelParameters> model_parameters;
            std::unique_ptr<UpdateFormalismParameters> update_parameters;
            std::unique_ptr<SiteUpdateFormalismParameters> site_update_parameters;
        };

        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
        class SiteSystem
                : public mcmc::simulation::SystemBase<SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > {
        public:
            explicit SiteSystem(
                    const SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &sp_)
                    : sp(sp_) {
                model = std::make_unique<typename ModelParameters::Model>(*sp.model_parameters);
                update_formalism = std::make_unique<typename UpdateFormalismParameters::MCMCUpdate>(
                        *sp.update_parameters, *model);
                site_update = std::make_unique<typename SiteUpdateFormalismParameters::UpdateDynamics>(
                        *sp.site_update_parameters);

                // Needs to stay here since other site_update or update_formalism use site as reference
                initialize_site();

                update_formalism->initialize(*this);
                site_update->initialize(*this);
            }

            SiteSystem(
                    std::unique_ptr<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > sp_ptr_)
                    : SiteSystem(*(sp_ptr_.release())) {}

            // Move constructor
            SiteSystem(SiteSystem &&site_system) :
            // sp_ptr(std::move(site_system.sp_ptr)),
            // meausures do also need to be moved, right? -> not implemented so far
                    sp(site_system.sp),
                    site(site_system.site),
                    model(std::move(site_system.model)),
                    update_formalism(std::move(site_system.update_formalism)),
                    site_update(std::move(
                            site_system.site_update)) {
                /* std::cout << "Move cosntructor is called" << std::endl; */ }

            // Move assignment
            SiteSystem &
            operator=(SiteSystem &site_system) // Changed on my own from no & to && (from DevDat other to &&other)
            {
                // std::cout << "Move assignment operator is called" << std::endl;
                using std::swap;
                swap(sp, site_system.sp);
                swap(site, site_system.site);
                model = std::move(site_system.model);
                update_formalism = std::move(site_system.update_formalism);
                site_update = std::move(site_system.site_update);
                return *this;
            }

            static SiteSystem from_json(const json params_) {
                // Object is generated on the heap
                auto sp_ptr_ = std::make_unique<
                        SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>>(params_);
                return SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>(
                        *(sp_ptr_.release()));
            }

            void update_step(uint measure_interval) {
                site_update->operator()(*this, measure_interval);
            }

            void initialize(std::string starting_mode) {
                this->generate_measures(sp.measures);

                if (starting_mode == "hot")
                    site = update_formalism->template random_state<T>();
                else
                    site = update_formalism->template cold_state<T>();
            }

            const uint get_size() const {
                return 1;
            }

            auto &at(int i) {
                return site;
            }

            auto at(int i) const {
                return site;
            }

            auto &get_system_representation() {
                return site;
            }

            auto get_system_representation() const {
                return site;
            }

            void generate_measures(const json &measure_names) override {
                auto lattice_related_measures = generate_site_system_measures(sp.measures);
                this->concat_measures(lattice_related_measures);

                auto model_related_measures = model->template generate_model_measures<SiteSystem,
                        SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>>(sp);
                this->concat_measures(model_related_measures);

                auto mcmc_update_related_measures = update_formalism->template generate_mcmc_update_measures<SiteSystem,
                        SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>>(sp);
                this->concat_measures(mcmc_update_related_measures);

                auto site_update_related_measures = site_update->template generate_update_dynamics_measures<SiteSystem,
                        SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>>(sp);
                this->concat_measures(site_update_related_measures);

                auto common_defined_measures = this->generate_systembase_measures(sp.measures);
                this->concat_measures(common_defined_measures);
            }

            auto &get_update_formalism() {
                return *update_formalism;
            }

            json get_params_json() const {
                return sp.get_json();
            }

            T energy() const {
                return model->get_potential(site);
            }

            T drift_term() const {
                return model->get_drift_term(site);
            }

            void normalize(T &site_elem) {
                site_elem = model->normalize(site_elem);
            }

            double normalization_factor() {
                return update_formalism->get_normalization_factor(site);
            }

            typedef T SiteType;
            typedef ModelParameters ModelType;

        protected:
            const SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &sp;

            T site;

            std::unique_ptr<typename ModelParameters::Model> model;
            std::unique_ptr<typename UpdateFormalismParameters::MCMCUpdate> update_formalism;
            std::unique_ptr<typename SiteUpdateFormalismParameters::UpdateDynamics> site_update;

            void initialize_site();

            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SiteSystem>>>
            generate_site_system_measures(const json &measure_names);
        };


        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
        void
        SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>::initialize_site() {
            site = update_formalism->template random_state<T>();
        }

        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
        std::ostream &operator<<(std::ostream &os,
                                 const SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &site) {
            std::cout << site(0) << std::endl;
            return os;
        }


        template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
        std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>>>>
        SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>::generate_site_system_measures(
                const json &measure_names) {
            typedef SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> SiteSys;
            std::vector<std::unique_ptr<mcmc::common_measures::MeasurePolicy<SiteSys>>> site_measures{};
            for (auto &measure_name :  measure_names)
                if (measure_name == "Energy")
                    site_measures.push_back(std::make_unique<util::lattice_system_model_measures::MeasureEnergyPolicy<SiteSys>>());
                else if (measure_name == "Drift")
                    site_measures.push_back(std::make_unique<util::lattice_system_model_measures::MeasureDriftPolicy<SiteSys>>());
            return site_measures;
        }

    }
}


#endif //MAIN_SITE_HPP
