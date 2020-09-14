//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_SITE_HPP
#define MAIN_SITE_HPP

#include "mcmc_simulation/header.hpp"

#include "../measures/lattice_measures.hpp"
#include "../update_dynamics/site_update_dynamics/site_update_formalisms/simple_update.hpp"

template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters=SiteSimpleUpdateParameters>
class SiteSystem;


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters=SiteSimpleUpdateParameters>
class SiteParameters : public SystemBaseParameters {
public:
    explicit SiteParameters(const json params_, const std::string rel_config_path_) : SystemBaseParameters(params_) {
        measures = get_value_by_key<json>("measures", {});

        model_parameters = std::make_unique<ModelParameters>(
                generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, ModelParameters> (
                        *this, model_parameters->param_file_name(), rel_config_path_));

        update_parameters = std::make_unique<UpdateFormalismParameters>(
                generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, UpdateFormalismParameters> (
                        *this, update_parameters->param_file_name(), rel_config_path_));

        site_update_parameters = std::make_unique<SiteUpdateFormalismParameters>(
                generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, SiteUpdateFormalismParameters> (
                        *this, site_update_parameters->param_file_name(), rel_config_path_));
    }

    typedef SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> System;

    const static std::string name() {
        return "Site";
    }

    std::unique_ptr<System> generate() {
        return std::make_unique<System>(*this);
    }

    void update_measures(const json& measures_)
    {
        measures = measures_;
    }

    void write_to_file(const std::string& rel_config_path) {
        std::string model_params_path = get_value_by_key<std::string>(model_parameters->param_file_name() + "_path", rel_config_path);
        model_parameters->write_to_file(model_params_path);

        std::string update_formalism_params_path = get_value_by_key<std::string>(update_parameters->param_file_name() + "_path", rel_config_path);
        update_parameters->write_to_file(update_formalism_params_path);

        std::string site_update_params_path = get_value_by_key<std::string>(site_update_parameters->param_file_name() + "_path", rel_config_path);
        site_update_parameters->write_to_file(site_update_params_path);

        json model_parameters_ = model_parameters->get_json();
        delete_entry(model_parameters->param_file_name());

        json update_parameters_ = update_parameters->get_json();
        delete_entry(update_parameters->param_file_name());

        json site_update_parameters_ = site_update_parameters->get_json();
        delete_entry(site_update_parameters->param_file_name());

        Parameters::write_to_file(rel_config_path, param_file_name());

        add_entry(model_parameters->param_file_name(), model_parameters_);
        add_entry(update_parameters->param_file_name(), update_parameters_);
        add_entry(site_update_parameters->param_file_name(), site_update_parameters_);
    }

    Parameters build_expanded_raw_parameters() const
    {
        Parameters parameters(params);
        parameters.add_entry(model_parameters->param_file_name(), model_parameters->get_json());
        parameters.add_entry(update_parameters->param_file_name(), update_parameters->get_json());
        parameters.add_entry(site_update_parameters->param_file_name(), site_update_parameters->get_json());
        return parameters;
    }

    json get_measure_names() const
    {
        return measures;
    }

    // Only used for testing in simulation.hpp
    typedef ModelParameters MP_;
    typedef UpdateFormalismParameters UP_;

private:
    template<typename, typename, typename, typename>
    friend class SiteSystem;

    json measures;

    std::unique_ptr<ModelParameters> model_parameters;
    std::unique_ptr<UpdateFormalismParameters> update_parameters;
    std::unique_ptr<SiteUpdateFormalismParameters> site_update_parameters;
};

template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
class SiteSystem : public SystemBase< SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > {
public:
    explicit SiteSystem(const SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &sp_) : sp(sp_) {
        model = std::make_unique<typename ModelParameters::Model>(*sp.model_parameters);
        update_formalism = std::make_unique<typename UpdateFormalismParameters::MCMCUpdate>(*sp.update_parameters, *model);

        initialize_site();

        site_update = std::make_unique<typename SiteUpdateFormalismParameters::UpdateDynamics>(*sp.site_update_parameters);

        update_formalism->initialize(*this);
        site_update->initialize(*this);

        // Needs to be called at the end so that update objects can already be used!
        this->generate_measures();
    }

    static SiteSystem from_json(const json params_, const std::string rel_config_path_)
    {
        SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> sp_(params_, rel_config_path_);
        return SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>(sp_);
    }

    void update_step(uint measure_interval)
    {
        site_update->operator()(*this, measure_interval);
    }

    T energy() const {
        return model->get_potential(site);
    }

    T drift_term() const {
        return model->get_drift_term(site);
    }

    void normalize(T &site_elem)
    {
        site_elem = model->normalize(site_elem);
    }

    double normalization_factor() {
        return update_formalism->get_normalization_factor(site);
    }

    /* struct MeasureNormalizationPolicy: public common_measures::MeasurePolicy< SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > {
    public:
        std::string measure(const SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &system) override {
            return std::to_string(system.normalization_factor());
        }

        std::string name()
        {
            return "Normalization";
        }
    }; */

    void generate_measures()
    {
        auto lattice_related_measures = generate_site_system_measures(sp.measures);
        this->concat_measures(lattice_related_measures);
        auto model_related_measures =  model->template generate_model_measures<SiteSystem>(sp.measures);
        this->concat_measures(model_related_measures);
        auto site_update_related_measures =  site_update->template generate_update_dynamics_measures<SiteSystem>(sp.measures);
        this->concat_measures(site_update_related_measures);
        auto common_defined_measures = common_measures::generate_measures<SiteSystem>(sp.measures);
        this->concat_measures(common_defined_measures);
    }

    const auto get_size() const
    {
        return 1;
    }

    auto& at(int i)
    {
        return site;
    }

    const auto at(int i) const
    {
        return site;
    }

    auto& get_system_representation()
    {
        return site;
    }

    const auto get_system_representation() const
    {
        return site;
    }

    auto& get_update_formalism()
    {
        return *update_formalism;
    }

    const static std::string get_name()
    {
        return SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>::name();
    }

    json get_params_json() const
    {
        return sp.get_json();
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

    std::vector< std::unique_ptr<common_measures::MeasurePolicy<SiteSystem>> > generate_site_system_measures(const json& measure_names);
};


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
void SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>::initialize_site() {
    site = model->template random_state<T>();
}

template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
std::ostream& operator<<(std::ostream &os, const SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &site) {
    std::cout << site(0) << std::endl;
    return os;
}


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
std::vector< std::unique_ptr<common_measures::MeasurePolicy<SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>>>>
SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>::generate_site_system_measures(const json& measure_names)
{
    typedef SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> SiteSys;
    std::vector< std::unique_ptr<common_measures::MeasurePolicy<SiteSys>> > site_measures {};
    for (auto& measure_name :  measure_names)
        if(measure_name == "Energy")
            site_measures.push_back(std::make_unique<lattice_model_measures::MeasureEnergyPolicy<SiteSys>>());
        else if(measure_name == "EnergyImag")
            site_measures.push_back(std::make_unique<lattice_model_measures::MeasureEnergyImagPolicy<SiteSys>>());
        else if(measure_name == "Drift")
            site_measures.push_back(std::make_unique<lattice_model_measures::MeasureDriftTermPolicy<SiteSys>>());
        else if(measure_name == "WilsonAction")
            site_measures.push_back(std::make_unique<lattice_model_measures::MeasureWilsonActionPolicy<SiteSys>>());
    return site_measures;
}


#endif //MAIN_SITE_HPP
