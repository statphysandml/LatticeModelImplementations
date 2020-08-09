//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_SITE_HPP
#define MAIN_SITE_HPP

#include "mcmc_simulation/header.hpp"


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
class SiteSystem;


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
class SiteParameters : public SystemBaseParameters {
public:
    explicit SiteParameters(const json params_, const std::string rel_config_path_) : SystemBaseParameters(params_) {
        measures = get_value_by_key<json>("measures");

        model_parameters = std::make_unique<ModelParameters>(
                generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, ModelParameters> (
                        *this, "model", rel_config_path_));

        update_parameters = std::make_unique<UpdateFormalismParameters>(
                generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, UpdateFormalismParameters> (
                        *this, "update", rel_config_path_));

        site_update_parameters = std::make_unique<SiteUpdateFormalismParameters>(
                generate_parameter_class_json<SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters>, SiteUpdateFormalismParameters> (
                        *this, "site_update", rel_config_path_));
    }

    typedef SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> System;

    static std::string name() {
        return "Site";
    }

    std::unique_ptr<System> generate() {
        return std::make_unique<System>(*this);
    }

    void update_measures(const json& measures_)
    {
        measures = measures_;
    }

    void write_to_file(const std::string& rel_config_path) const {
        std::string model_params_path = get_value_by_key<std::string>("model_params_path", rel_config_path);
        model_parameters->write_to_file(model_params_path);

        std::string update_formalism_params_path = get_value_by_key<std::string>("update_params_path", rel_config_path);
        update_parameters->write_to_file(update_formalism_params_path);

        std::string lattice_update_params_path = get_value_by_key<std::string>("site_update_params_path", rel_config_path);
        site_update_parameters->write_to_file(lattice_update_params_path);

        Parameters::write_to_file(rel_config_path, "systembase_params");
    }

    Parameters build_expanded_raw_parameters() const
    {
        Parameters parameters(params);
        parameters.add_entry("model", model_parameters->get_json());
        parameters.add_entry("update", update_parameters->get_json());
        parameters.add_entry("site_update", site_update_parameters->get_json());
        return parameters;
    }

    // Only used for testing in simulation.hpp
    typedef ModelParameters MP_;
    typedef UpdateFormalismParameters UP_;

private:
    template<typename, typename, typename>
    friend class SiteSystem;

    json measures;

    std::unique_ptr<ModelParameters> model_parameters;
    std::unique_ptr<UpdateFormalismParameters> update_parameters;
    std::unique_ptr<SiteUpdateFormalismParameters> site_update_parameters;
};

template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename SiteUpdateFormalismParameters>
class SiteSystem : SystemBase< SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > {
public:
    explicit SiteSystem(const SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &sp_) : sp(sp_) {
        this->generate_measures();

        model = std::make_unique<typename ModelParameters::Model>(*sp.model_parameters);
        update_formalism = std::make_unique<typename UpdateFormalismParameters::UpdateFormalism>(*sp.update_parameters, *model);

        initialize_site();

        site_update = std::make_unique<typename SiteUpdateFormalismParameters::SiteUpdate>(*sp.site_update_parameters);
        site_update->initialize(*this);
    }

    void update_step(uint measure_interval)
    {
        site_update->operator()(*this, measure_interval);
    }

    T energy() {
        return model->get_potential(site);
    }

    T drift_term() {
        return model->get_drift_term(site);
    }

    double normalization_factor() {
        return update_formalism->get_normalization_factor(site);
    }

    struct MeasureEnergyPolicy: public MeasurePolicy< SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > {
    public:
        std::string measure(SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &system) override {
            return std::to_string(TypePolicy<T>::realv(system.energy()));
        }

        std::string name()
        {
            return "Energy";
        }
    };

    struct MeasureDriftTermPolicy: public MeasurePolicy< SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > {
    public:
        std::string measure(SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &system) override {
            auto drift_term = system.drift_term();
            return std::to_string(TypePolicy<decltype(drift_term)>::realv(drift_term)) + " " + std::to_string(TypePolicy<decltype(drift_term)>::imagv(drift_term));
        }

        std::string name()
        {
            return "Drift";
        }
    };

    struct MeasureEnergyImagPolicy: public MeasurePolicy< SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > {
    public:
        std::string measure(SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &system) override {
            return std::to_string(TypePolicy<T>::imagv(system.energy()));
        }

        std::string name()
        {
            return "EnergyImag";
        }
    };

    struct MeasureNormalizationPolicy: public MeasurePolicy< SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> > {
    public:
        std::string measure(SiteSystem<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &system) override {
            return std::to_string(system.normalization_factor());
        }

        std::string name()
        {
            return "Normalization";
        }
    };

    void generate_measures()
    {
        measures = std::vector< std::unique_ptr<MeasurePolicy<SiteSystem>> > {};
        for (auto& element :  sp.measures)
        {
            if(element == "Energy")
                measures.push_back(std::make_unique<MeasureEnergyPolicy>());
            else if(element == "EnergyImag")
                measures.push_back(std::make_unique<MeasureEnergyImagPolicy>());
            else if(element == "Drift")
                measures.push_back(std::make_unique<MeasureDriftTermPolicy>());
            else if(element == "Normalization")
                measures.push_back(std::make_unique<MeasureNormalizationPolicy>());
            else
                measures.push_back(unique_measure_factory< SiteSystem>(element));
        }
    }

    std::vector<std::string> perform_measure()
    {
        std::vector<std::string> results;
        for(auto const& element: measures) {
            results.push_back(element->measure(*this));
        }
        return results;
    }

    std::vector<std::string> get_measure_names()
    {
        std::vector<std::string> results;
        for(auto const& element: measures) {
            results.push_back(element->name());
        }
        return results;
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

    auto& get_site()
    {
        return site;
    }

    const auto get_site() const
    {
        return site;
    }

    auto& get_update_formalism()
    {
        return *update_formalism;
    }

    typedef T SiteType;
private:
    const SiteParameters<T, ModelParameters, UpdateFormalismParameters, SiteUpdateFormalismParameters> &sp;

    T site;
    std::vector< std::unique_ptr<MeasurePolicy<SiteSystem>> > measures;

    std::unique_ptr<typename ModelParameters::Model> model;
    std::unique_ptr<typename UpdateFormalismParameters::UpdateFormalism> update_formalism;
    std::unique_ptr<typename SiteUpdateFormalismParameters::LatticeUpdate> site_update;

    // For adaptive step size
    double KExpectation = 0;
    int thermalization_counter = 0;

    void initialize_site();
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

#endif //MAIN_SITE_HPP
