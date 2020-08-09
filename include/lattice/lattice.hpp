//
// Created by lukas on 09.10.19.
//

#ifndef MAIN_LATTICE_HPP
#define MAIN_LATTICE_HPP

#include "mcmc_simulation/header.hpp"


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
class LatticeSystem;

template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
class LatticeParameters : public SystemBaseParameters {
public:
    explicit LatticeParameters(const json params_, const std::string rel_config_path_) : SystemBaseParameters(params_) {

        measures = get_value_by_key<json>("measures", {"Mean"});
        dimensions = get_value_by_key< std::vector<int> >("dimensions", std::vector<int> {4, 4});
        lattice_action_type = get_value_by_key<std::string>("lattice_action_type", "nearest_neighbour");

        n_sites = 1;
        dim_mul.push_back(n_sites);
        for(auto dim: dimensions) {
            n_sites *= dim;
            dim_mul.push_back(n_sites);
        }
        dimension = dimensions.size();

        if(lattice_action_type == "nearest_neighbour")
            elem_per_site = 1;
        else
            elem_per_site = dimension;
        size = n_sites * elem_per_site;

        model_parameters = std::make_unique<ModelParameters>(
                generate_parameter_class_json<LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>, ModelParameters> (
                        *this, "model", rel_config_path_));

        update_parameters = std::make_unique<UpdateFormalismParameters>(
                generate_parameter_class_json<LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>, UpdateFormalismParameters> (
                        *this, "update", rel_config_path_));

        lattice_update_parameters = std::make_unique<LatticeUpdateFormalismParameters>(
                generate_parameter_class_json<LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>, LatticeUpdateFormalismParameters> (
                        *this, "lattice_update", rel_config_path_));
    }

    const std::vector<int>& get_dimensions() const
    {
        return dimensions;
    }

    const uint& get_elems_per_site() const
    {
        return elem_per_site;
    }

    typedef LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> System;

    const static std::string name() {
        return "Lattice";
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

        std::string lattice_update_params_path = get_value_by_key<std::string>("lattice_update_params_path", rel_config_path);
        lattice_update_parameters->write_to_file(lattice_update_params_path);

        Parameters::write_to_file(rel_config_path, "systembase_params");
    }

    Parameters build_expanded_raw_parameters() const
    {
        Parameters parameters(params);
        parameters.add_entry("model", model_parameters->get_json());
        parameters.add_entry("update", update_parameters->get_json());
        parameters.add_entry("lattice_update", lattice_update_parameters->get_json());
        return parameters;
    }

    // Only used for testing in simulation.hpp
    typedef ModelParameters MP_;
    typedef UpdateFormalismParameters UP_;

private:
    template<typename, typename, typename, typename>
    friend class LatticeSystem;

    json measures;

    std::unique_ptr<ModelParameters> model_parameters;
    std::unique_ptr<UpdateFormalismParameters> update_parameters;
    std::unique_ptr<LatticeUpdateFormalismParameters> lattice_update_parameters;

    uint16_t n_sites; // Total number sites
    uint16_t size; // Total number of elements on lattice
    std::vector<int> dimensions; // Different dimensions
    std::vector<int> dim_mul; // Accumulated different dimensions (by product)
    int dimension; // Number of dimensions
    uint elem_per_site; // Number of elements per site (is equal to dimension), but only for the link model

    std::string lattice_action_type;
};


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
class LatticeSystem : public SystemBase< LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> > {
public:
    explicit LatticeSystem(const LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &lp_) : lp(lp_) {
        this->generate_measures();

        model = std::make_unique<typename ModelParameters::Model>(*lp.model_parameters);
        update_formalism = std::make_unique<typename UpdateFormalismParameters::UpdateFormalism>(*lp.update_parameters, *model);

        initialize_lattice();
        if(lp.lattice_action_type == "plaquette")
            set_plaquette_neighbours();
        else
            set_nearest_neighbours();

        lattice_update = std::make_unique<typename LatticeUpdateFormalismParameters::LatticeUpdate>(*lp.lattice_update_parameters);
        lattice_update->initialize(*this);

    }

    void update_step(uint measure_interval)
    {
        lattice_update->operator()(*this, measure_interval);
    }

    T energy() const {
        T energy(0);
        for(uint i = 0; i < get_size(); i++) {
            energy += model->get_potential(lattice[i], neighbours[i]);
        }
        // return  2 * energy / (double(size()) * 6.0);
        return 0.5 * energy / double(get_size()); // 2 *
    }

    T drift_term() const {
        T drift_term(0);
        for(uint i = 0; i < get_size(); i++) {
            // ToDo: Think about how this can be integrated!
            // drift_term += model->get_drift_term(lattice[i], neighbours[i]);
        }
        // return  2 * energy / (double(get_size()) * 6.0);
        return 0.5 * drift_term / double(get_size()); // 2 *
    }


    struct MeasureEnergyPolicy: public MeasurePolicy< LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> > {
        public:
        std::string measure(const LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &system) override {
                return std::to_string(TypePolicy<double>::realv(system.energy()));
        }

        std::string name()
        {
            return "Energy";
        }
    };

    struct MeasureEnergyImagPolicy: public MeasurePolicy< LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> > {
        public:
        std::string measure(const LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &system) override {
                return std::to_string(TypePolicy<double>::imagv(system.energy()));
        }

        std::string name()
        {
            return "EnergyImag";
        }
    };

    struct MeasureDriftTermPolicy: public MeasurePolicy< LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> > {
    public:
        std::string measure(const LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &system) override {
            auto drift_term = system.drift_term();
            return std::to_string(TypePolicy<decltype(drift_term)>::realv(drift_term)) + " " + std::to_string(TypePolicy<decltype(drift_term)>::imagv(drift_term));
        }

        std::string name()
        {
            return "Drift";
        }
    };

    struct MeasureWilsonActionPolicy: public MeasurePolicy< LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> > {
    public:
        std::string measure(const LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &system) override {
            return std::to_string(TypePolicy<double>::realv(system.energy()));
        }

        std::string name()
        {
            return "WilsonAction";
        }
    };

    void generate_measures()
    {
        measures = std::vector< std::unique_ptr<MeasurePolicy<LatticeSystem>> > {};
        for (auto& element :  lp.measures)
        {
            if(element == "Energy")
                measures.push_back(std::make_unique<MeasureEnergyPolicy>());
            else if(element == "EnergyImag")
                measures.push_back(std::make_unique<MeasureEnergyImagPolicy>());
            else if(element == "Drift")
                measures.push_back(std::make_unique<MeasureDriftTermPolicy>());
            else if(element == "WilsonAction")
                measures.push_back(std::make_unique<MeasureWilsonActionPolicy>());
            else {
                /* auto model_policy_measure = std::make_unique(model->template model_measure_factory<LatticeSystem<T, ModelParameters, UpdateFormalismParameters>, LatticeParameters<T, ModelParameters, UpdateFormalismParameters> >(
                        element, lp); */
                // if (model_policy_measure == nullptr)
                    measures.push_back(unique_measure_factory< LatticeSystem>(element));
                /* else
                    measures.push-_back(model_policy_measure); */
            }
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

    // Returns the total number of sites
    const auto get_size() const
    {
        return lp.size;
    }

    const auto elem_per_site() const
    {
        return lp.elem_per_site;
    }

    auto& at(int i)
    {
        return lattice[i];
    }

    const auto at(int i) const
    {
        return lattice[i];
    }

    auto& neighbours_at(int i)
    {
        return neighbours[i];
    }

    const auto neighbours_at(int i) const
    {
        return neighbours[i];
    }

    const auto get_neighbours() const
    {
        return neighbours;
    }

    auto& get_update_formalism()
    {
        return *update_formalism;
    }

    auto& get_lattice_grid()
    {
        return lattice;
    }

    typedef T SiteType;
private:
    const LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &lp;

    std::vector<T> lattice;
    std::vector< std::vector < T* > > neighbours;
    std::vector< std::unique_ptr<MeasurePolicy<LatticeSystem>> > measures;

    std::unique_ptr<typename ModelParameters::Model> model;
    std::unique_ptr<typename UpdateFormalismParameters::UpdateFormalism> update_formalism;
    std::unique_ptr<typename LatticeUpdateFormalismParameters::LatticeUpdate> lattice_update;

    void initialize_lattice();
    int neigh_dir(int n, int d, bool dir, int mu) const;
    void set_nearest_neighbours();
    void set_plaquette_neighbours();
};


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
void LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::initialize_lattice() {
    lattice = std::vector<T> (get_size(), T(0));
    for(auto &site : lattice)
        site = model->template random_state<T>();
}


//site, moving dimension, direction, index
template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
int LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::neigh_dir(int n, int d, bool dir, int mu) const {
    if(dir)
        return (n-n%(lp.dim_mul[d]*lp.dimensions[d])+(n+lp.dim_mul[d])%(lp.dim_mul[d]*lp.dimensions[d]))*lp.elem_per_site+mu;
    else
        return (n-n%(lp.dim_mul[d]*lp.dimensions[d])+(n-lp.dim_mul[d]+lp.dim_mul[d]*lp.dimensions[d])%(lp.dim_mul[d]*lp.dimensions[d]))*lp.elem_per_site+mu;
}


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
void LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::set_nearest_neighbours() {
    int offset;
    for(uint i = 0; i < get_size(); i++) {
        // offset = 1;
        std::vector < T* > nn_of_site;
        // std::cout << "i: " << i << std::endl;
        for(uint d = 0; d < lp.dimensions.size(); d++) {
            //std::cout << i-i%(offset*lp.dimensions[d])+(i+offset)%(offset*lp.dimensions[d]) << " - " << i-i%(offset*lp.dimensions[d])+(i-offset+offset*lp.dimensions[d])%(offset*lp.dimensions[d]) << std::endl;
            // x + nu
            nn_of_site.push_back(&lattice[neigh_dir(i, d, true, 0)]); // i-i%(offset*lp.dimensions[d])+(i+offset)%(offset*lp.dimensions[d])]);
            // nn_of_site.push_back(&lattice[i-i%(offset*lp.dimensions[d])+(i+offset)%(offset*lp.dimensions[d])]);
            // x - nu
            nn_of_site.push_back(&lattice[neigh_dir(i, d, false, 0)]);
            // nn_of_site.push_back(&lattice[i-i%(offset*lp.dimensions[d])+(i-offset+offset*lp.dimensions[d])%(offset*lp.dimensions[d])]);
            // offset = offset*lp.dimensions[d];
        }
        neighbours.push_back(nn_of_site);
    }
}


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
void LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::set_plaquette_neighbours() {
    // Loop over all sites
    for(uint n = 0; n < get_size()/lp.elem_per_site; n++) {
        // Loop over all links for a given site
        for(auto mu = 0; mu < lp.dimension; mu++) {
            std::vector < T* > nn_of_site;
            // Loop over possible plaquettes (left and right to the respective dimension)
            for(auto nu = 0; nu < lp.dimension; nu++) {
                if (nu != mu) {
                    /*lat.cumneigh_[lat.neigh_dir(n,mu,true,nu)]++;
                    lat.cumneigh_[lat.neigh_dir(n,nu,true,mu)]++;
                    lat.cumneigh_[n*lat.dim()+nu]++;
                    lat.cumneigh_[lat.neigh_dir(lat.neigh_dir(n,mu,true,0)/lat.dim(),nu,false,nu)]++;
                    lat.cumneigh_[lat.neigh_dir(n,nu,true,mu)]++;
                    lat.cumneigh_[lat.neigh_dir(n,nu,false,nu)]++;*/

                    nn_of_site.push_back(&lattice[neigh_dir(n, mu,true, nu)]);
                    nn_of_site.push_back(&lattice[neigh_dir(n, nu, true, mu)]);
                    nn_of_site.push_back(&lattice[n * lp.elem_per_site + nu]);
                    nn_of_site.push_back(&lattice[neigh_dir(neigh_dir(n, mu, true, 0) / lp.elem_per_site, nu, false, nu)]);
                    nn_of_site.push_back(&lattice[neigh_dir(n,nu,false,mu)]);
                    nn_of_site.push_back(&lattice[neigh_dir(n,nu,false,nu)]);
//                    A += lat(lat.neigh_dir(n,mu,true,nu))*((lat(lat.neigh_dir(n,nu,true,mu))).adjungate())*(lat(n*lat.dim()+nu).adjungate());
//                    if(both_orientations) A += ((lat(lat.neigh_dir(lat.neigh_dir(n,mu,true,0)/lat.dim(),nu,false,nu))).adjungate())*((lat(lat.neigh_dir(n,nu,false,mu))).adjungate())*lat(lat.neigh_dir(n,nu,false,nu));
                }
            }
            neighbours.push_back(nn_of_site);
        }
    }
}


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
std::ostream& operator<<(std::ostream &os, const LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &lattice) {
    for(uint i = 0; i < lattice.length(); i++) {
        std::cout << lattice(i) << " ";
        //if((i+1)%30 == 0) std::cout << std::endl;
    }
    std::cout<<std::endl;
    return os;
}


#endif //MAIN_LATTICE_HPP
