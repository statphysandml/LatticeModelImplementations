//
// Created by lukas on 09.10.19.
//

#ifndef MAIN_LATTICE_HPP
#define MAIN_LATTICE_HPP

#include "mcmc_simulation/header.hpp"

#include "../util/measures/lattice_measures.hpp"


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
                        *this, model_parameters->param_file_name(), rel_config_path_));

        update_parameters = std::make_unique<UpdateFormalismParameters>(
                generate_parameter_class_json<LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>, UpdateFormalismParameters> (
                        *this, update_parameters->param_file_name(), rel_config_path_));

        lattice_update_parameters = std::make_unique<LatticeUpdateFormalismParameters>(
                generate_parameter_class_json<LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>, LatticeUpdateFormalismParameters> (
                        *this, lattice_update_parameters->param_file_name(), rel_config_path_));
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
        params["measures"] = measures_;
    }

    void write_to_file(const std::string& rel_config_path) {
        // ToDo: Write function in parameters that accepts various objects and does the same operation as is done here for an arbitrary number of parameters -> parameters can even be collected in a vector of std::vector<*Parameteres>

        std::string model_params_path = get_value_by_key<std::string>(model_parameters->param_file_name() + "_path", rel_config_path);
        model_parameters->write_to_file(model_params_path);

        std::string update_formalism_params_path = get_value_by_key<std::string>(update_parameters->param_file_name() + "_path", rel_config_path);
        update_parameters->write_to_file(update_formalism_params_path);

        std::string lattice_update_params_path = get_value_by_key<std::string>(lattice_update_parameters->param_file_name() + "_path", rel_config_path);
        lattice_update_parameters->write_to_file(lattice_update_params_path);

        json model_parameters_ = model_parameters->get_json();
        delete_entry(model_parameters->param_file_name());

        json update_parameters_ = update_parameters->get_json();
        delete_entry(update_parameters->param_file_name());

        json lattice_update_parameters_ = lattice_update_parameters->get_json();
        delete_entry(lattice_update_parameters->param_file_name());

        Parameters::write_to_file(rel_config_path, param_file_name());

        add_entry(model_parameters->param_file_name(), model_parameters_);
        add_entry(update_parameters->param_file_name(), update_parameters_);
        add_entry(lattice_update_parameters->param_file_name(), lattice_update_parameters_);
    }

    Parameters build_expanded_raw_parameters() const
    {
        Parameters parameters(params);
        parameters.add_entry(model_parameters->param_file_name(), model_parameters->get_json());
        parameters.add_entry(update_parameters->param_file_name(), update_parameters->get_json());
        parameters.add_entry(lattice_update_parameters->param_file_name(), lattice_update_parameters->get_json());
        return parameters;
    }

    // ToDo: Needed? Add measures also as variable to SystemBaseParameters instead?
    json get_measure_names() const
    {
        return measures;
    }

    // Only used for testing in simulation.hpp
    typedef ModelParameters MP_;
    typedef UpdateFormalismParameters UP_;

protected:
    template<typename, typename, typename, typename>
    friend class LatticeSystem;

    json measures;

    std::unique_ptr<ModelParameters> model_parameters;
    std::unique_ptr<UpdateFormalismParameters> update_parameters;
    std::unique_ptr<LatticeUpdateFormalismParameters> lattice_update_parameters;

    uint16_t n_sites; // Total number of sites
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
        model = std::make_unique<typename ModelParameters::Model>(*lp.model_parameters);
        update_formalism = std::make_unique<typename UpdateFormalismParameters::MCMCUpdate>(*lp.update_parameters, *model);

        initialize_lattice();
        if(lp.lattice_action_type == "plaquette")
            set_plaquette_neighbours();
        else
            set_nearest_neighbours();

        lattice_update = std::make_unique<typename LatticeUpdateFormalismParameters::UpdateDynamics>(*lp.lattice_update_parameters);

        update_formalism->initialize(*this);
        lattice_update->initialize(*this);

        // Needs to be called at the end so that update objects can already be used!
        this->generate_measures();
    }

    void update_step(uint measure_interval)
    {
        lattice_update->operator()(*this, measure_interval);
    }

    void initialize(std::string starting_mode)
    {
        std::cout << "Inializte lattice not implemented" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    const T energy() const {
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

    void normalize(std::vector<T> &lattice_grid)
    {
        for(auto& elem : lattice_grid)
            elem = model->normalize(elem);
    }

    void generate_measures()
    {
        auto lattice_related_measures = generate_lattice_system_measures(lp.measures);
        this->concat_measures(lattice_related_measures);
        auto model_related_measures =  model->template generate_model_measures<LatticeSystem>(lp.measures);
        this->concat_measures(model_related_measures);
        auto lattice_update_related_measures =  lattice_update->template generate_update_dynamics_measures<LatticeSystem>(lp.measures);
        this->concat_measures(lattice_update_related_measures);
        auto common_defined_measures = common_measures::generate_measures<LatticeSystem>(lp.measures);
        this->concat_measures(common_defined_measures);
    }

    // Returns the total number of elements of the lattice - not the total number of sites
    const auto get_size() const
    {
        return lp.size;
    }

    const auto get_elem_per_site() const
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

    auto& get_neighbours()
    {
        return neighbours;
    }

    const auto get_neighbours() const
    {
        return neighbours;
    }

    auto& get_update_formalism()
    {
        return *update_formalism;
    }

    const auto get_system_representation() const
    {
        return lattice;
    }

    auto& get_system_representation()
    {
        return lattice;
    }

    const static std::string get_name()
    {
        return LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::name();
    }

    typedef T SiteType;
protected:
    const LatticeParameters<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> &lp;

    std::vector<T> lattice;
    std::vector< std::vector < T* > > neighbours;

    std::unique_ptr<typename ModelParameters::Model> model;
    std::unique_ptr<typename UpdateFormalismParameters::MCMCUpdate> update_formalism;
    std::unique_ptr<typename LatticeUpdateFormalismParameters::UpdateDynamics> lattice_update;

    void initialize_lattice();
    int neigh_dir(int n, int d, bool dir, int mu) const;
    void set_nearest_neighbours();
    void set_plaquette_neighbours();

    std::vector< std::unique_ptr<common_measures::MeasurePolicy<LatticeSystem>> > generate_lattice_system_measures(const json& measure_names);
};


template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
void LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::initialize_lattice() {
    lattice = std::vector<T> (get_size(), T(0));
    for(auto &site : lattice)
        site = update_formalism->template random_state<T>();
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
    // int offset;
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

template<typename T, typename ModelParameters, typename UpdateFormalismParameters, typename LatticeUpdateFormalismParameters>
std::vector< std::unique_ptr<common_measures::MeasurePolicy<LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>>>>
LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters>::generate_lattice_system_measures(const json& measure_names)
{
    typedef LatticeSystem<T, ModelParameters, UpdateFormalismParameters, LatticeUpdateFormalismParameters> LatSys;
    std::vector< std::unique_ptr<common_measures::MeasurePolicy<LatSys>> > lattice_measures {};
    for (auto& measure_name :  measure_names)
        if(measure_name == "Energy")
            lattice_measures.push_back(std::make_unique<lattice_model_measures::MeasureEnergyPolicy<LatSys>>());
        else if(measure_name == "EnergyImag")
            lattice_measures.push_back(std::make_unique<lattice_model_measures::MeasureEnergyImagPolicy<LatSys>>());
        else if(measure_name == "Drift")
            lattice_measures.push_back(std::make_unique<lattice_model_measures::MeasureDriftTermPolicy<LatSys>>());
        else if(measure_name == "WilsonAction")
            lattice_measures.push_back(std::make_unique<lattice_model_measures::MeasureWilsonActionPolicy<LatSys>>());
    return lattice_measures;
}


#endif //MAIN_LATTICE_HPP
