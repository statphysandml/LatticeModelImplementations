//
// Created by lukas on 06.08.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP
#define LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP


#include "../global_update.hpp"

template<typename ModelParameters>
class HybridMonteCarloUpdate;


template<typename ModelParameters>
class HybridMonteCarloUpdateParameters : public GlobalUpdateFormalismParameters {
public:
    explicit HybridMonteCarloUpdateParameters(const json params_) : GlobalUpdateFormalismParameters(params_),
                                                                    epsilon(get_value_by_key<double>("epsilon")),
                                                                    n(get_value_by_key<int>("n"))

    {}

    explicit HybridMonteCarloUpdateParameters(const double epsilon_, const int n_
    ) : HybridMonteCarloUpdateParameters(json{
        {"epsilon", epsilon_},
        {"n", n_}})
    {}

    static std::string name() {
        return "HybridMonteCarloUpdate";
    }

    typedef HybridMonteCarloUpdate<ModelParameters> UpdateFormalism;

protected:
    friend class HybridMonteCarloUpdate<ModelParameters>;

    const double epsilon;
    const int n;
};


template<typename ModelParameters>
class HybridMonteCarloUpdate // : public SingleGlobalUpdateFormalism< HybridMonteCarloUpdate<ModelParameters> >
{
public:
    explicit HybridMonteCarloUpdate(const HybridMonteCarloUpdateParameters<ModelParameters> &up_, typename ModelParameters::Model & model_) : up(up_), model(model_)
    {
        normal = std::normal_distribution<double>(0.0, 1.0);
    }

    template<typename Lattice>
    void initialize_update(const Lattice& lattice)
    {
        momenta = std::vector<double> (lattice.get_elems_per_site(), 0.0);
    }

    template<typename T>
    void operator() (std::vector<T>& lattice, const std::vector< std::vector<T*> >& neighbours)
    {
        // Sample momenta
        std::generate(momenta.begin(), momenta.end(), [this] () { return normal(gen); });

        std::cout << "ToDo: implement Leapfrog update with boost::odeint" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // using SingleGlobalUpdateFormalism< HybridMonteCarloUpdate<ModelParameters> >::operator();

protected:
    const HybridMonteCarloUpdateParameters<ModelParameters> & up;
    typename ModelParameters::Model & model;

    std::vector<double> momenta;
    std::normal_distribution<double> normal;
};

#endif //LATTICEMODELIMPLEMENTATIONS_HYBRID_MONTE_CARLO_HPP
