//
// Created by lukas on 11.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_DUMMY_MCMC_UPDATE_HPP
#define LATTICEMODELIMPLEMENTATIONS_DUMMY_MCMC_UPDATE_HPP


#include "mcmc_update_base.hpp"

class DummyMCMCUpdate;


class DummyMCMCUpdateParameters : public MCMCUpdateBaseParameters {
public:
    explicit DummyMCMCUpdateParameters(const json params_) : MCMCUpdateBaseParameters(params_)
    {}

    explicit DummyMCMCUpdateParameters() : DummyMCMCUpdateParameters(json{})
    {}

    static std::string name() {
        return "DummyMCMCUpdate";
    }

    typedef DummyMCMCUpdate MCMCUpdate;

protected:
    friend class DummyMCMCUpdate;
};


class DummyMCMCUpdate : public MCMCUpdateBase< DummyMCMCUpdate >
{
public:
    template<typename Model>
    explicit DummyMCMCUpdate(const DummyMCMCUpdateParameters &up_, Model & model_) : up(up_)
    {}

    template<typename T>
    T operator() (const T site)
    {
        return site;
    }

    template<typename T>
    T operator() (const T site, const std::vector< T* > neighbours)
    {
        return site;
    }

protected:
    const DummyMCMCUpdateParameters & up;
};

#endif //LATTICEMODELIMPLEMENTATIONS_DUMMY_MCMC_UPDATE_HPP
