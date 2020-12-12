//
// Created by lukas on 12.12.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_ON_HPP
#define LATTICEMODELIMPLEMENTATIONS_ON_HPP


#include "../../link_lattice/link.hpp"
#include "mcmc_simulation/util/random.hpp"

template<typename T, uint N>
class ON : public Link<T>
{
protected:
    using Link<T>::x_;
public:
    explicit ON(const std::vector<T> x)
    {
        if(x.size() != N)
        {
            std::cout << "Invalid size of ON" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        x_ = x;
    }

    explicit ON(const double val)
    {
        x_ = std::vector<T>(dim(), val);
    }

    ON(const ON<T,N> &on)
    {
        x_ = on.x_;
    }

    uint dim() const
    {
        return N;
    }
};

template<typename T, uint N>
T operator*(const ON<T,N>& x, const ON<T,N>& y)
{
    T sum = 0;
    for(auto i = 0; i < N; i++)
        sum += x(i) * y(i);
    return sum;
}

template<typename T, uint N>
T operator*(const ON<T,N>& x, const double& y)
{
    T sum = 0;
    for(auto i = 0; i < N; i++)
        sum += x(i) * y;
    return sum;
}

namespace std
{
    template<typename T, uint N>
    double fabs(ON<T, N> arg)
    {
        return 0.0;
    }

    template<typename T, uint N>
    double abs(ON<T, N> arg)
    {
        return 0.0;
    }

    template<typename T, uint N>
    const double abs(const ON<T, N> arg)
    {
        return 0.0;
    }
}

template<typename T, uint N>
ON<T, N> operator+(const ON<T, N>& a,const ON<T, N>& b)
{
    auto val(a);
    val += b;
    return val;
}

namespace std
{
    template<typename T, uint N>
    std::string to_string(ON<T, N> x)
    {
        std::string conf = "";
        for(auto j = 0; j < x.dim(); j++)
            conf += std::to_string(x(j)) + " ";
        conf = conf.substr(0, conf.size() -1);
        return conf;
    }
}

/* namespace common_measures
{
    template<typename T>
    struct TypePolicy;

    template<typename T, uint N>
    struct TypePolicy< ON<T,N> > {
    public:
        static double realv(const ON<T,N> state) {
            return state.reduce();
        }

        static double imagv(const ON<T,N> state) {
            return 0.0; //state(0).imag();
        }

        static std::string conf(const ON<T, N> state) {
            std::string conf = "";
            for(auto j = 0; j < state.dim(); j++)
                conf += std::to_string(state(j)) + " ";
            conf = conf.substr(0, conf.size() -1);
            return conf;
        }
    };

    template<typename T, uint N>
    struct TypePolicy< ON<const T, N> > {
    public:
        static double realv(const ON<T,N> state) {
            return state.reduce();
        }

        static double imagv(const ON<T,N> state) {
            return 0.0; //state(0).imag();
        }

        static std::string conf(const ON<T, N> state) {
            std::string conf = "";
            for(auto j = 0; j < state.dim(); j++)
                conf += std::to_string(state(j)) + " ";
            conf = conf.substr(0, conf.size() -1);
            return conf;
        }
    };

    template<typename T, uint N>
    struct TypePolicy< const ON<T,N> > {
    public:
        static double realv(const ON<T,N> state) {
            return state.reduce();
        }

        static double imagv(const ON<T,N> state) {
            return 0.0; //state(0).imag();
        }

        static std::string conf(const ON<T, N> state) {
            std::string conf = "";
            for(uint j = 0; j < state.dim(); j++)
                conf += std::to_string(state(j)) + " ";
            conf = conf.substr(0, conf.size() -1);
            return conf;
        }
    };

    template<typename T, uint N>
    struct TypePolicy< const ON<const T, N> > {
    public:
        static double realv(const ON<T,N> state) {
            return state.reduce();
        }

        static double imagv(const ON<T,N> state) {
            return 0.0; //state(0).imag();
        }

        static std::string conf(const ON<T, N> state) {
            std::string conf = "";
            for(auto j = 0; j < state.dim(); j++)
                conf += std::to_string(state(j)) + " ";
            conf = conf.substr(0, conf.size() -1);
            return conf;
        }
    };
} */

#endif //LATTICEMODELIMPLEMENTATIONS_ON_HPP
