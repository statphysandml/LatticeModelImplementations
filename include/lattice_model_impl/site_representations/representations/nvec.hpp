//
// Created by lukas on 01.10.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_HPP


#include <cmath>

#include "../../link_lattice/link.hpp"
#include "mcmc_simulation/util/random.hpp"

template<typename T, uint N>
class NVec : public Link<T>
{
protected:
    using Link<T>::x_;
    std::uniform_int_distribution<int> uniform;
    std::normal_distribution<double> normal;
public:
    explicit NVec(const std::vector<T> x)
    {
        if(x.size() != N)
        {
            std::cout << "Invalid size of NVec" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        x_ = x;
        uniform = std::uniform_int_distribution<int>(0, dim() - 1);
        normal = std::normal_distribution<double>(0.0, 2.0);
    }

    NVec(const NVec<T,N> &nvec)
    {
        x_ = nvec.x_;
        uniform = std::uniform_int_distribution<int>(0, dim() - 1);
        normal = std::normal_distribution<double>(0.0, 2.0);
    }

    explicit NVec(const T x)
    {
        x_ = std::vector<T> (N, 0.0);
        x_[0] = x;
        uniform = std::uniform_int_distribution<int>(0, dim() - 1);
        normal = std::normal_distribution<double>(0.0, 2.0);
    }

    explicit NVec()
    {
        x_ = std::vector<T> (N, 0.0);
        uniform = std::uniform_int_distribution<int>(0, dim() - 1);
        normal = std::normal_distribution<double>(0.0, 2.0);
    }

    NVec(const T x, const std::vector<T> x_old)
    {
        if(x_old.size() != N)
        {
            std::cout << "Invalid size of NVec" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        x_ = x_old;
        x_[0] = x;
        uniform = std::uniform_int_distribution<int>(0, dim() - 1);
        normal = std::normal_distribution<double>(0.0, 2.0);
    }

    NVec(const T x, const NVec<T, N> x_old)
    {
        *this = x_old;
        x_[0] = x;
        uniform = std::uniform_int_distribution<int>(0, dim() - 1);
        normal = std::normal_distribution<double>(0.0, 2.0);
    }

    /* NVec& operator= (const NVec<T, N> &nvec)
    {
        x_ = nvec.x_;
        return *this;
    } */

    T reduce() const
    {
        double sum = 0;
        for (size_t i = 0; i < x_.size(); i++) sum += x_[i];
        return sum;
    }

    void rescale(const T val)
    {
        /* if(fabs(this->orig()) > val)
        {
            x_[0] = this->reduce();
            for(auto i = 1; i < this->dim(); i++)
                x_[i] = 0.0;
        } */
        auto resc = false;
        for (size_t i = 0; i < x_.size(); i++)
            if(fabs(x_[i]) > val)
                resc = true;
        if(resc)
        {
            remap();
        }
    }

    void remap()
    {
        /* auto j = 0; // uniform(gen);
        x_[j] = this->reduce();
        for (size_t i = 0; i < x_.size(); i++)
        {
            if(i != j)
                x_[i] = 0.0;
        } */

        auto j = uniform(gen);
        auto val = this->reduce();
        x_[j] = normal(gen);
        for (size_t i = 0; i < x_.size(); i++)
        {
            if(i != j)
                x_[i] = val - x_[j];
        }
    }

    T derivative() const
    {
        return 1.0;
    }

    T orig() const
    {
        return x_[0];
    }

    explicit operator double() const
    {
        return this->reduce();
    }

    uint dim() const
    {
        return N;
    }

    typedef T Ttype;
};

#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_HPP
