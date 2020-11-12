//
// Created by lukas on 12.11.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_NVEC_COMPLEX_HPP
#define LATTICEMODELIMPLEMENTATIONS_NVEC_COMPLEX_HPP

#include "nvec.hpp"

template<typename T, uint N, std::vector<bool>& real_>
class NVecComplex : public NVec<T, N>
{
private:
    using NVec<T, N>::x_;
public:
    explicit NVecComplex(const std::vector<T> x)
    {
        if(x.size() != N)
        {
            std::cout << "Invalid size of NVec" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        x_ = x;
    }

    NVecComplex(const NVecComplex<T,N, real_> &nvec)
    {
        x_ = nvec.x_;
    }

    explicit NVecComplex(const T x)
    {
        x_ = std::vector<T> (N, 0.0);
        x_[0] = x;
    }

    explicit NVecComplex()
    {
        x_ = std::vector<T> (N, 0.0);
    }

    NVecComplex(const T x, const std::vector<T> x_old)
    {
        if(x_old.size() != N)
        {
            std::cout << "Invalid size of NVec" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        x_ = x_old;
        x_[0] = x;
    }

    NVecComplex(const T x, const NVec<T, N> x_old)
    {
        *this = x_old;
        x_[0] = x;
    }

    /* NVec& operator= (const NVec<T, N> &nvec)
    {
        x_ = nvec.x_;
        return *this;
    } */

    std::complex<T> reduce() const
    {
        std::complex<T>sum(0.0, 0.0);
        for (size_t i = 0; i < x_.size(); i++)
        {
            if(real_[i])
                sum.real(sum.real() + x_[i]);
            else
                sum.imag(sum.imag() + x_[i]);
        }
        return sum;
    }

    void rescale(const T val)
    {
        if(fabs(x_[0]) > 10)
        {
            x_[0] = this->reduce().real();
            // for(auto i = 2; i < this->dim(); i++)
            x_[2] = 0.0;
        }

        if(fabs(x_[1]) > 10)
        {
            x_[1] = this->reduce().imag();
            for(auto i = 2; i < this->dim(); i++)
                x_[i] = 0.0;
        }
    }

    operator std::complex<T>() const
    {
        return this->reduce();
    }

    typedef T Ttype;
};


#endif //LATTICEMODELIMPLEMENTATIONS_NVEC_COMPLEX_HPP
