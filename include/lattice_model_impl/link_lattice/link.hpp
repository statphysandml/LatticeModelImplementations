//
// Created by lukas on 04.11.19.
//

#ifndef MAIN_LINK_HPP
#define MAIN_LINK_HPP


#include <iomanip>
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>

#include "mcmc_simulation/util/random.hpp"


template<typename T>
class Link
{
protected:
    //Link();
    //Link(double a, double b, double c, double d);
    //Link(double epsilon);
    std::vector< T > x_;
public:
    void Print() const;

    T operator()(int i) const;
    T& operator()(int i);

    virtual Link& operator+=(const Link& x);
    virtual Link& operator-=(const Link& x);

    /* virtual size_t size() const
    {
        return x_.size();
    } */

    //virtual T trace();

    virtual Link& adjungate();
};

// Link operator*(const Link& x, const Link& y);

template<typename T>
Link<T> operator-(const Link<T>& a,const Link<T>& b)
{
    auto val(a);
    val -= b;
    return val;
}

//std::ostream& operator<<(std::ostream &os, const Link& x);


template<typename T>
void Link<T>::Print() const
{
    std::cout<<"Link = " << x_[0];
    for(auto i = 1; i < x_.size(); i++)
        std::cout << ", " << x_[i];
    std::cout << std::endl;
}

template<typename T>
T Link<T>::operator()(int i) const
{
    return x_[i];
}

template<typename T>
T& Link<T>::operator()(int i)
{
    return x_[i];
}

template<typename T>
Link<T>& Link<T>::operator+=(const Link<T>& x) {
    for(size_t i = 0; i < x_.size(); i ++)
        x_[i] += x(i);
    return *this;
}

template<typename T>
Link<T>& Link<T>::operator-=(const Link<T>& x) {
    for(size_t i = 0; i < x_.size(); i ++)
        x_[i] -= x(i);
    return *this;
}

/*template<typename T>
T Link<T>::trace() {
    return x_[0]*2.0;
}*/

template<typename T>
Link<T>& Link<T>::adjungate() {
    x_[3] = -x_[3];
    x_[2] = -x_[2];
    x_[1] = -x_[1];
    return *this;
}

#endif //MAIN_LINK_HPP
