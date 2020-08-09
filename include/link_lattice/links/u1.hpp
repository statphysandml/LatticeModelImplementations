//
// Created by lukas on 04.11.19.
//

#ifndef MAIN_U1_HPP
#define MAIN_U1_HPP

#include <iomanip>
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>

#include "../link.hpp"

template<typename T>
class U1 : public Link<T>
{
private:
    using Link<T>::x_;
public:
    //U1();
    U1(T a);
    U1(double epsilon);
    U1(std::string init = "random");

    //void Print() const;

    //double operator()(int i) const;

    // U1& operator+=(const U1& x);
    // U1& operator-=(const U1& x);

    T trace();

    U1& adjungate();

    /* size_t size() const
    {
        return 1;
    } */

    //private:
    //  std::vector< double > x_;
};

template<typename T>
U1<T> operator*(const U1<T>& x, const U1<T>& y);

template<typename T>
U1<T> operator-(const U1<T>& a,const U1<T>& b);

template<typename T>
std::ostream& operator<<(std::ostream &os, const U1<T>& x);

template <typename T>
U1<T> operator*(const U1<T>& x, const U1<T>& y)
{
    return U1<T>(x(0)*y(0));
}

template <typename T>
U1<T> operator-(const U1<T>& a,const U1<T>& b)
{
    U1<T> temp(a);
    temp -= b;
    return temp;
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const U1<T>& x) {
    os << x(0);
    return os;
}

/*std::complex<double> operator-(const int& a,const std::complex<double>& b) {
    return std::complex<double> (a,0) - b;
}*/


template <typename T>
U1<T>::U1(std::string init)
{
    if(init == "null") {
        x_.push_back(T(0,0));
    }
    else if(init == "identity") {
        x_.push_back(T(1,0));
    }
    else {
        std::uniform_real_distribution<double> uniform(-1,1);
        x_.push_back(T(uniform(gen),uniform(gen)));
        x_[0] = x_[0]/std::abs(x_[0]);
    }
}

/* template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
} */

template <typename T>
U1<T>::U1(double epsilon) {//######################################################
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> uniform(-epsilon,epsilon);

    x_.push_back(T(1,uniform(gen)));

    x_[0] = x_[0]/std::abs(x_[0]);

    /*double length;
    for(auto i = 1; i < 4; i++) {
        x_.push_back(uniform(gen));
        length += x_[i]*x_[i];
    }
    length = sqrt(length);
    for(int i = 1; i < 4; i++) x_[i] = epsilon*x_[i]/length;*/
};

template <typename T>
U1<T>::U1(T a) {
    x_.push_back(a);
}

/* template <typename T>
U1<T>& U1<T>::operator+=(const U1<T>& x) {
    x_[0] += x(0);
    return *this;
}

template <typename T>
U1<T>& U1<T>::operator-=(const U1<T>& x) {
    x_[0] -= x(0);
    return *this;
} */

template <typename T>
T U1<T>::trace() {
    return x_[0];
}

template <typename T>
U1<T>& U1<T>::adjungate() {
    x_[0] = std::conj(x_[0]);
    return *this;
}


#endif //MAIN_U1_HPP
