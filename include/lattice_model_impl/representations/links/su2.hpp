//
// Created by lukas on 05.11.19.
//

#ifndef MAIN_SU2_HPP
#define MAIN_SU2_HPP

#include <iomanip>
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>

#include "../link.hpp"


double compute_x0(double beta, double a);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T>
class SU2 : public Link<T>
{
private:
    using Link<T>::x_;
public:
    SU2(double a, double b, double c, double d);
    explicit SU2(double epsilon);
    explicit SU2(std::string init = "random");
    SU2(const SU2& A, double beta);

    double trace();
    double det();
    SU2& adjungate();
};

template <typename T>
SU2<T>::SU2(double a, double b, double c, double d) {
    x_.push_back(a);
    x_.push_back(b);
    x_.push_back(c);
    x_.push_back(d);
}

template <typename T>
SU2<T>::SU2(double epsilon) {
    // ( proposal state according to equation (4.30))
    std::uniform_real_distribution<double> distribution(-0.5,0.5);

    x_.push_back(sgn(distribution(mcmc::util::gen))*sqrt(1-epsilon*epsilon));

    double length = 0;
    for(auto i = 1; i < 4; i++) {
        x_.push_back(distribution(mcmc::util::gen));
        length += x_[i]*x_[i];
    }
    length = sqrt(length);
    for(int i = 1; i < 4; i++) x_[i] = epsilon*x_[i]/length;
};

template <typename T>
SU2<T>::SU2(std::string init)
{
    if(init == "null") {
        x_.push_back(0);
        x_.push_back(0);
        x_.push_back(0);
        x_.push_back(0);
    }
    else if(init == "identity") {
        x_.push_back(1);
        x_.push_back(0);
        x_.push_back(0);
        x_.push_back(0);
    }
    else {
        std::uniform_real_distribution<double> distribution(-1,1);
        double length;
        while(true) {
            x_.clear();
            length = 0;
            for(int i = 0; i < 4; i++) {
                x_.push_back(distribution(mcmc::util::gen));
                length += x_[i]*x_[i];
            }
            length = sqrt(length);
            if(length <= 1)
                break;
        }
        for(int i = 0; i < 4; i++) x_[i] = x_[i]/length;
    }
}

template <typename T>
SU2<T>::SU2(const SU2<T>& A, double beta) {
    // Proposal according to heatbath!!
    *this = A;
    double determinante = this->det();
    if(determinante == 0)
        *this = SU2<T>();
    else {
        double a = sqrt(determinante);
        *this /= a; // = V
        double x0 = compute_x0(beta,a);
        double length_x = sqrt(1-x0*x0);
        std::uniform_real_distribution<double> distribution(-1,1);
        double x1,x2,x3;
        double length_r;
        while(true) {
            x1 = distribution(mcmc::util::gen);
            x2 = distribution(mcmc::util::gen);
            x3 = distribution(mcmc::util::gen);
            length_r = sqrt(x1*x1+x2*x2+x3*x3);
            if(length_r <= 1)
                break;
        }
        SU2<T> X(x0,x1*length_x/length_r,x2*length_x/length_r,x3*length_x/length_r);
        *this = X*(this->adjungate());
    }
}

template <typename T>
double SU2<T>::trace() {
    return x_[0]*2.0;
}

template <typename T>
double SU2<T>::det() {
    double determinante = 0;
    for(auto i = 0; i < 4; i++) determinante += x_[i]*x_[i];
    return determinante;
}

template <typename T>
SU2<T>& SU2<T>::adjungate() {
    x_[3] = -x_[3];
    x_[2] = -x_[2];
    x_[1] = -x_[1];
    return *this;
}


template<typename T>
SU2<T> operator*(const SU2<T>& x, const double& y)
{
    SU2<T>x_(x);
    for(auto i = 0; i < 4; i++)
        x_(i) *= y;
    return x_;
}

template <typename T>
SU2<T> operator-(const SU2<T>& a,const SU2<T>& b)
{
    SU2<T> temp(a);
    temp -= b;
    return temp;
}

template <typename T, typename T2>
SU2<T> operator/(const SU2<T>& x, const T2& y)
{
    SU2<T> temp(x);
    temp /= y;
    return temp;
}

template <typename T>
SU2<T> operator*(const SU2<T>& x, const SU2<T>& y)
{
    return SU2<T>(x(0)*y(0)-x(1)*y(1)-x(2)*y(2)-x(3)*y(3),y(0)*x(1)+x(0)*y(1)-x(2)*y(3)+x(3)*y(2),y(0)*x(2)+x(0)*y(2)-x(3)*y(1)+x(1)*y(3),y(0)*x(3)+x(0)*y(3)-x(1)*y(2)+x(2)*y(1));
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const SU2<T>& x) {
    os << "(("<<x(0)<<"+I*"<<x(3)<<","<<x(2)<<"+I*"<<x(1)<<"),("<<-x(2)<<"+I*"<<x(1)<<","<<x(0)<<"+I*"<<-x(3)<<"))";
    return os;
}


namespace std {
    template<typename T>
    std::string to_string(SU2<T> x)
    {
        std::string conf = "";
        for(auto j = 0; j < 4; j++)
            conf += std::to_string(x(j)) + " ";
        conf = conf.substr(0, conf.size() -1);
        return conf;
    }

    template<typename T>
    double fabs(SU2<T> x)
    {
        // ToDo!! Incorrect return type for SU2
        return std::fabs(x(0));
    }
}

#endif //MAIN_SU2_HPP
