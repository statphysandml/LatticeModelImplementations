//
// Created by lukas on 23.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONS_COMPLEX_GAUSSIAN_DISTRIBUTION_FROM_FILE_HPP
#define LATTICEMODELIMPLEMENTATIONS_COMPLEX_GAUSSIAN_DISTRIBUTION_FROM_FILE_HPP

#include <fstream>
#include <complex>
#include <random>
#include <sstream>

#include "param_helper/fileos.hpp"
#include "param_helper/pathfinder.hpp"
#include "boost/filesystem.hpp"

class complex_gaussian_distribution_from_file
{
public:
    complex_gaussian_distribution_from_file(std::complex<double> a_, const int N=100000, const std::string files_dir_="ThrustGaussianDistribution") : a(a_), files_dir(files_dir_)
    {
        std::string rel_data_path = "/data/" + files_dir + "_" + std::to_string(a.real()) + "_" + std::to_string(a.imag()) + "/";
        std::string filename = "expectation_value";
        std::cout << gcp() + rel_data_path + "/" +  filename + ".dat" << std::endl;
        if(!boost::filesystem::is_regular_file(gcp() + rel_data_path + filename + ".dat"))
        {
            std::cout << "Data file not found" << std::endl;
        }

        std::ifstream infile(gcp() + rel_data_path + filename + ".dat");

        double real_part, imag_part;
        std::string e, f;
        infile >> e >> f;
        auto c = 0;
        while (infile >> real_part >> imag_part)
        {
            real_random_numbers.push_back(real_part);
            imag_random_numbers.push_back(imag_part);
            c++;
            if(c % int(N * 1.0/10) == 0)
                std::cout << "Loaded " << c * 1.0 / N * 100 << " % of the random numbers" << std::endl;
            if(c == N)
                break;
        }

        real_it = real_random_numbers.begin();
        imag_it = imag_random_numbers.begin();
    }

    template< class Generator >
    double operator()( Generator& g )
    {
        if(real_it == real_random_numbers.end())
        {
            std::cout << "Out of random numbers" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        double val = *real_it;
        real_it++;
        return val;
    }

    std::complex<double> sample_complex()
    {
        if(real_it == real_random_numbers.end())
        {
            std::cout << "Out of random numbers" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        double real_val = *real_it;
        double imag_val = *imag_it;
        real_it++;
        imag_it++;
        return {real_val, imag_val};
    }

    static std::string name()
    {
        return "complex_gaussian_distribution_from_file";
    }

private:
    const std::complex<double> a;
    const std::string files_dir;
    std::vector<double> real_random_numbers;
    std::vector<double> imag_random_numbers;

    std::vector<double>::iterator real_it;
    std::vector<double>::iterator imag_it;
};

#endif //LATTICEMODELIMPLEMENTATIONS_COMPLEX_GAUSSIAN_DISTRIBUTION_FROM_FILE_HPP
