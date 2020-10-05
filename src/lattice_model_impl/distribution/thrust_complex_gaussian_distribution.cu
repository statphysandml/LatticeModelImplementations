#include "../../../include/lattice_model_impl/thrust/thrust_complex_gaussian_distribution.hpp"

#ifdef THRUST

// https://docs.nvidia.com/cuda/curand/device-api-overview.html#thrust-and-curand-example

struct perform_gaussian_complex_langevin_update
{
    const double a_real;
    const double a_imag;
    const uint autocorrelation;
    const double epsilon;
    const double sqrt2epsilon;
    const unsigned long long offset;

    perform_gaussian_complex_langevin_update(const double a_real, const double a_imag, const uint autocorrelation, const double epsilon, const unsigned long long offset) :
            a_real(a_real), a_imag(a_imag), autocorrelation(autocorrelation), epsilon(epsilon), sqrt2epsilon(std::sqrt(2.0 * epsilon)), offset(offset)
    {
        std::cout << "Epsilon: " << epsilon << std::endl;
        std::cout << "Sqrt2Epsilon: " << sqrt2epsilon << std::endl;
    }

    template <typename Tuple>
    __device__
    void operator()(Tuple t) {
        unsigned int seed = thrust::get<0>(t);

        curandState s;

        // seed a random number generator
        curand_init(seed, 0, offset, &s);

        for (unsigned int i = 0; i < autocorrelation; ++i) {
            double x = thrust::get<1>(t);
            thrust::get<1>(t) = thrust::get<1>(t) - epsilon * (a_real * thrust::get<1>(t) - a_imag * thrust::get<2>(t)) + sqrt2epsilon * curand_normal(&s); // normal(rng);
            thrust::get<2>(t) = thrust::get<2>(t) - epsilon * (a_real * thrust::get<2>(t) + a_imag * x);
        }
    }
};

struct second_moment_function : public thrust::binary_function<cudaT, cudaT, float>
{
    __device__
    thrust::complex<double> operator() (cudaT &real_val, cudaT &imag_val)
    {
        return {real_val * real_val - imag_val * imag_val, 2.0 * real_val * imag_val};
    }
};

struct fourth_moment_function : public thrust::binary_function<cudaT, cudaT, float>
{
    __device__
    thrust::complex<double> operator() (cudaT &real_val, cudaT &imag_val)
    {
        return {real_val * real_val * real_val * real_val - 6.0 * real_val * real_val * imag_val * imag_val + imag_val * imag_val * imag_val * imag_val,
                4.0 * real_val * imag_val * (real_val * real_val - imag_val * imag_val)};
    }
};

struct print_functor
{
    explicit print_functor(std::ofstream & os_) : os(os_)
    {}

    template <typename Tuple>
    __host__
    void operator() (Tuple t)
    {
        os << thrust::get<0>(t) << "\t" << thrust::get<1>(t) << std::endl;
    }

    std::ofstream & os;
};

thrust_complex_gaussian_distribution::thrust_complex_gaussian_distribution(std::complex<double> a_, uint n_autocorrelation_,
                                         uint n_initialization, const double epsilon_, const int M_, const std::string files_dir_) :
            a(a_), epsilon(epsilon_), n_autocorrelation(n_autocorrelation_), M(M_), total_updates(0), files_dir(files_dir_)
    {
        std::cout << "a_real: " << a.real() << "\ta_imag: " << a.imag() << std::endl;
        // allocate storage
        real_random_numbers = dev_vec (M, 0);
        imag_random_numbers = dev_vec (M, 0);
        // indices = dev_vec_int (M, 0);
        indices = thrust::host_vector<double> (M, 0);
        host_real_random_numbers = thrust::host_vector<int> (M, 0);
        host_imag_random_numbers = thrust::host_vector<int> (M, 0);

        dist = thrust::uniform_real_distribution<float>(0.0,1.0);

        if(files_dir != "None") {
            std::string rel_data_path =
                    "/data/" + files_dir + "_" + std::to_string(a.real()) + "_" + std::to_string(a.imag()) + "/";
            std::string filename = "expectation_value";
            if(boost::filesystem::is_directory(gcp() + rel_data_path)) {
                std::cout << "Clear data directory" << std::endl;
                boost::filesystem::path path_to_remove(gcp() + rel_data_path);
                for (boost::filesystem::directory_iterator end_dir_it, it(path_to_remove); it != end_dir_it; ++it) {
                    boost::filesystem::remove_all(it->path());
                }
            }
            std::cout << "Create and prepare directory new for storage of random numbers" << std::endl;
            boost::filesystem::create_directories(gcp() + rel_data_path);
            Fileos fileos (gcp() + rel_data_path + "/" +  filename + ".dat");
            auto& os = fileos.get();
            os << "StateReal\tStateImag" << std::endl;
        }

        update_random_numbers(n_initialization);
    }

void thrust_complex_gaussian_distribution::write_data_to_file(const std::string files_dir) const
{
    std::string rel_data_path = "/data/" + files_dir + "/";

    std::string filename = "expectation_value";
    std::cout << gcp() + rel_data_path + "/" +  filename + ".dat" << std::endl;
    Fileos fileos (gcp() + rel_data_path + "/" +  filename + ".dat", true);
    auto& os = fileos.get();
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(host_real_random_numbers.begin(), host_imag_random_numbers.begin())),
                     thrust::make_zip_iterator(thrust::make_tuple(host_real_random_numbers.end(), host_imag_random_numbers.end())),
                     print_functor(os));

}

void thrust_complex_gaussian_distribution::update_random_numbers(const uint n_updates)
{
    if(n_updates > 1000000)
    {
        std::cout << "Start evolution (total updates - " << n_updates << ")" << std::endl;
        auto c = 0;
        while(c < n_updates)
        {
            thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<int>(0),
                                                                          real_random_numbers.begin(),
                                                                          imag_random_numbers.begin())),
                             thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<int>(M),
                                                                          real_random_numbers.end(),
                                                                          imag_random_numbers.end())),
                             perform_gaussian_complex_langevin_update(a.real(), a.imag(), 1000000, epsilon, total_updates));
            std::cout << "At " << 100.0 * c / n_updates << " % " << std::endl;
            total_updates += 1000000;
            c += 1000000;
        }

    }
    else
    {
        std::cout << "Start evolution" << std::endl;
        thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<int>(0),
                                                                  real_random_numbers.begin(),
                                                                  imag_random_numbers.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<int>(M),
                                                         real_random_numbers.end(),
                                                         imag_random_numbers.end())),
            perform_gaussian_complex_langevin_update(a.real(), a.imag(), n_updates, epsilon, total_updates));

        total_updates += n_updates;

    }


    std::cout << "Finished evolution, starting to sort" << std::endl;

    // Compute second moment
    thrust::device_vector<thrust::complex<double>> second_moment (M, 0);
    thrust::transform(real_random_numbers.begin(), real_random_numbers.end(), imag_random_numbers.begin(), second_moment.begin(), second_moment_function());
    thrust::complex<double> result = thrust::reduce(second_moment.begin(), second_moment.end());
    result = result / M;
    std::cout << "Second moment = " << result.real() << " + i " << result.imag() << std::endl;

    // Compute fourth moment
    thrust::device_vector<thrust::complex<double>> fourth_moment (M, 0);
    thrust::transform(real_random_numbers.begin(), real_random_numbers.end(), imag_random_numbers.begin(), fourth_moment.begin(), fourth_moment_function());
    thrust::complex<double> result_fourth_momentum = thrust::reduce(fourth_moment.begin(), fourth_moment.end());
    result_fourth_momentum = result_fourth_momentum / M;
    std::cout << "Fourth moment = " << result_fourth_momentum.real() << " + i " << result_fourth_momentum.imag() << std::endl;

    // Copy data to host vectors
    thrust::copy(real_random_numbers.begin(), real_random_numbers.end(), host_real_random_numbers.begin());
    thrust::copy(imag_random_numbers.begin(), imag_random_numbers.end(), host_imag_random_numbers.begin());

    // Write data to file
    if(files_dir != "None")
        write_data_to_file(files_dir + "_" + std::to_string(a.real()) + "_" + std::to_string(a.imag()));

    // Prepare permutation by assigning random uniform values to permutation_keys
    thrust::host_vector<double> permutation_keys(M, 0);
    thrust::sequence(indices.begin(), indices.end());
    thrust::generate(thrust::host, permutation_keys.begin(), permutation_keys.end(), [this] () { return dist(rng); });

    // Sort indices to obtain a permutation in indices
    thrust::host_vector<double> permutation_keys_dev(permutation_keys);
    thrust::sort_by_key(permutation_keys_dev.begin(), permutation_keys_dev.end(), indices.begin());

    // Generate iterator over random_numbers based on permutation in indices
    iter_real_random_numbers = thrust::permutation_iterator<thrust::host_vector<cudaT>::iterator, thrust::host_vector<int>::iterator> (
            host_real_random_numbers.begin(), indices.begin());
    iter_imag_random_numbers = thrust::permutation_iterator<thrust::host_vector<cudaT>::iterator, thrust::host_vector<int>::iterator> (
            host_imag_random_numbers.begin(), indices.begin());

    i = 0;

    /* print_range("Random Numbers", real_random_numbers.begin(), real_random_numbers.end());
    print_range("Permutation", indices.begin(), indices.end());
    print_range("Iter", iter_random_numbers, iter_random_numbers + 10); */
}

double thrust_complex_gaussian_distribution::get_random_number()
{
    if(i == M - 1) {
        update_random_numbers(n_autocorrelation);
        return iter_real_random_numbers[i];
    }
    else {
        i += 1;
        return iter_real_random_numbers[i];
    }
}

std::complex<double> thrust_complex_gaussian_distribution::get_complex_random_number()
{
    if(i == M - 1) {
        update_random_numbers(n_autocorrelation);
        return std::complex<double> {iter_real_random_numbers[i], iter_imag_random_numbers[i]};
    }
    else {
        i += 1;
        return std::complex<double> {iter_real_random_numbers[i], iter_imag_random_numbers[i]};
    }
}

#endif