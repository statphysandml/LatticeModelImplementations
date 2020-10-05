//
// Created by lukas on 30.09.20.
//

#ifndef LATTICEMODELIMPLEMENTATIONSEXAMPLES_IMAGINARY_GAUSSIAN_DISTRIBUTION_HPP
#define LATTICEMODELIMPLEMENTATIONSEXAMPLES_IMAGINARY_GAUSSIAN_DISTRIBUTION_HPP

// ToDo

#ifdef THRUST
// print_system_info();

// 1.0 + 1.0 I -> 0.5 - 0.5 I / 0. - 1.5 I ### 10000 10000 0.0005
// 1.0 + 0.2 I -> 0.961538 - 0.192308 I / 2.66272 - 1.10947 I ### 100000, 200000, 0.00005
// 0.2 + 1.0 I -> 0.192308 - 0.961538 I / -2.66272 - 1.10947 I
// 0.2 + 0.2 I -> 2.5 - 2.5 I / 0. - 37.5 I ### 1000000, 1000000, 0.00002
// 0.1 + 1.0 I -> 50000, 10000000, 0.000002
// 0.05 + 1.0 I -> 0.049875 - 0.997506 I / -2.97759 -0.298506 I
// 0.02 + 1.0 I -> 0.0199920 - 0.999600 I / -2.99640 -0.119904 I
/* auto M = 1000000; // 50000, 20000000, 0.000001
thrust_complex_gaussian_distribution complex_gaussian_distribution({0.1, 1.0}, 50000, 1000000, 0.00002, M, "ThrustGaussianDistribution"); // 10000000

for (auto i = 0; i < 100000 - 1 ; i++) {
    // std::cout << thrust_complex_gaussian_distribution(gen) << std::endl;
    complex_gaussian_distribution(gen);
    if (i % 10000 == 0)
        std::cout << i << "\t" << complex_gaussian_distribution.get_total_updates() << std::endl;
} */

/* complex_gaussian_distribution_from_file complex_gaussian_distribution_from_file({1.0, 0.2}, 100000, "ThrustGaussianDistribution");

for (auto i = 0; i < 100000 - 1 ; i++) {
    // std::cout << thrust_complex_gaussian_distribution(gen) << std::endl;
    complex_gaussian_distribution_from_file(gen);
    if (i % 10000 == 0)
        std::cout << i << "\t" << std::endl;
} */

/* std::normal_distribution<double> norm(0.0, 1.0);
for (auto i = 0; i < M; i++) {
    norm(gen);
    if (i % 10000 == 0)
        std::cout << i << "\t" << std::endl;
} */

#endif

#endif //LATTICEMODELIMPLEMENTATIONSEXAMPLES_IMAGINARY_GAUSSIAN_DISTRIBUTION_HPP
