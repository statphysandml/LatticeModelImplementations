//
// Created by lukas on 06.10.20.
//

#ifndef COMPLEXMONTECARLO_ALGORITHM_TEST_RUNS_HPP
#define COMPLEXMONTECARLO_ALGORITHM_TEST_RUNS_HPP

#include "algorithms.hpp"
#include "execution_templates.hpp"

namespace algorithm_test_runs
{
    using namespace memory_site_algorithms;
    using namespace site_execution_templates;

    void test_run_complex_langevin() {
        std::string model_name = "TestRunCubicModelComplexLangevinB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});

        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexLangevinDynamics<ComplexPolynomialModelParameters> >(
                model_name, model_parameters, {{"epsilon", 0.02}, {"eps", 0.02}});
    }

    void test_run_second_order_complex_langevin() {
        std::string model_name = "TestRunCubicModelSecondOrderComplexLangevinB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});

        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, SecondOrderComplexLangevinDynamics<ComplexPolynomialModelParameters> >(
                model_name, model_parameters, {{"epsilon", 0.02}, {"eps", 0.02}});
    }

    void test_run_complex_gaussian_metropolis() {
        std::string model_name = "TestRunCubicModelComplexGaussianMetropolisB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});

        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexGaussianMetropolis<ComplexPolynomialModelParameters> >(
                model_name, model_parameters, {{"eps", 0.05}});
    }

    void test_run_complex_hat_function_metropolis() {
        std::string model_name = "TestRunCubicModelComplexHatFunctionMetropolisB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});

        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexHatFunctionMetropolis<ComplexPolynomialModelParameters> >(
                model_name, model_parameters, {{"eps", 0.4}});
    }

    void test_run_complex_implicit_gaussian() {
        std::string model_name = "TestRunCubicModelComplexImplicitGaussianB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});

        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexImlicitGaussianAlgorithm<ComplexPolynomialModelParameters> >(
                model_name, model_parameters, {{"eps", 0.02}});
    }

    void test_run_complex_implicit_hat_function() {
        std::string model_name = "TestRunCubicModelComplexImplicitHatFunctionB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});

        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexImlicitHatFunctionAlgorithm<ComplexPolynomialModelParameters> >(
                model_name, model_parameters, {{"eps", 0.5}});
    }

    void test_run_complex_uniform_metropolis_with_uniform_sampler() {
        std::string model_name = "TestRunCubicModelComplexUniformMetropolisUniformSamplerB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});
        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexUniformMetropolis<ComplexPolynomialModelParameters, UniformSampler> >(
                model_name, model_parameters, {{"eps", 0.2}});
    }

    void test_run_complex_uniform_metropolis_with_gaussian_sampler() {
        std::string model_name = "TestRunCubicModelComplexUniformMetropolisGaussianSamplerB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});
        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexUniformMetropolis<ComplexPolynomialModelParameters, GaussianSampler> >(
                model_name, model_parameters, {{"eps", 0.02}});
    }

    void test_run_complex_implicit_uniform_with_uniform_sampler() {
        std::string model_name = "TestRunCubicModelComplexImplicitUniformUniformSamplerB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});
        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexImlicitUniformAlgorithm<ComplexPolynomialModelParameters, UniformSampler> >(
                model_name, model_parameters, {{"eps", 0.5}});
    }

    void test_run_complex_implicit_uniform_with_hat_function_sampler() {
        std::string model_name = "TestRunCubicModelComplexImplicitUniformHatFunctionSamplerB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});
        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexImlicitUniformAlgorithm<ComplexPolynomialModelParameters, HatFunctionSampler> >(
                model_name, model_parameters, {{"eps", 0.5}});
    }

    /* void test_run_complex_implicit_uniform_with_gaussian_sampler() {
        std::string model_name = "TestRunCubicModelComplexImplicitUniformGaussianSamplerB";
        ComplexPolynomialModelParameters model_parameters({1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0});
        run_expectation_value_and_plot_site_distribution< ComplexPolynomialModelParameters, ComplexImlicitUniformAlgorithm<ComplexPolynomialModelParameters, GaussianSampler> >(
                model_name, model_parameters, {{"eps", 0.005}});
    } */

    void run_test_algorithms()
    {
        test_run_complex_langevin();
        test_run_second_order_complex_langevin();
        test_run_complex_gaussian_metropolis();
        test_run_complex_hat_function_metropolis();
        test_run_complex_implicit_gaussian();
        test_run_complex_implicit_hat_function();
        test_run_complex_uniform_metropolis_with_uniform_sampler();
        test_run_complex_uniform_metropolis_with_gaussian_sampler();
        test_run_complex_implicit_uniform_with_uniform_sampler();
        test_run_complex_implicit_uniform_with_hat_function_sampler();
        // test_run_complex_implicit_uniform_with_gaussian_sampler();  -> Does not work
    }
}

#endif //COMPLEXMONTECARLO_ALGORITHM_TEST_RUNS_HPP
