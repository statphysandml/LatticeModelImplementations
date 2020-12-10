//
// Created by lukas on 05.10.20.
//

#ifndef COMPLEXMONTECARLO_EXECUTION_TEMPLATES_HPP
#define COMPLEXMONTECARLO_EXECUTION_TEMPLATES_HPP

#include "execution/executer.hpp"

namespace site_execution_templates
{
    // measure_interval=0 <-> run based on correlation_time_results_path
    template<typename ModelParameters, typename Algorithm>
    void run_expectation_value(const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
                               mcmc::execution::Executer::RunningMode running_mode=mcmc::execution::Executer::local, bool execute_code=true, const std::vector<std::string> additional_args={},
                               uint measure_interval=0, uint number_of_measurements=100000, uint start_measuring=10000,
                      json measures={"Mean", "ComplexConfig", "AbsoluteDetailedBalanceAccuracy", "DetailedBalanceAccuracy",
                                     "RealStepSize", "ImagStepSize", "Energy", "EnergyImag", "Drift", "DriftImag"},
                               json post_measures={"2ndMoment"}, uint n_means_bootstrap=50,
                               const std::string correlation_time_results_path_="", const uint maximum_measure_interval=0) {

        std::string rel_config_path = "/configs/" + target_name + "/";
        std::string rel_data_path = "/data/" + target_name + "/";

        std::string correlation_time_results_path;
        if(correlation_time_results_path_ == "")
            correlation_time_results_path = "/results/" + target_name + "/";
        else
            correlation_time_results_path = correlation_time_results_path_;

        if(maximum_measure_interval != 0)
        {
            // Load correlation time and adapt it, if it is larger then maximum_measure_interval
            auto correlation_time_results = param_helper::fs::read_parameter_file(correlation_time_results_path, "correlation_time_results");
            auto correlation_time = param_helper::params::entry_by_key<uint>(correlation_time_results["CorrelationTime"], "default");
            if(correlation_time > maximum_measure_interval)
                measure_interval = maximum_measure_interval;
        }

        auto mcmc_update_parameters = Algorithm::generate_update_parameters(mcmc_update_params);

        typename Algorithm::SystemBaseParams site_parameters(
                json{
                        {"measures", measures},
                        {ModelParameters::param_file_name(), model_parameters.get_json()},
                        {Algorithm::MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()}},
                rel_config_path
        );

        if(measure_interval == 0)
        {
            mcmc::execution::ExpectationValueParameters expectation_value_parameters(correlation_time_results_path, number_of_measurements, start_measuring,
                                                                    {}, // optional additional measures
                                                                    post_measures, n_means_bootstrap);

            // Generates all simulation parameters
            mcmc::simulation::SimulationParameters<typename Algorithm::SystemBaseParams, mcmc::execution::ExpectationValueParameters>::generate_traceable_simulation(
                    site_parameters, expectation_value_parameters, rel_config_path, rel_data_path);

            mcmc::execution::execute<typename Algorithm::SystemBaseParams>(mcmc::execution::ExpectationValueParameters::name(), target_name, "/./", true, running_mode, execute_code, additional_args);
        }
        else
        {
            mcmc::execution::ExpectationValueParameters expectation_value_parameters(measure_interval, number_of_measurements, start_measuring,
                                                                    {}, // optional additional measures
                                                                    post_measures, n_means_bootstrap);
            // expectation_value_parameters.write_to_file(rel_data_path);

            mcmc::simulation::SimulationParameters<typename Algorithm::SystemBaseParams, mcmc::execution::ExpectationValueParameters>::generate_traceable_simulation(
                    site_parameters, expectation_value_parameters, rel_config_path, rel_data_path);

            mcmc::execution::execute<typename Algorithm::SystemBaseParams>(mcmc::execution::ExpectationValueParameters::name(), target_name, "/./", true, running_mode, execute_code, additional_args);
        }
    }

    template<typename ModelParameters, typename Algorithm>
    void run_correlation_time(
            const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
            mcmc::execution::Executer::RunningMode running_mode=mcmc::execution::Executer::local, bool execute_code=true, const std::vector<std::string> additional_args={},
            uint max_correlation_time=400, uint minimum_sample_size=40000, uint start_measuring=10000)
    {
        std::string rel_config_path = "/configs/" + target_name + "/";
        std::string rel_data_path = "/data/" + target_name + "/";
        std::string correlation_time_results_path = "/results/" + target_name + "/";

        auto mcmc_update_parameters = Algorithm::generate_update_parameters(mcmc_update_params);

        typename Algorithm::SystemBaseParams site_parameters(
                json{
                    {ModelParameters::param_file_name(), model_parameters.get_json()},
                    {Algorithm::MCMCUpdateParams::param_file_name(), mcmc_update_parameters.get_json()}},
                rel_config_path
        );

        mcmc::execution::CorrelationTimeParameters correlation_time_parameters(minimum_sample_size, max_correlation_time, start_measuring, {"Mean"});
        // correlation_time_parameters.write_to_file(rel_data_path);

        mcmc::simulation::SimulationParameters< typename Algorithm::SystemBaseParams , mcmc::execution::CorrelationTimeParameters >::generate_traceable_simulation(
                site_parameters, correlation_time_parameters, rel_config_path, rel_data_path);

        mcmc::execution::execute< typename Algorithm::SystemBaseParams > (mcmc::execution::CorrelationTimeParameters::name(), target_name, "/./", true, running_mode, execute_code, additional_args);
    }

    template<typename Algorithm>
    void run_plot_site_distribution(const std::string target_name,
                                    const double rmin_x=-5.0, const double rmax_x=5.0, const double rmin_y=-5.0, const double rmax_y=5.0,
                                    mcmc::execution::Executer::RunningMode running_mode=mcmc::execution::Executer::local, bool execute_code=true, const std::vector<std::string> additional_args={})
    {
        std::string rel_config_path = "/configs/" + target_name + "/";
        std::string rel_data_path = "/data/" + target_name + "/";

        mcmc::execution::PlotSiteDistributionParameters site_distribution_parameters(
                json {{"xkey", "StateReal"},
                      {"ykey", "StateImag"},
                      {"rmin_x", rmin_x},
                      {"rmax_x", rmax_x},
                      {"rmin_y", rmin_y},
                      {"rmax_y", rmax_y}});
        site_distribution_parameters.write_to_file(rel_config_path);

        mcmc::simulation::SimulationParameters< typename Algorithm::SystemBaseParams , mcmc::execution::PlotSiteDistributionParameters >::generate_simulation_from_mode(
                rel_config_path, rel_data_path, mcmc::execution::PlotSiteDistributionParameters::name());

        mcmc::execution::execute< typename Algorithm::SystemBaseParams > (mcmc::execution::PlotSiteDistributionParameters::name(), target_name, "/./", true, running_mode, execute_code, additional_args);
    }

    // The following execution templates should be executed with separate target_names
    
    template<typename ModelParameters, typename Algorithm>
    void run_correlation_time_and_expectation_value(
            const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
            uint max_correlation_time=400, uint minimum_sample_size=40000, uint number_of_measurements=100000, uint start_measuring=10000,
            json measures={"Mean", "ComplexConfig", "AbsoluteDetailedBalanceAccuracy", "DetailedBalanceAccuracy", "RealStepSize", "ImagStepSize"}, json post_measures={"2ndMoment"})
    {
        run_correlation_time<ModelParameters, Algorithm>(target_name, model_parameters, mcmc_update_params, mcmc::execution::Executer::local, true, {}, max_correlation_time, minimum_sample_size, start_measuring);
        run_expectation_value<ModelParameters, Algorithm>(target_name, model_parameters, mcmc_update_params, mcmc::execution::Executer::local, true, {}, 0, number_of_measurements, start_measuring, measures, post_measures);
    }

    template<typename ModelParameters, typename Algorithm>
    void run_expectation_value_and_plot_site_distribution(
            const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
            uint measure_interval=1, uint number_of_measurements=100000, uint start_measuring=10000,
            const double rmin_x=-2.0, const double rmax_x=2.0, const double rmin_y=-0.8, const double rmax_y=0.8)
    {
        run_expectation_value<ModelParameters, Algorithm>(target_name, model_parameters, mcmc_update_params, mcmc::execution::Executer::local, true, {},
                                                          measure_interval, number_of_measurements, start_measuring, {"Mean", "ComplexConfig"}, {});
        run_plot_site_distribution<Algorithm>(target_name, rmin_x, rmax_x, rmin_y, rmax_y, mcmc::execution::Executer::local, true, {});
    }
    
    template<typename ModelParameters, typename Algorithm>
    void run_complex_step_size_estimation(const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
            uint number_of_measurements=200000, uint start_measuring=100000, mcmc::execution::Executer::RunningMode running_mode = mcmc::execution::Executer::local, bool execute=true, std::vector<std::string> additional_args={})
    {
        run_expectation_value<ModelParameters, Algorithm>(target_name, model_parameters, mcmc_update_params, running_mode, execute, additional_args, 1, number_of_measurements, start_measuring,
                              {"ComplexStepSize", "RealStepSize", "ImagStepSize"}, {});
    }
}

#endif //COMPLEXMONTECARLO_EXECUTION_TEMPLATES_HPP
