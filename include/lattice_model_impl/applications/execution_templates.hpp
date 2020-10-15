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
                               Executer::RunningMode running_mode=Executer::local, bool execute_code=true, const std::vector<std::string> additional_args={},
                               uint measure_interval=0, uint number_of_measurements=100000, uint start_measuring=10000,
                      json measures={"Mean", "ComplexConfig", "AbsoluteDetailedBalanceAccuracy", "DetailedBalanceAccuracy", "RealStepSize", "ImagStepSize"}, json post_measures={"2ndMoment"}) {

        std::string rel_config_path = "/configs/" + target_name + "/";
        std::string rel_data_path = "/data/" + target_name + "/";
        std::string correlation_time_results_path = "/results/" + target_name + "/";

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
            ExpectationValueParameters expectation_value_parameters(correlation_time_results_path, number_of_measurements, start_measuring,
                                                                    {}, // optional additional measures
                                                                    post_measures);

            // Generates all simulation parameters
            SimulationParameters<typename Algorithm::SystemBaseParams, ExpectationValueParameters>::generate_traceable_simulation(
                    site_parameters, expectation_value_parameters, rel_config_path, rel_data_path);

            execute<typename Algorithm::SystemBaseParams>(ExpectationValueParameters::name(), target_name, "/./", true, running_mode, execute_code, additional_args);
        }
        else
        {
            ExpectationValueParameters expectation_value_parameters(measure_interval, number_of_measurements, start_measuring,
                                                                    {}, // optional additional measures
                                                                    post_measures);
            // expectation_value_parameters.write_to_file(rel_data_path);

            SimulationParameters<typename Algorithm::SystemBaseParams, ExpectationValueParameters>::generate_traceable_simulation(
                    site_parameters, expectation_value_parameters, rel_config_path, rel_data_path);

            execute<typename Algorithm::SystemBaseParams>(ExpectationValueParameters::name(), target_name, "/./", true, running_mode, execute_code, additional_args);
        }
    }

    template<typename ModelParameters, typename Algorithm>
    void run_correlation_time(
            const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
            Executer::RunningMode running_mode=Executer::local, bool execute_code=true, const std::vector<std::string> additional_args={},
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

        CorrelationTimeParameters correlation_time_parameters(minimum_sample_size, max_correlation_time, start_measuring, {"Mean"});
        // correlation_time_parameters.write_to_file(rel_data_path);

        SimulationParameters< typename Algorithm::SystemBaseParams , CorrelationTimeParameters >::generate_traceable_simulation(
                site_parameters, correlation_time_parameters, rel_config_path, rel_data_path);

        execute< typename Algorithm::SystemBaseParams > (CorrelationTimeParameters::name(), target_name, "/./", true, running_mode, execute_code, additional_args);
    }

    template<typename Algorithm>
    void run_plot_site_distribution(const std::string target_name,
                                    const double rmin_x=-5.0, const double rmax_x=5.0, const double rmin_y=-5.0, const double rmax_y=5.0,
                                    Executer::RunningMode running_mode=Executer::local, bool execute_code=true, const std::vector<std::string> additional_args={})
    {
        std::string rel_config_path = "/configs/" + target_name + "/";
        std::string rel_data_path = "/data/" + target_name + "/";

        PlotSiteDistributionParameters site_distribution_parameters(
                json {{"xkey", "StateReal"},
                      {"ykey", "StateImag"},
                      {"rmin_x", rmin_x},
                      {"rmax_x", rmax_x},
                      {"rmin_y", rmin_y},
                      {"rmax_y", rmax_y}});
        site_distribution_parameters.write_to_file(rel_config_path);

        SimulationParameters< typename Algorithm::SystemBaseParams , PlotSiteDistributionParameters >::generate_simulation_from_mode(
                rel_config_path, rel_data_path, PlotSiteDistributionParameters::name());

        execute< typename Algorithm::SystemBaseParams > (PlotSiteDistributionParameters::name(), target_name, "/./", true, running_mode, execute_code, additional_args);
    }

    // The following execution templates should be executed with separate target_names
    
    template<typename ModelParameters, typename Algorithm>
    void run_correlation_time_and_expectation_value(
            const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
            uint max_correlation_time=400, uint minimum_sample_size=40000, uint number_of_measurements=100000, uint start_measuring=10000,
            json measures={"Mean", "ComplexConfig", "AbsoluteDetailedBalanceAccuracy", "DetailedBalanceAccuracy", "RealStepSize", "ImagStepSize"}, json post_measures={"2ndMoment"})
    {
        run_correlation_time<ModelParameters, Algorithm>(target_name, model_parameters, mcmc_update_params, Executer::local, true, {}, max_correlation_time, minimum_sample_size, start_measuring);
        run_expectation_value<ModelParameters, Algorithm>(target_name, model_parameters, mcmc_update_params, Executer::local, true, {}, 0, number_of_measurements, start_measuring, measures, post_measures);
    }

    template<typename ModelParameters, typename Algorithm>
    void run_expectation_value_and_plot_site_distribution(
            const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
            uint measure_interval=1, uint number_of_measurements=100000, uint start_measuring=10000,
            const double rmin_x=-2.0, const double rmax_x=2.0, const double rmin_y=-0.8, const double rmax_y=0.8)
    {
        run_expectation_value<ModelParameters, Algorithm>(target_name, model_parameters, mcmc_update_params, Executer::local, true, {},
                                                          measure_interval, number_of_measurements, start_measuring, {"Mean", "ComplexConfig"}, {});
        run_plot_site_distribution<Algorithm>(target_name, rmin_x, rmax_x, rmin_y, rmax_y, Executer::local, true, {});
    }
    
    template<typename ModelParameters, typename Algorithm>
    void run_complex_step_size_estimation(const std::string target_name, ModelParameters model_parameters, const json mcmc_update_params,
            uint number_of_measurements=200000, uint start_measuring=100000, Executer::RunningMode running_mode = Executer::local, bool execute=true, std::vector<std::string> additional_args={})
    {
        run_expectation_value<ModelParameters, Algorithm>(target_name, model_parameters, mcmc_update_params, running_mode, execute, additional_args, 1, number_of_measurements, start_measuring,
                              {"ComplexStepSize", "RealStepSize", "ImagStepSize"}, {});
    }
}

#endif //COMPLEXMONTECARLO_EXECUTION_TEMPLATES_HPP
