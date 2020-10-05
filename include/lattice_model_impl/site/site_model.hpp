//
// Created by lukas on 11.10.19.
//

#ifndef MAIN_SITE_MODEL_HPP
#define MAIN_SITE_MODEL_HPP


#include "param_helper/params.hpp"


class SiteModelParameters : public Parameters {
public:
    using Parameters::Parameters;

    void write_to_file(const std::string& root_dir) {
        Parameters::write_to_file(root_dir, param_file_name());
    }

    static const std::string param_file_name()
    {
        return "model_params";
    }
};


template <typename ModelCL>
class SiteModel
{
public:
    template<typename T>
    T normalize(T state)
    {
        return state;
    }

    template<typename SB>
    std::vector< std::unique_ptr<common_measures::MeasurePolicy<SB>> > generate_model_measures(const json& measure_names)
    {
        return std::vector< std::unique_ptr<common_measures::MeasurePolicy<SB>> > {};
    }

    template<typename func>
    std::complex<double> get_potential(const std::complex<double> site, func& transformer) const
    {
        std::complex<double> site_{transformer(site.real()), site.imag()};
        return site_model().get_potential(site);
    }

    template<typename func>
    std::complex<double> get_drift_term(const std::complex<double> site, func& transformer) const {
        std::complex<double> site_{transformer(site.real()), site.imag()};
        return site_model()(site);
    }

private:
    ModelCL& site_model() {
        return *static_cast<ModelCL*>(this);
    }

    const ModelCL& site_model() const {
        return *static_cast<const ModelCL*>(this);
    }
};

#endif //MAIN_SITE_MODEL_HPP
