//
// Created by lukas on 28.02.20.
//

#ifndef MAIN_SIMULATION_HEADER_HPP
#define MAIN_SIMULATION_HEADER_HPP

#include <Python.h>

#include "mcmc_simulation/header.hpp"

#include "../../include/site/site_header.hpp"
#include "../../include/site_update/site_update_header.hpp"

#include "../../include/lattice/lattice_header.hpp"
#include "../../include/lattice_update/lattice_update_header.hpp"

#include "../../include/link_lattice/link_lattice_header.hpp"

#include "../../include/local_update/local_update_header.hpp"
#include "../../include/global_update/global_update_header.hpp"

#include "../../include/lattice/lattice_models/ising_model.hpp"


typedef int BasicType;
typedef IsingModelParameters ModelParams;
typedef MetropolisUpdateParameters<ModelParams> UpdateParams;
typedef LatticeParameters< BasicType, ModelParams, UpdateParams, SequentialUpdateParameters> SystemBaseParams;

// typedef

// extern std::vector<bool> real_assign;

/* typedef std::complex<double> BasicType;
typedef CubicGaussianModelParameters<GaussianSampler> ModelParams;
typedef ComplexLangevinUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters< BasicType, ModelParams, UpdateParams > SystemBaseParams; */

/* typedef std::complex<double> BasicType;
typedef CubicGaussianModelParameters<GaussianSampler> ModelParams;
typedef SecondOrderComplexLangevinUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters< BasicType, ModelParams, UpdateParams > SystemBaseParams; */

/*typedef std::complex<double> BasicType;
typedef CubicGaussianModelParameters<GaussianSampler> ModelParams;
typedef ComplexImplicitGaussianUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters< BasicType, ModelParams, UpdateParams > SystemBaseParams; */

/* typedef std::complex<double> BasicType;
typedef CubicGaussianModelParameters<GaussianSampler> ModelParams;
typedef ComplexGaussianMetropolisUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters< BasicType, ModelParams, UpdateParams > SystemBaseParams; */

/* template<typename Derived>
struct Model<Derived>
{
    void call(filename, dir, mode_name)
    {
        if(SimulationParameters<Derived::SystemBaseParams, ExpectationValueParameters> simparams(filename, dir, mode_name))
            json_file_name == Derived::SystemBaseParams

            return true;
        else
            return false;

    }
}; */

/* template<typename Ts...>
struct Tfs
{
    using variant_type = std::variant<Ts...>;
};

Tfs<ModelA, ModelB, ModelC> a;
a::variant_type */

// struct ModelA{
/* typedef std::complex<double> BasicType;
typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParameters;
typedef ComplexLangevinUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters< BasicType, ModelParameters, UpdateParams > SystemBaseParams; */
// };

/* typedef std::complex<double> BasicType;
typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParameters;
typedef SecondOrderComplexLangevinUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters< BasicType, ModelParameters, UpdateParams > SystemBaseParams; */

/* typedef std::complex<double> BasicType;
typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParameters;
typedef ComplexImplicitGaussianUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters< BasicType, ModelParameters, UpdateParams > SystemBaseParams; */

/* typedef std::complex<double> BasicType;
typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParameters;
typedef ComplexInfinitesimalMixedGaussianUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters< BasicType, ModelParameters, UpdateParams > SystemBaseParams; */

/* typedef std::complex<double> BasicType;
typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParameters;
typedef ComplexMixedGaussianUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters< BasicType, ModelParameters, UpdateParams > SystemBaseParams; */

/* typedef std::complex<double> BasicType;
typedef ComplexPolynomialModelParameters<GaussianSampler> ModelParameters;
typedef ComplexImplicitHatFunctionUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters< BasicType, ModelParameters, UpdateParams > SystemBaseParams; */

/* typedef double T;
typedef PolynomialModelParameters<UniformSampler> ModelParameters;
typedef MetropolisUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters<T, ModelParameters, UpdateParams> SystemBaseParams; */

/* typedef double T;
typedef PolynomialModelParameters<UniformSampler> ModelParameters;
typedef LangevinUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters<T, ModelParameters, UpdateParams> SystemBaseParams; */

/* typedef double T;
typedef PolynomialModelParameters<UniformSampler> ModelParameters;
typedef SecondOrderLangevinUpdateParameters<ModelParameters> UpdateParams;
typedef SiteParameters<T, ModelParameters, UpdateParams> SystemBaseParams; */

/* typedef NVec<double, 2> T;
typedef RepaPolynomialModelParameters<GaussianSampler> ModelParams;
typedef RepaLangevinUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters<T, ModelParams, UpdateParams> SystemBaseParams; */

// typedef NVecComplex<double, 2, real_assign> T;
//typedef NVecComplex<double, 3, real_assign> T;
/* typedef NVec<double, 2> T;
typedef RepaGaussianModelParameters<GaussianSampler> ModelParams;
typedef RepaLangevinUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters<T, ModelParams, UpdateParams> SystemBaseParams; */

/* typedef std::complex<double> T;
typedef ComplexScalarGaussianModelParameters<UniformSampler> ModelParams;
typedef ComplexLangevinUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters<T, ModelParams, UpdateParams> SystemBaseParams; */

/* typedef double T;
typedef PolynomialModelParameters<UniformSampler> ModelParams;
typedef MetropolisUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters< T, ModelParams, UpdateParams > SystemBaseParams; */

/* typedef double T;
typedef PolynomialModelParameters<UniformSampler> ModelParams;
typedef LangevinUpdateParameters<ModelParams> UpdateParams;
typedef SiteParameters< T, ModelParams, UpdateParams > SystemBaseParams; */

#endif //MAIN_SIMULATION_HEADER_HPP
