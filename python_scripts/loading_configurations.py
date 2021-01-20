from mcmctools.loading.loading import ConfigurationLoader


if __name__ == '__main__':
    # configuration_loader_args = {
    #     "path": "../simulations/ComplexPolynomialModel/data/ComplexPolynomialModelComplexLangevin/",
    #     "total_number_of_data_per_file": 100,
    #     # "identifier": "expectation_value",
    #     # "running_parameter": "default",
    #     "chunksize": 21,
    #     "drop_last": True
    # }

    # configuration_loader_args = {
    #     "path": "../data/ComplexXYModelComplexLangevin/",
    #     "total_number_of_data_per_file": 100,
    #     "identifier": "expectation_value",
    #     "running_parameter": "mu",
    #     "chunksize": 100
    # }

    # loader = ConfigurationLoader(**configuration_loader_args)
    #
    # for i in range(6):
    #     data = loader.get_next_chunk_collection(resample=False)
    #     print(i)

    data, filenames = ConfigurationLoader.load_all_configurations(
        path="../simulations/PolynomialModel/data/PolynomialModelLangevin/",
        identifier="expectation_value",
        running_parameter="default")
    # #
    # # data, filenames = ConfigurationLoader.load_all_configurations(
    # #     path="../data/XYModelHMC/",
    # #     identifier="expectation_value",
    # #     running_parameter="default")
    # #
    # data, filenames = ConfigurationLoader.load_all_configurations(
    #     path="../simulations/ONModel/data/ONModelMetropolis/",
    #     identifier="expectation_value",
    #     running_parameter="kappa")
    #
    # data, filenames = ConfigurationLoader.load_all_configurations(
    #     path="../simulations/ComplexPolynomialModel/data/ComplexPolynomialModelComplexLangevin/",
    #     identifier="expectation_value",
    #     running_parameter="default")
    #
    # data, filenames = ConfigurationLoader.load_all_configurations(
    #     path="../simulations/ComplexXYModel/data/ComplexXYModelComplexLangevin/",
    #     identifier="expectation_value",
    #     running_parameter="mu")

    # data, filenames = ConfigurationLoader.load_all_configurations(
    #     path="../simulations/ComplexONModel/data/ComplexONModelComplexLangevin/",
    #     identifier="expectation_value",
    #     running_parameter="beta")

    # data, filenames = ConfigurationLoader.load_all_configurations(
    #     path="../simulations/SU2Model/data/SU2ModelMetropolis/",
    #     identifier="expectation_value",
    #     running_parameter="beta")

    pass