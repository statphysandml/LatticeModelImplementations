import os


from mcmctools.loading.loading import ConfigurationLoader


if __name__ == '__main__':
    # To ensure to run code from this file
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    configuration_loader_args = {
        "path": "../simulations/ONModel/data/ONModelMetropolis/",
        "total_number_of_data_per_file": 10000,
        "identifier": "expectation_value",
        "running_parameter": "kappa",
        "chunksize": 400,
        "drop_last": True
    }

    loader = ConfigurationLoader(**configuration_loader_args)

    # Loading chunk by chunk
    for i in range(6):
        data = loader.get_next_chunk_collection(resample=False)
        print(i)

    # Loading the entire dataset
    data, filenames = ConfigurationLoader.load_all_configurations(
        path="../simulations/ONModel/data/ONModelMetropolis/",
        identifier="expectation_value",
        running_parameter="kappa")
