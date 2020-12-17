

if __name__ == '__main__':

    ''' Normal DataLoader and ConfigDataGenerator with an OnTheFlyDataset '''

    from mcmctools.pytorch.data_generation.configdatagenerator import ConfigDataGenerator

    data_generator_args = {
        "path": "../data/XYModelHMC/",
        "total_number_of_data_per_file": 1000,
        "identifier": "expectation_value",
        "running_parameter": "beta",
        "chunksize": 400
    }

    data_generator = ConfigDataGenerator(
        # ConfigDataGenerator Args
        data_type="target_param",
        # DataGeneratoreBaseClass Args
        seed=None,
        set_seed=True,
        # ConfigurationLoader Args
        **data_generator_args)

    # Possible usage:
    # sample = data_generator.sampler()

    from pystatplottools.ppd_pytorch_data_generation.pytorch_data_utils.datasets import OnTheFlyDataset
    dataset = OnTheFlyDataset(datagenerator=data_generator, n=len(data_generator))

    from pystatplottools.ppd_pytorch_data_generation.pytorch_data_utils.dataloaders import DataLoader
    data_loader_params = {
        'batch_size': 89,
        'shuffle': True,  # Where used in dataset -> not at all, but by the Dataloader??
        'num_workers': 0}
    data_loader = DataLoader(dataset=dataset, **data_loader_params)

    # Load training data
    for batch_idx, batch in enumerate(data_loader):
        data, target = batch
        print(batch_idx, len(data))

    # Load training data - Second epoch
    for batch_idx, batch in enumerate(data_loader):
        data, target = batch
        print(batch_idx, len(data))


    ''' BatchDataLoader and BatchConfigDataGenerator with an OnTheFlyDataset '''

    # The difference to the ConfigDataGenerator is that batches are generated already by the dataloader
    # This leads to a performance boost since the underlying data frame is accessed via slicing
    # i.e: batch = self.data.iloc[i:i+batch_size].values instead of
    #    : batch = np.stack([self.data.iloc[j] for j in range(i, i+batch_size])

    from mcmctools.pytorch.data_generation.batchconfigdatagenerator import BatchConfigDataGenerator

    data_generator_args = {
        "path": "../data/XYModelHMC/",
        "total_number_of_data_per_file": 1000,
        "identifier": "expectation_value",
        "running_parameter": "beta",
        "chunksize": 400  # If no chunksize is given, all data is loaded at once
    }


    # It is also possible to use generate_data_loader here
    data_generator = BatchConfigDataGenerator(
        # BatchConfigDataGenerator Args
        batch_size=89,
        # ConfigDataGenerator Args
        data_type="target_param",
        # DataGeneratoreBaseClass Args
        seed=None,
        set_seed=True,
        # ConfigurationLoader Args
        **data_generator_args)

    # Possible usage:
    # sample = data_generator.sampler() # return batchs

    from pystatplottools.ppd_pytorch_data_generation.pytorch_data_utils.datasets import OnTheFlyDataset

    dataset = OnTheFlyDataset(datagenerator=data_generator,
                              n=len(data_generator))

    from pystatplottools.ppd_pytorch_data_generation.pytorch_data_utils.dataloaders import BatchDataLoader

    data_loader_params = {
        # 'status': "raw",  # Returns the sample from the BatchDataGenerator as it is. Can be used if the BatchDataGenerator produces already the correct datatype for training
        'shuffle': True,  # Where used in dataset -> not at all, but by the Dataloader??
        'num_workers': 0}
    data_loader = BatchDataLoader(dataset=dataset, **data_loader_params)

    # Load training data
    for batch_idx, batch in enumerate(data_loader):
        data, target = batch
        print(batch_idx, len(data))

    # Load training data - Second epoch
    for batch_idx, batch in enumerate(data_loader):
        data, target = batch
        print(batch_idx, len(data))


    ''' Same as above with using the load_on_the_fly_data_laoder function '''

    data_generator_args = {
        # ConfigDataGenerator Args
        "data_type": "target_param",
        # Args for ConfigurationLoader
        "path": "../data/XYModelHMC/",
        "total_number_of_data_per_file": 1000,
        "identifier": "expectation_value",
        "running_parameter": "beta",
        "chunksize": 400  # If no chunksize is given, all data is loaded at once
    }

    data_loader_params = {
        # 'status': "raw",  # Returns the sample from the BatchDataGenerator as it is. Can be used if the BatchDataGenerator produces already the correct datatype for training
        'shuffle': True,  # Used correctly by the Dataloader??
        'num_workers': 0}

    from pystatplottools.ppd_pytorch_data_generation.data_generation.datagenerationroutines import load_on_the_fly_data_loader
    from mcmctools.pytorch.data_generation.datagenerationroutines import data_generator_factory
    data_loader = load_on_the_fly_data_loader(
        batch_size=89,
        data_generator_args=data_generator_args,
        data_generator_name="BatchConfigDataGenerator",
        data_generator_factory=data_generator_factory,
        data_loader_params=data_loader_params,
        data_loader_name="BatchDataLoader",
    )

    # Load training data
    for batch_idx, batch in enumerate(data_loader):
        data, target = batch
        print(batch_idx, len(data))

    # Load training data - Second epoch
    for batch_idx, batch in enumerate(data_loader):
        data, target = batch
        print(batch_idx, len(data))

    # import matplotlib.pyplot as plt
    #
    # fig, axes = plt.subplots(3, 7, figsize=(14, 7))
    # ConfigDataGenerator.inspect_magnetization(axes=axes, data_loader=data_loader)
    # plt.show()
    #
    # import matplotlib.pyplot as plt
    #
    # fig, axes = plt.subplots(2, 2)
    #
    # ConfigDataGenerator.inspect_observables(axes=axes, data_loader=data_loader)
    #
    # plt.show()
