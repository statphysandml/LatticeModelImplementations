import torch


if __name__ == '__main__':
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    mcmc_model = "XYModelHMC"
    # mcmc_model = "ComplexONModelMetropolis"

    path = "../data/" + mcmc_model
    root = "../data/" + mcmc_model

    from pystatplottools.ppd_utils.utils import set_up_directories
    data_root, results_root = set_up_directories(data_dir=mcmc_model, results_dir=mcmc_model,
                                                 data_base_dir="../data/", results_base_dir="../results/")

    ''' Data generation with storage of a permament file '''

    data_generator_args = {
        # ConfigDataGenerator Args
        "data_type": "target_param",
        "complex_config": True,
        # Args for ConfigurationLoader
        "path": path,
        "total_number_of_data_per_file": 10000,
        "identifier": "expectation_value",
        "running_parameter": "beta",
        "chunksize": 400  # If no chunksize is given, all data is loaded at once
    }

    # Prepare in memory dataset
    from pystatplottools.ppd_pytorch_data_generation.data_generation.datagenerationroutines import prepare_in_memory_dataset

    prepare_in_memory_dataset(
        batch_size=89,
        root=root,
        data_loader_name="BatchDataLoader",
        data_generator_name="BatchConfigDataGenerator",
        data_generator_args=data_generator_args
    )

    # Load in memory dataset
    from pystatplottools.ppd_pytorch_data_generation.data_generation.datagenerationroutines import load_in_memory_dataset
    from mcmctools.pytorch.data_generation.datagenerationroutines import data_generator_factory
    data_loader = load_in_memory_dataset(
        batch_size=89, root=root, slices=None, shuffle=True, num_workers=0,
        rebuild=False, data_generator_name="ConfigDataGenerator", data_generator_factory=data_generator_factory)

    # # Load training data
    # for batch_idx, batch in enumerate(data_loader):
    #     data, target = batch
    #     print(batch_idx, len(data))
    #
    # # Load training data - Second epoch
    # for batch_idx, batch in enumerate(data_loader):
    #     data, target = batch
    #     print(batch_idx, len(data))

    dataset_inspector = data_loader.get_dataset_inspector()
    config_dim = (4, 4)
    import numpy as np
    ab = (0.0, 2 * np.pi)

    # config, label = dataset_inspector.sampler()
    config, label = data_loader.dataset.get_random_sample()

    from pystatplottools.ppd_visualization.sample_visualization import SampleVisualization

    SampleVisualization.im_single_sample(sample=config, config_dim=config_dim, ab=ab)

    # batch, batch_label = next(iter(data_loader))
    batch, batch_label = data_loader.dataset.get_random_batch(108)

    SampleVisualization.im_batch(batch, num_samples=36, dim=(6, 6), config_dim=config_dim, ab=ab)
    SampleVisualization.im_batch_grid(batch, config_dim=config_dim, ab=ab)