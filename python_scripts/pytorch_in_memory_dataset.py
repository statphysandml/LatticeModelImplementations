import os
import torch


if __name__ == '__main__':
    # To ensure to run code from this file
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    mcmc_model = "XYModelHMC"

    path = "../simulations/XYModel/data/" + mcmc_model
    root = "../simulations/XYModel/data/" + mcmc_model

    from pystatplottools.utils.utils import set_up_directories
    data_root, results_root = set_up_directories(data_dir=mcmc_model, results_dir=mcmc_model,
                                                 data_base_dir="../simulations/XYModel/data/",
                                                 results_base_dir="../simulations/XYModel/results/")

    ''' Data generation with storage of a permament file '''

    data_generator_args = {
        # ConfigDataGenerator Args
        "data_type": "target_param",
        # Args for ConfigurationLoader
        "path": path,
        "total_number_of_data_per_file": 1000,
        "identifier": "expectation_value",
        "running_parameter": "beta",
        "chunksize": 400  # If no chunksize is given, all data is loaded at once
    }

    # Prepare in memory dataset
    from pystatplottools.pytorch_data_generation.data_generation.datagenerationroutines import prepare_in_memory_dataset
    from mcmctools.pytorch.data_generation.datagenerationroutines import data_generator_factory

    prepare_in_memory_dataset(
        root=root,
        batch_size=128,
        data_generator_args=data_generator_args,
        data_generator_name="BatchConfigDataGenerator",
        data_generator_factory=data_generator_factory
    )

    # Load in memory dataset
    from pystatplottools.pytorch_data_generation.data_generation.datagenerationroutines import load_in_memory_dataset

    # The dataset is generated and stored as a .pt file in the data_dir/data directory the first time this function is called. Otherwise the .pt is loaded.
    data_loader = load_in_memory_dataset(
        root=root, batch_size=128, data_generator_factory=data_generator_factory, slices=None, shuffle=True,
        num_workers=0, rebuild=False
        # sample_data_generator_name="ConfigDataGenerator"  # optional: for a generation of new samples
    )

    # Load training data
    for batch_idx, batch in enumerate(data_loader):
        data, target = batch
        print(batch_idx, len(data))

    # Load training data - Second epoch
    for batch_idx, batch in enumerate(data_loader):
        data, target = batch
        print(batch_idx, len(data))

    dataset_inspector = data_loader.get_dataset_inspector()
    config_dim = (4, 4)
    import numpy as np
    ab = (0.0, 2 * np.pi)

    # config, label = dataset_inspector.sampler()
    config, label = data_loader.dataset.get_random_sample()

    from pystatplottools.visualization import sample_visualization

    sample_visualization.fd_im_single_sample(sample=config, config_dim=config_dim, ab=ab)

    # batch, batch_label = next(iter(data_loader))
    batch, batch_label = data_loader.dataset.get_random_batch(108)

    sample_visualization.fd_im_batch(batch, num_samples=36, dim=(6, 6), config_dim=config_dim, ab=ab)
    sample_visualization.fd_im_batch_grid(batch, config_dim=config_dim, ab=ab)