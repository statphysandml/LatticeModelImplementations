import numpy as np


from util import get_neighbour_index


def compute_on_model_measures(data, measure_name, sim_params):
    if measure_name == "TwoPointCorrelation":
        return compute_two_point_correlation(data, sim_params)


def compute_two_point_correlation(data, sim_params):
    dimensions = sim_params["systembase_params"]["dimensions"]
    dim_mul = np.cumprod([1] + dimensions)
    lattice_configs = data["Config"].values
    elem_per_site = 1
    n_sites = np.size(lattice_configs, axis=1)

    two_point_correlation = np.zeros(len(data))
    N = int(n_sites / np.prod(dimensions))

    # Using translation symmetry
    for site in range(n_sites):
        # Using rotation symmetry
        for dim in range(len(dimensions)):
            neighbour = get_neighbour_index(n=site, dim=dim, direction=True, mu=0, dimensions=dimensions, dim_mul=dim_mul,
                                            elem_per_site=elem_per_site)
            # Using symmetries in exchanging two components of the n-vector
            two_point_correlation += np.sum(lattice_configs[:, site * N:(site + 1) * N] * lattice_configs[:, neighbour * N:(neighbour + 1) * N], axis=1)

    data.insert(len(data.columns), "TwoPointCorrelation", two_point_correlation / (N * n_sites * len(dimensions)))

    return ["TwoPointCorrelation"], data

