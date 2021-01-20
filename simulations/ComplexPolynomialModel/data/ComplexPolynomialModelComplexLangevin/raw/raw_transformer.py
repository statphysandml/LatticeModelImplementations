def transformer(data):
    import numpy as np
    polynomial_action = lambda x: (1.0 + 0.0j) * np.power(x, 4.0) + (1.0 + 1.0j) * np.power(x, 2.0)
    configs = data["Config"].values[:, 0] + 1.0j * data["Config"].values[:, 1]
    action = polynomial_action(configs)
    data.insert(len(data.columns), ("Action","", 0), np.real(action))
    data.insert(len(data.columns), ("Action", "", 1), np.imag(action))
    return data

