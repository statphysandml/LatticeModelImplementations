

def compute_measures(data, measure_name, sim_params):
    model_name = sim_params["systembase_params"]["model_params"]["model_params_name"]
    if model_name == "IsingModel":
        from ising_model_measures import compute_ising_model_measures
        return compute_ising_model_measures(data=data, measure_name=measure_name, sim_params=sim_params)
    elif model_name == "SUTwoModel":
        from su_two_model_measures import compute_su_two_model_measures
        return compute_su_two_model_measures(data=data, measure_name=measure_name, sim_params=sim_params)