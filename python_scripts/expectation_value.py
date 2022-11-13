import os

from mcmctools.modes.expectation_value import expectation_value, load_expectation_value_results

if __name__ == '__main__':
    # To ensure to run code from this file
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # os.chdir("../ONModel/")

    # rel_data_dir = "./data/ONModelMetropolis/"
    # rel_results_dir = "./results/ONModelMetropolis/"

    # import numpy as np
    # expectation_values = expectation_value(
    #     number_of_measurements=1000, measures=["Mean", "AbsMean"],
    #     error_type="statistical",
    #     running_parameter="kappa", rp_values=list(np.linspace(-1.0, 1.0, 21)),
    #     rel_data_dir=rel_data_dir,
    #     rel_results_dir=rel_results_dir,
    #     sim_base_dir="./"
    # )

    os.chdir("../PolynomialModel/")

    rel_data_dir = "./data/PolynomialModelLangevin/"
    rel_results_dir = "./results/PolynomialModelLangevin/"

    import numpy as np
    expectation_values = expectation_value(
        number_of_measurements=1000000, measures=["Mean", "AbsMean"],
        error_type="statistical",
#        running_parameter="kappa", rp_values=list(np.linspace(-1.0, 1.0, 21)),
        rel_data_dir=rel_data_dir,
        rel_results_dir=rel_results_dir,
        sim_base_dir="./"
    )