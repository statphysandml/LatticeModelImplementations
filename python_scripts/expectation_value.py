import os

from mcmctools.modes.expectation_value import expectation_value, load_expectation_value_results

if __name__ == '__main__':
    # To ensure to run code from this file
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    os.chdir("../simulations/ONModel/")
    expectation_value(files_dir="ONModelMetropolis")
    expectation_values = load_expectation_value_results(files_dir="ONModelMetropolis")
