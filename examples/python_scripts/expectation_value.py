import os

from mcmctools.modes.expectation_value import expectation_value, load_expectation_value_results

if __name__ == '__main__':
        # os.chdir("../simulations/SU2Model/")
        # expectation_value(files_dir="SU2ModelMetropolis")
        # expectation_value = load_expectation_value_results(files_dir="SU2ModelMetropolis")
        # expectation_value = load_expectation_value_results(files_dir="SU2ModelMetropolis")

        os.chdir("../")

        expectation_value(files_dir="ComplexONModelComplexLangevin")

        # expectation_value(files_dir="PolynomialModelLangevin")
        # expectation_value(files_dir="XYModelHMC")
        # expectation_value(files_dir="ONModelMetropolis")
        #
        # expectation_value(files_dir="CubicSiteModelComplexLangevin")
        # expectation_value(files_dir="ComplexXYModelComplexLangevin")
        # expectation_value(files_dir="ComplexONModelMetropolis")

        # load_expectation_value_results(files_dir="PolynomialModelLangevin")
        # load_expectation_value_results(files_dir="XYModelHMC")
        # load_expectation_value_results(files_dir="ONModelMetropolis")
        #
        # load_expectation_value_results(files_dir="CubicSiteModelComplexLangevin")
        # load_expectation_value_results(files_dir="ComplexXYModelComplexLangevin")
        # load_expectation_value_results(files_dir="ComplexONModelMetropolis")
