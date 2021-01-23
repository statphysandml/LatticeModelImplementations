LatticeModelImplementions
=================

The repository contains example code for an MCMC simulation and evaluation of different interesting model in physics. The code makes use of the C++ LatticeModelSimulationLib (https://github.com/statphysandml/LatticeModelSimulationLib). The resulting simulation data is evaluated in Python with modules of the MCMCEvaluationLib (https://github.com/statphysandml/MCMCEvaluationLib) and of the pystatplottools library  (https://github.com/statphysandml/pystatplottools). Besides the evaluation, the sampled configurations can be loaded with the help of the libraries into PyTorch. The generation of a respective data loader is shown for each of the models. The training of achine learning algorithms with Monte Carlo samples is therefore also straightforward

In the long run, the project should serve as an easy-to-use simulator for various models which is equipped in addition with powerful libraries for an convenient evaluation of the simulation data in Python. The implemented models can be simulated and evaluated right away.

The C++ code for the different implemented models can be found in the simulations/ folder. A possible computation of interesting observables and figures is implemented in jupyter. The respective notebooks are stored in the jupyter_notebooks/ directory.

The code of the used C++ LatticeModelSimulationLib is integrated into the project as a submodule in the external_submodules/ directory.

Implemented Models
------------------

The repository provides code for the simulation of the following models:

- the Ising model. See: https://en.wikipedia.org/wiki/Ising_model.
- the classical XY model. See: https://en.wikipedia.org/wiki/Classical_XY_model.
- the O(n) model. See chapter 3 in "Introduction to quantum fields on a lattice: A robust mate" from Jan Smit (https://inspirehep.net/literature/601059).
- the U(1) model. See chapter 4.2 in "Quantum Chromodynamics on the Lattice" from Christof Gattringer and Christian B. Lang (https://www.springer.com/de/book/9783642018497) or here: https://www.sciencedirect.com/science/article/abs/pii/0010465581901405.
- the SU(2) model. See also chapter 4.2 in "Quantum Chromodynamics on the Lattice" from Christof Gattringer and Christian B. Lang (https://www.springer.com/de/book/9783642018497) or here: https://www.sciencedirect.com/science/article/abs/pii/0010465581900862.
- the complex and real polynomial/quartic model. Introduced, for example, here: https://www.sciencedirect.com/science/article/pii/S0003491613001516?via%3Dihub.
- the complex XY model. See, for example, here: https://arxiv.org/abs/1009.5838.
- the complex ON model. An implementation of the ON model with complex coupling paramters.

All lattice models can be simulated in arbitrary dimensions.

Build
-----

The following steps are necessary for a possible execution of the code. The steps are very similar to the one explained for building the MCMCSimulationLib (see: https://github.com/statphysandml/MCMCSimulationLib) and for the LatticeModelSimulationLib (see: https://github.com/statphysandml/LatticeModelSimulationLib).

Two config files: config.sh and project_config.sh need to be generated in the bash_scripts/ directory. Templates for both files can be found there already and can be copied.

The defined parameters in the project_config.sh file are by default the correct one for the beginning. Additional information about the two parameters is given in the link to MCMCSimulationLib.

The config.sh contains information about the used virtual environment. The parameters need to be adapted for a successful execution of the data evaluation in Python.

Having these two files, all examples can be by executing the bash script build_examples.sh in the bash_scripts/ directory:

```bash
cd bash_scripts
bash build_examples.sh
```

Running
-------

The generated executables can be found after the building step in the debug/ and release/ directories of the different examples (see simulation folders). It is important that the virtual environment is activated before the simulation. The two Python packages: MCMCEvaluationLib (https://github.com/statphysandml/MCMCEvaluationLib) and pystatplottools (https://github.com/statphysandml/pystatplottools) need to be installed in the virtual environment for a successful exeuction of the code.

The different jupyter notebooks in the jupyter_notebooks/ directory only run successfully after the respective simulation in C++ has generated the data.

The main.cpp file of the project contains based on the example of a O(n) simulation additional information on the different possiblities to execute the code (locally, on a cluster, etc.). The main project is build in the same way if one generates a sample project with the help of the LatticeModelSimulationLib (https://github.com/statphysandml/LatticeModelSimulationLib) and can be used as a starting point for your own project.

Support and Development
----------------------

For bug reports/suggestions/complaints please file an issue on GitHub.

Or start a discussion on our mailing list: statphysandml@thphys.uni-heidelberg.de
