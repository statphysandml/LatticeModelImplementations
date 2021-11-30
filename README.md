LatticeModelImplementions
=================

The repository contains example code for MCMC simulations and evaluations of different interesting model in physics. The code makes use of the C++ LatticeModelSimulationLib (https://github.com/statphysandml/LatticeModelSimulationLib). The resulting simulation data is evaluated in Python with modules of the MCMCEvaluationLib (https://github.com/statphysandml/MCMCEvaluationLib) and of the pystatplottools library  (https://github.com/statphysandml/pystatplottools). Besides the evaluation, the sampled configurations can be loaded with the help of the libraries into PyTorch. The generation of a respective data loader is shown for each of the models. The training of achine learning algorithms with Monte Carlo samples is therefore also straightforward

In the long run, the project should serve as an easy-to-use simulator for various models which is equipped in addition with powerful libraries for an convenient evaluation of the simulation data in Python. The implemented models can be simulated and evaluated right away.

The C++ code for the different implemented models can be found in the examples/ folder. A possible computation of interesting observables and figures is implemented in jupyter. The respective notebooks are stored in the jupyter_notebooks/ directory.

The code is also integrated as a submodule in the examples/ directory of the LatticeModelSimulationLib.

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

The building process is equivalent to the one described in the Examples section of the LatticeModelSimulationLib. The only difference is that the path to the lattice model simulation lib needs to be set explicitly by the CMake variable `PATH_TO_LATTICE_MODEL_SIMULATION_LIB` if the respository is not integrated as a directory in LatticeModelSimulationLib (default is "../").

After a navigation into the top-level directory, the following sequence of commands builds the respective example. Note that if you use a virtual envirnonment, it is important that it is activate for the building process to find the correct python version.

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPATH_TO_LATTICE_MODEL_SIMULATION_LIB=[path_to_the_lattice_model_simulation_lib]/LatticeModelSimulationLib .. 
cmake --build .
```

The build process can be customized with the following CMake variables,
which can be set by adding `-D<var>={ON, OFF} / {other}` to the `cmake` call:

* `CLUSTER_MODE`: Indicates whether the code is executed on the cluster (`on_cluster`) or locally (`local`) (default: `local`). The "local" option can be used to test whether the library prepares a computation on a cluster correctly. In this case, the simulation runs locally on your machine. The option can be changed to "on_cluster". In this case the jobs are sent to the cluster. There are two functions "prepare_execution_on_cpu_cluster" and "run_execution_on_cpu_cluster" that take care of this. The functions can be found in the file ext/MCMCSimulationLib/src/execution/execution.cpp and need to be adapted according to the used cluster. More details on how to execute a simulation on the cluster can be found in the main.cpp file of the SimulateAndExecute example (https://github.com/statphysandml/MCMCSimulationLib/blob/master/examples/SimulateAndExecute//src/main.cpp) or in the main.cpp file of the template project (see Template Project).
* `PYTHON_SCRIPTS_PATH`: Path to a directory including additional python files for a possible execution of code of custom written functions and modules. (default: `./python_scripts`). The path is appended by the programming to sys.path. It needs to be defined relative to the project root path.


Running
-------

The generated executables can be found after the building step in the debug/ and release/ directories of the different examples (see simulation folders). It is important that the virtual environment is activated before the simulation. The two Python packages: MCMCEvaluationLib (https://github.com/statphysandml/MCMCEvaluationLib) and pystatplottools (https://github.com/statphysandml/pystatplottools) need to be installed in the virtual environment for a successful exeuction of the code.

The different jupyter notebooks in the jupyter_notebooks/ directory only run successfully after the respective simulation in C++ has generated the data.

The main.cpp file of the project contains based on the example of a O(n) simulation additional information on the different possiblities to execute the code (locally, on a cluster, etc.). The main project is build in the same way if one generates a sample project with the help of the LatticeModelSimulationLib (https://github.com/statphysandml/LatticeModelSimulationLib) and can be used as a starting point for your own project.

Support and Development
----------------------

For bug reports/suggestions/complaints please file an issue on GitHub.

Or start a discussion on our mailing list: statphysandml@thphys.uni-heidelberg.de
