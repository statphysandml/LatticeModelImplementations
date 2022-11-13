# OneLinkSU3Model

The generated project serves as a template project for running your own Lattice simulation.



# Prerequisites

Building OneLinkSU3Model requires the following software installed:

* A C++14-compliant compiler
* CMake `>= 3.15`
* LatticeModelSimulationLib

# Building OneLinkSU3Model

The following sequence of commands builds OneLinkSU3Model. If you use a virtual envirnonment, it is important that it is activated for the building process to find the correct python version. The sequence assumes that your current working directory is the top-level directory
of the freshly cloned repository:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

The build process can be again customized with the following CMake variables,
which can be set by adding `-D<var>={ON, OFF} / {other}` to the `cmake` call:

* `CLUSTER_MODE`: Indicates whether the code is executed on the cluster (`on_cluster`) or locally (`local`) (default: `local`). The "local" option can be used to test whether the library prepares a computation on a cluster correctly. In this case, the simulation runs locally on your machine. The option can be changed to "on_cluster". In this case the jobs are sent to the cluster. There are two functions "prepare_execution_on_cpu_cluster" and "run_execution_on_cpu_cluster" that take care of this. The functions can be found in the file src/execution/execution.cpp and need to be adapted according to the used cluster. More details on how to execute a simulation on the cluster can be found in the main.cpp file of the SimulateAndExecute example (https://github.com/statphysandml/MCMCSimulationLib/blob/master/examples/SimulateAndExecute//src/main.cpp) or in the main.cpp file of a template project (see Template Project).
* `PYTHON_SCRIPTS_PATH`: Path to a directory including additional python files for a possible execution of code of custom written functions and modules. (default: `./python_scripts`). The path is appended by the programming to sys.path. It needs to be defined relative to the project root path.