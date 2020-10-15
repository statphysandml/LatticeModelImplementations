cat >$project_path/cmake/CMakeLists.txt <<EOL
cmake_minimum_required(VERSION 3.9)
project(${project_name})

set(CMAKE_CXX_STANDARD 14)

# Boost
EOL
if [ -v path_to_boost ]; then
cat >>$project_path/cmake/CMakeLists.txt <<EOL
set(BOOST_ROOT "${path_to_boost}")
FIND_PACKAGE( Boost REQUIRED COMPONENTS filesystem)
EOL
else
cat >>$project_path/cmake/CMakeLists.txt <<EOL
FIND_PACKAGE( Boost 1.67 REQUIRED COMPONENTS filesystem)
if(Boost_FOUND)
    include_directories(\${Boost_INCLUDE_DIRS})
    message("Boost = \${Boost_INCLUDE_DIRS}")
endif()
EOL
fi
cat >>$project_path/cmake/CMakeLists.txt <<EOL

# Ceres -> not needed at the moment
# find_package(Ceres REQUIRED)
# include_directories(\${CERES_INCLUDE_DIRS})


set(PYTHON_LIBRARIES "${path_to_python3}lib/libpython3.7m.so")
set(PYTHON_EXECUTABLE "${path_to_python3}bin/python3.7m")
set(Python3_ROOT_DIR "$path_to_python3")
include_directories("${path_to_python3}include/python3.7m")
# find_package(PythonInterp 3 REQUIRED) # Might be necessary to uncomment on your system
find_package(PythonLibs 3 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(\${PYTHON_INCLUDE_DIRS})
message("\${PYTHON_EXECUTABLE}")


find_library(ParamHelper NAMES libparamhelper.a PATHS ${path_to_param_helper}lib)
message("ParamHelper = \${ParamHelper}")
include_directories(${path_to_param_helper}include/)


find_library(MCMCSimulationLib NAMES libmcmcsimulationlib.a PATHS ${path_to_mcmc_simulation_lib}lib)
message("MCMCSimulationLib = \${MCMCSimulationLib}")
include_directories(${path_to_mcmc_simulation_lib}include/)

set(PYTHON_SCRIPTS_PATH "${path_to_mcmc_simulation_lib}python_scripts/")

if(NOT CLUSTER_MODE)
    set(CLUSTER_MODE "${cluster_mode}") # else local
endif()

if(NOT CONDA_ACTIVATE_PATH)
    set(CONDA_ACTIVATE_PATH "${path_to_conda_activate}")
endif()

if(NOT VIRTUAL_ENV)
    set(VIRTUAL_ENV "${virtual_env}")
endif()

EOL
if [ "$project_path" = "../" ]; then
cat >>$project_path/cmake/CMakeLists.txt <<EOL
configure_file(./../include/config.h.in ./../include/config.h @ONLY)
EOL
else
cat >>$project_path/cmake/CMakeLists.txt <<EOL
configure_file(./../config.h.in ./../config.h @ONLY)
EOL
fi
cat >>$project_path/cmake/CMakeLists.txt <<EOL

SET(CudaUsage "None" CACHE STRING "Some user-specified option")

if( CudaUsage MATCHES "GPU" OR CudaUsage MATCHES "CPU" )
    find_library(LatticeModelImplementations NAMES liblatticemodelimplementations.a PATHS ${path_to_lattice_model_implementations}libgpu)
    message("LatticeModelImplementations = \${LatticeModelImplementations}")
    include_directories(${path_to_lattice_model_implementations}include/)

    option( THRUST "Enable Thrust" ON)
    message("Thrust = \${THRUST}")

    # Cuda
EOL
if [ -v path_to_cuda ]; then
cat >>$project_path/cmake/CMakeLists.txt <<EOL
    set(CUDA_TOOLKIT_ROOT_DIR "${path_to_cuda}")
EOL
else
cat >>$project_path/cmake/CMakeLists.txt <<EOL
    FIND_PACKAGE(CUDA QUIET REQUIRED)
    if(CUDA_FOUND)
        message("Cuda = \${CUDA_INCLUDE_DIRS}")
    endif()
EOL
fi
cat >>$project_path/cmake/CMakeLists.txt <<EOL

    file(GLOB CUDA_FILES "src/" *.cu)

    CUDA_COMPILE(CU_O \${CUDA_FILES})


    if( CudaUsage MATCHES "GPU" )
        option( GPU "Enable GPU" ON )
        message("GPU = \${GPU}")
        list( APPEND CUDA_NVCC_FLAGS "-gencode arch=$nvcc_flag_gencode_arch, code=$nvcc_flag_gencode_code; --expt-extended-lambda; --expt-relaxed-constexpr") #  -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp"
    else()
        option( GPU "Enable GPU" OFF )
        message("GPU = \${GPU}")
        list( APPEND CUDA_NVCC_FLAGS "-Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp; --expt-extended-lambda; --expt-relaxed-constexpr") #  -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp"
    endif()

    set(CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++ -lboost_system -lboost_filesystem")

    if(CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -Werror")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    endif()

    if( CMAKE_BUILD_TYPE MATCHES "RELEASE" )
        set(CMAKE_EXE_LINKER_FLAGS "-s")  # Strip binary
    endif()

    cuda_add_executable(
            ${project_name}
EOL
if [ "$project_path" = "../" ]; then
cat >>$project_path/cmake/CMakeLists.txt <<EOL
            ./../src/main.cpp
EOL
else
cat >>$project_path/cmake/CMakeLists.txt <<EOL
            ./../main.cpp
EOL
fi
cat >>$project_path/cmake/CMakeLists.txt <<EOL
    )
    set_property(TARGET ${project_name} PROPERTY CUDA_STANDARD 14)

    target_link_libraries( ${project_name} \${LatticeModelImplementations} \${MCMCSimulationLib} \${ParamHelper} \${Boost_LIBRARIES} \${CERES_LIBRARIES} \${PYTHON_LIBRARIES})

    target_compile_definitions(${project_name} PUBLIC -D GPU -D THRUST)
    set_target_properties(${project_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
else()
    set(CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++ -lboost_system -lboost_filesystem")

    find_library(LatticeModelImplementations NAMES liblatticemodelimplementations.a PATHS ${path_to_lattice_model_implementations}lib)
    message("LatticeModelImplementations = \${LatticeModelImplementations}")
    include_directories(${path_to_lattice_model_implementations}include/)

    if(CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -Werror")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    endif()

    if( CMAKE_BUILD_TYPE MATCHES "RELEASE" )
        set(CMAKE_EXE_LINKER_FLAGS "-s")  # Strip binary
    endif()

    add_executable(
            ${project_name}
EOL
if [ "$project_path" = "../" ]; then
cat >>$project_path/cmake/CMakeLists.txt <<EOL
            ./../src/main.cpp
EOL
else
cat >>$project_path/cmake/CMakeLists.txt <<EOL
            ./../main.cpp
EOL
fi
cat >>$project_path/cmake/CMakeLists.txt <<EOL
    )

    target_link_libraries( ${project_name} \${LatticeModelImplementations} \${MCMCSimulationLib} \${ParamHelper} \${Boost_LIBRARIES}  \${CERES_LIBRARIES} \${PYTHON_LIBRARIES})
endif()

# Go to build directory and call
# cmake ../cmake/ -DCMAKE_BUILD_TYPE=Release
# make -j9
EOL
