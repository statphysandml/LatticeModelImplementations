cat >../CMakeLists.txt <<EOL
cmake_minimum_required(VERSION 3.13)
project(LatticeModelImplementations)

set(CMAKE_CXX_STANDARD 14)

# Boost
EOL
if [ -v path_to_boost ]; then
cat >>../CMakeLists.txt <<EOL
set(BOOST_ROOT "${path_to_boost}")
EOL
fi
cat >>../CMakeLists.txt <<EOL
FIND_PACKAGE( Boost 1.67 REQUIRED COMPONENTS filesystem)
if(Boost_FOUND)
    include_directories(\${Boost_INCLUDE_DIRS})
    message("Boost = \${Boost_INCLUDE_DIRS}")
endif()

# Ceres -> not needed at the moment
# find_package(Ceres REQUIRED)
# include_directories(\${CERES_INCLUDE_DIRS})

# Python
set(PYTHON_LIBRARIES "${path_to_python3}lib/libpython${python_version}m.so")
set(PYTHON_EXECUTABLE "${path_to_python3}bin/python${python_version}m")
set(Python3_ROOT_DIR "$path_to_python3")
include_directories("${path_to_python3}include/python${python_version}m")
# find_package(PythonInterp 3 REQUIRED) # Not working at the university
find_package(PythonLibs 3 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(\${PYTHON_INCLUDE_DIRS})
message("Python executable = \${PYTHON_EXECUTABLE}")

find_library(ParamHelper NAMES libparamhelper.a PATHS ${path_to_param_helper}lib)
message("ParamHelper = \${ParamHelper}")
include_directories(${path_to_param_helper}include/)

find_library(MCMCSimulationLib NAMES libmcmcsimulationlib.a PATHS ${path_to_mcmc_simulation_lib}lib)
message("MCMCSimulationLib = \${MCMCSimulationLib}")
include_directories(${path_to_mcmc_simulation_lib}include/)

SET(CudaUsage "None" CACHE STRING "Some user-specified option")

if( CudaUsage MATCHES "GPU" OR CudaUsage MATCHES "CPU" )
    option( THRUST "Enable Thrust" ON)
    message("Thrust = ${THRUST}")

    # Cuda
EOL
if [ -v path_to_cuda ]; then
cat >>../CMakeLists.txt <<EOL
    set(CUDA_TOOLKIT_ROOT_DIR "${path_to_cuda}")
EOL
else
cat >>../CMakeLists.txt <<EOL
    FIND_PACKAGE(CUDA QUIET REQUIRED)
    if(CUDA_FOUND)
        message("Cuda = \${CUDA_INCLUDE_DIRS}")
    endif()
EOL
fi
cat >>../CMakeLists.txt <<EOL

    file(GLOB CUDA_FILES "src/" *.cu)

    CUDA_COMPILE(CU_O ${CUDA_FILES})

    if( CudaUsage MATCHES "GPU" )
        option( GPU "Enable GPU" ON )
        message("GPU = \${GPU}")
        list( APPEND CUDA_NVCC_FLAGS "-gencode arch=$nvcc_flag_gencode_arch,code=$nvcc_flag_gencode_code; --expt-extended-lambda; --expt-relaxed-constexpr") #  -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp"
    else()
        option( GPU "Enable GPU" OFF )
        message("GPU = \${GPU}")
        list( APPEND CUDA_NVCC_FLAGS "-Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp; --expt-extended-lambda; --expt-relaxed-constexpr") #  -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp"
    endif()

    set(CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++")

    if(CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS_DEBUG "\${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -Werror")
      set(CMAKE_CXX_FLAGS_RELEASE "\${CMAKE_CXX_FLAGS_RELEASE} -O3")
      set(CMAKE_EXE_LINKER_FLAGS "-s")  # Strip binary
    endif()

    cuda_add_library(
            latticemodelimplementations STATIC src/main.cpp
            src/lattice_model_impl/distribution/thrust_complex_gaussian_distribution.cu
            src/lattice_model_impl/thrust/thrust_header.cu
            src/lattice_model_impl/thrust/thrust_finite_integration.cu
            src/lattice_model_impl/thrust/thrust_integration.cu
    )
    set_property(TARGET latticemodelimplementations PROPERTY CUDA_STANDARD 14)

    target_link_libraries( latticemodelimplementations \${MCMCSimulationLib} \${ParamHelper} \${Boost_LIBRARIES} \${CERES_LIBRARIES} \${PYTHON_LIBRARIES} )

    target_compile_definitions(latticemodelimplementations PUBLIC -D GPU -D THRUST)
    set_target_properties(latticemodelimplementations PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
else()
    set(CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++")
    
    if(CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS_DEBUG "\${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -Werror")
      set(CMAKE_CXX_FLAGS_RELEASE "\${CMAKE_CXX_FLAGS_RELEASE} -O3")
      set(CMAKE_EXE_LINKER_FLAGS "-s")  # Strip binary
    endif()
    
    add_library(latticemodelimplementations STATIC src/main.cpp)

    target_link_libraries( latticemodelimplementations \${MCMCSimulationLib} \${ParamHelper} \${Boost_LIBRARIES} \${CERES_LIBRARIES} \${PYTHON_LIBRARIES})
endif()

EOL
