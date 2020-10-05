cmake_minimum_required(VERSION 3.16)
project(LatticeModelImplementations)

set(CMAKE_CXX_STANDARD 14)

# https://codeyarns.com/2013/09/13/how-to-build-cuda-programs-using-cmake/

# Uni
# set(BOOST_ROOT "/opt/boost_1_70_0")

# Boost
message("${BOOST_ROOT}")
FIND_PACKAGE( Boost 1.67 REQUIRED COMPONENTS filesystem)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message("Boost = ${Boost_INCLUDE_DIRS}")
endif()

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})


# Python

# Home
set(PYTHON_LIBRARIES "/home/lukas/.miniconda3/envs/flowequation/lib/libpython3.7m.so")
set(PYTHON_EXECUTABLE "/home/lukas/.miniconda3/envs/flowequation/bin/python3.7m")
include_directories("/home/lukas/.miniconda3/envs/flowequation/include/python3.7m")
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Uni
# include_directories("/home/kades/.conda/envs/pytorchlocal3/include/python3.7m")
# set(PYTHON_LIBRARIES "/home/kades/.conda/envs/pytorchlocal3/lib/libpython3.7m.so")
# set(PYTHON_EXECUTABLE "/home/kades/.conda/envs/pytorchlocal3/bin/python3.7m")
# find_package(PythonInterp 3 REQUIRED)

message("Python executable = ${PYTHON_EXECUTABLE}")

find_library(ParamHelper NAMES libparamhelper.a PATHS /home/lukas/ParamHelper/lib)
message("ParamHelper = ${ParamHelper}")
include_directories(/home/lukas/ParamHelper/include/)

find_library(MCMCSimulationLib NAMES libmcmcsimulationlib.a PATHS /home/lukas/MCMCSimulationLib/lib)
message("MCMCSimulationLib = ${MCMCSimulationLib}")
include_directories(/home/lukas/MCMCSimulationLib/include/)

SET(CudaUsage "None" CACHE STRING "Some user-specified option")

if( CudaUsage MATCHES "GPU" OR CudaUsage MATCHES "CPU" )
    option( THRUST "Enable Thrust" ON)
    message("Thrust = ${THRUST}")

    # Cuda
    # set(CUDA_TOOLKIT_ROOT_DIR "/opt/cuda-10.1")
    FIND_PACKAGE(CUDA QUIET REQUIRED)
    if(CUDA_FOUND)
        message("Cuda = ${CUDA_INCLUDE_DIRS}")
    endif()

    file(GLOB CUDA_FILES "src/" *.cu)

    CUDA_COMPILE(CU_O ${CUDA_FILES})


    if( CudaUsage MATCHES "GPU" )
        option( GPU "Enable GPU" ON )
        message("GPU = ${GPU}")
        list( APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_75,code=sm_75; --expt-extended-lambda; --expt-relaxed-constexpr") #  -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp"
    else()
        option( GPU "Enable GPU" OFF )
        message("GPU = ${GPU}")
        list( APPEND CUDA_NVCC_FLAGS "-Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp; --expt-extended-lambda; --expt-relaxed-constexpr") #  -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp"
    endif()
    # ToDo: Add -O3 for release!!

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++")

    cuda_add_library(
            latticemodelimplementations STATIC src/main.cpp
            src/lattice_model_impl/distribution/thrust_complex_gaussian_distribution.cu
            src/lattice_model_impl/thrust/thrust_header.cu
            src/lattice_model_impl/thrust/thrust_finite_integration.cu
            src/lattice_model_impl/thrust/thrust_integration.cu
    )
    set_property(TARGET latticemodelimplementations PROPERTY CUDA_STANDARD 14)

    target_link_libraries( latticemodelimplementations ${MCMCSimulationLib} ${ParamHelper} ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${PYTHON_LIBRARIES} )

    target_compile_definitions(latticemodelimplementations PUBLIC -D GPU -D THRUST)
    set_target_properties(latticemodelimplementations PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++")

    add_library(latticemodelimplementations STATIC src/main.cpp)

    target_link_libraries( latticemodelimplementations ${MCMCSimulationLib} ${ParamHelper} ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${PYTHON_LIBRARIES})
endif()






# CMake instructions to test using the static lib

# SET( APP_EXE StaticTest )

# ADD_EXECUTABLE( ${APP_EXE} )

# TARGET_LINK_LIBRARIES( ${APP_EXE}
#         latticemodelimplementationslib)
