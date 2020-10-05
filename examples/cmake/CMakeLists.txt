cmake_minimum_required(VERSION 3.9)
project(LatticeModelImplementationsExamples)

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

set(PYTHON_LIBRARIES "/home/lukas/.miniconda3/envs/flowequation/lib/libpython3.7m.so")
set(PYTHON_EXECUTABLE "/home/lukas/.miniconda3/envs/flowequation/bin/python3.7m")
set(Python3_ROOT_DIR "/home/lukas/.miniconda3/envs/flowequation")
include_directories("/home/lukas/.miniconda3/envs/flowequation/include/python3.7m")
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(${PYTHON_INCLUDE_DIRS})
message("${PYTHON_EXECUTABLE}")


find_library(ParamHelper NAMES libparamhelper.a PATHS /home/lukas/ParamHelper/lib)
message("ParamHelper = ${ParamHelper}")
include_directories(/home/lukas/ParamHelper/include/)


find_library(MCMCSimulationLib NAMES libmcmcsimulationlib.a PATHS /home/lukas/MCMCSimulationLib/lib)
message("MCMCSimulationLib = ${MCMCSimulationLib}")
include_directories(/home/lukas/MCMCSimulationLib/include/)


SET(CudaUsage "None" CACHE STRING "Some user-specified option")

if( CudaUsage MATCHES "GPU" OR CudaUsage MATCHES "CPU" )
    find_library(LatticeModelImplementations NAMES liblatticemodelimplementations.a PATHS /home/lukas/LatticeModelImplementations/libgpu)
    message("LatticeModelImplementations = ${LatticeModelImplementations}")
    include_directories(/home/lukas/LatticeModelImplementations/include/)

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

    cuda_add_executable(
            LatticeModelImplementationsExamples
            ./../src/main.cpp
            ./../src/simulation_header.cpp
            ./../src/examples/implicit_solver/implicit_solver.cpp
            ./../src/examples/integration/integration.cu
            ./../src/examples/implicit_solver/implicit_solver.cu
    )
    set_property(TARGET LatticeModelImplementationsExamples PROPERTY CUDA_STANDARD 14)

    target_link_libraries( LatticeModelImplementationsExamples ${LatticeModelImplementations} ${MCMCSimulationLib} ${ParamHelper} ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${PYTHON_LIBRARIES})

    target_compile_definitions(LatticeModelImplementationsExamples PUBLIC -D GPU -D THRUST)
    set_target_properties(LatticeModelImplementationsExamples PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
else()
    find_library(LatticeModelImplementations NAMES liblatticemodelimplementations.a PATHS /home/lukas/LatticeModelImplementations/lib)
    message("LatticeModelImplementations = ${LatticeModelImplementations}")
    include_directories(/home/lukas/LatticeModelImplementations/include/)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++")

    add_executable(
            LatticeModelImplementationsExamples ./../src/main.cpp
            ./../src/simulation_header.cpp
            ./../src/examples/integration/integration.cpp
            ./../src/examples/implicit_solver/implicit_solver.cpp
    )

    target_link_libraries( LatticeModelImplementationsExamples ${LatticeModelImplementations} ${MCMCSimulationLib} ${ParamHelper} ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${PYTHON_LIBRARIES})
endif()

# Go to build directory and call:df
# cmake ../cmake/ -DCMAKE_BUILD_TYPE=Release
# make -j6


# if(CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -Werror")
#    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp")
#    set(CMAKE_CXX_FLAGS "-fopenmp")        ## Optimize
#    #    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
# endif()