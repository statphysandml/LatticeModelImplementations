cmake_minimum_required(VERSION 3.9)
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


option( GPU "Enable GPU" OFF )

if( GPU )
    # Cuda
    # set(CUDA_TOOLKIT_ROOT_DIR "/opt/cuda-10.1")
    FIND_PACKAGE(CUDA QUIET REQUIRED)
    if(CUDA_FOUND)
        message("Cuda = ${CUDA_INCLUDE_DIRS}")
    endif()

    file(GLOB CUDA_FILES "src/" *.cu)
    list( APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_75,code=sm_75; --expt-extended-lambda; --expt-relaxed-constexpr") #  -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp"
    # list( APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_60,code=sm_60; --expt-extended-lambda; --expt-relaxed-constexpr") #  -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp"
    CUDA_COMPILE(CU_O ${CUDA_FILES})

    # -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp

    # ToDo: Add -O3 for release!!

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++")

    # set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})

    # SET(CUDA_SEPARABLE_COMPILATION ON)

    cuda_add_executable(
            LatticeModelImplementations
#             ./../src/thrust/thrust_integration.cu
#             ./../src/thrust/thrust_finite_intergration.cu
            ./../../src/main.cpp
            ./../src/simulation_header.cpp
    )
    set_property(TARGET LatticeModelImplementations PROPERTY CUDA_STANDARD 14)

    target_link_libraries( LatticeModelImplementations ${MCMCSimulationLib} ${ParamHelper} ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${PYTHON_LIBRARIES})

    target_compile_definitions(LatticeModelImplementations PUBLIC -D GPU)
    set_target_properties(LatticeModelImplementations PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++")

    add_executable(
            LatticeModelImplementations ./../src/main.cpp
            ./../src/simulation_header.cpp
    )

    target_link_libraries( LatticeModelImplementations ${MCMCSimulationLib} ${ParamHelper} ${Boost_LIBRARIES} ${CERES_LIBRARIES} ${PYTHON_LIBRARIES})
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