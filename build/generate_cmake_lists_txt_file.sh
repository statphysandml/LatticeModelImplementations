target_link_libraries_appendix="\${MCMCSimulationLib} \${ParamHelper}"
cat >../CMakeLists.txt <<EOL
cmake_minimum_required(VERSION 3.13)
project(LatticeModelImplementations)

set(CMAKE_CXX_STANDARD 14)

# Boost
EOL
if [ -v path_to_boost ]; then
target_link_libraries_appendix="${target_link_libraries_appendix} \${Boost_LIBRARIES}"
cat >>../CMakeLists.txt <<EOL
set(BOOST_ROOT "${path_to_boost}")
FIND_PACKAGE( Boost REQUIRED COMPONENTS filesystem)
EOL
else
target_link_libraries_appendix="${target_link_libraries_appendix} \${Boost_LIBRARIES}"
cat >>../CMakeLists.txt <<EOL
FIND_PACKAGE( Boost 1.67 REQUIRED COMPONENTS filesystem)
if(Boost_FOUND)
    include_directories(\${Boost_INCLUDE_DIRS})
    message("Boost = \${Boost_INCLUDE_DIRS}")
endif()
EOL
fi
cat >>../CMakeLists.txt <<EOL

# Ceres -> not needed at the moment
# find_package(Ceres REQUIRED)
# include_directories(\${CERES_INCLUDE_DIRS})

EOL
if [ -v path_to_python3 ]; then
target_link_libraries_appendix="${target_link_libraries_appendix} \${PYTHON_LIBRARIES}"
cat >>../CMakeLists.txt <<EOL
# Python
set(PYTHON_LIBRARIES "${path_to_python3}lib/libpython3.7m.so")
set(PYTHON_EXECUTABLE "${path_to_python3}bin/python3.7m")
set(Python3_ROOT_DIR "${path_to_python3}")
set(PYTHON_INCLUDE_DIRS "${path_to_python3}include/python3.7m")
include_directories(\${PYTHON_INCLUDE_DIRS})
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
message("Python executable = \${PYTHON_EXECUTABLE}")

option(PYTHON "Enable Python" ON)

EOL
else
cat >>../CMakeLists.txt <<EOL
option(PYTHON "Disable Python" OFF)
EOL
fi
cat >>../CMakeLists.txt <<EOL

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

    set(CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++ -lboost_system -lboost_filesystem")

    if(CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS_DEBUG "\${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -Werror")
      set(CMAKE_CXX_FLAGS_RELEASE "\${CMAKE_CXX_FLAGS_RELEASE} -O3")
      set(CMAKE_EXE_LINKER_FLAGS "-s")  # Strip binary
    endif()
        
    cuda_add_library(latticemodelimplementations STATIC src/main.cpp)
    set_property(TARGET latticemodelimplementations PROPERTY CUDA_STANDARD 14)

    if (PYTHON)
      target_compile_definitions(latticemodelimplementations PUBLIC -D PYTHON)
    endif()

    if (THRUST)
      target_compile_definitions(latticemodelimplementations PUBLIC -D THRUST)
    endif()

    if (GPU)
      target_compile_definitions(latticemodelimplementations PUBLIC -D GPU)
    endif()

    target_link_libraries(latticemodelimplementations ${target_link_libraries_appendix})
    set_target_properties(latticemodelimplementations PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
else()
    set(CMAKE_CXX_FLAGS "\${CMAKE_CXX_FLAGS} -std=c++14 -static-libstdc++ -lboost_system -lboost_filesystem")
    
    if(CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -Werror")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    endif()

    if( CMAKE_BUILD_TYPE MATCHES "RELEASE" )
        set(CMAKE_EXE_LINKER_FLAGS "-s")  # Strip binary
    endif()
    
    add_library(latticemodelimplementations STATIC src/main.cpp)

    if (PYTHON)
      target_compile_definitions(latticemodelimplementations PUBLIC -D PYTHON)
    endif()

    target_link_libraries( latticemodelimplementations ${target_link_libraries_appendix})
endif()

SET( APP_EXE StaticTest )

ADD_EXECUTABLE( \${APP_EXE}
        src/main.cpp )
        
TARGET_LINK_LIBRARIES( \${APP_EXE}
        latticemodelimplementations)

EOL
