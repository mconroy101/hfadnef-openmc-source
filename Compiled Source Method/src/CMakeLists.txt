cmake_minimum_required(VERSION 3.3 FATAL_ERROR)               
project(openmc_sources CXX)
find_package(HDF5 REQUIRED COMPONENTS C HL)
add_library(source SHARED n_source.cpp)
find_package(OpenMC REQUIRED HINTS $ENV{CONDA_PREFIX}/lib64/cmake/OpenMC)               
target_link_libraries(source OpenMC::libopenmc)
target_compile_definitions(source PRIVATE 
    DATA_FILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/"
)
