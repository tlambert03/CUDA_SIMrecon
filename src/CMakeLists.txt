# ##############################################################################
# Build instruction for CUDA SIMrecon project
#
# 1. Use conda to set up dependencies: conda env create -f environment.yml  # for
#   linux conda env create -f environment_win.yml  # for windows conda activate
#   simbuild
#
# 1. CUDA: The linux environment comes with the cudatoolkit-dev package that
#   includes the nvcc compiler For windows, you will need to install the CUDA
#   toolkit manually
#
# 1. MRC/DV support (optional) If you want to build with MRC support add the IVE
#   libraries to the src directory with the following paths: For Linux:
#   src/IVE/linux64/INCLUDE src/IVE/linux64/LIB For Windows:
#   src/IVE/win64/INCLUDE src/IVE/win64/LIB
#
# If you DON'T want to build with MRC or do not have access to the IVE libraries
# please run cmake with the `-DBUILD_MRC=OFF` option
#
# To build the otfviewer, use `-DBUILD_OTF_VIEWER=ON`
#
# 1. run build.sh (linux) or bld.bat (windows) inspect those scripts for more
#   details
#
# 1. If building is successful, an executable cudaSireconDriver is generated in
#   cmake_build/cudaSirecon. Run the test to make sure it works: $
#   test_data/run.sh
#
# ##############################################################################

option(BUILD_MRC "Build cudasirecon with MRC support" ON)
option(BUILD_OTF_VIEWER "Build OTF viewer (requires X11)" OFF)

cmake_minimum_required(VERSION 3.14)

# This is the one place where version is set
project(cudasirecon VERSION 1.2.0)
configure_file("${PROJECT_SOURCE_DIR}/cudasireconConfig.h.in" "${PROJECT_SOURCE_DIR}/cudasireconConfig.h")
include_directories("${PROJECT_SOURCE_DIR}")

enable_language(CUDA)
enable_testing()


if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES "52;60;61;70;72;75")
endif()
find_package(CUDA REQUIRED)
find_package(CUDAToolkit REQUIRED)
message(STATUS "FOUND CUDA: " ${CUDAToolkit_BIN_DIR})
include_directories(${CUDAToolkit_INCLUDE_DIRS})

add_definitions(-Dcimg_display=0)

find_package(OpenMP)
message(STATUS "OpenMP FOUND: " ${OPENMP_FOUND})
if(OPENMP_FOUND)
  # nvcc and openmp with CImg aren't working at the moment
  # https://github.com/dtschump/CImg/issues/286
  add_definitions(-Dcimg_use_openmp=0)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(TIFF)
find_package(
  Boost
  COMPONENTS program_options filesystem system
  REQUIRED)
link_directories(${Boost_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

if(WIN32)
  if(${MSVC_VERSION} GREATER 1800)
    find_library(LEG_STDIO legacy_stdio_definitions)
    message(STATUS "adding legacy stdio definitions for >VS14.0 " ${LEG_STDIO})
  endif()
endif(WIN32)

message(STATUS "CUDA_NVCC_FLAGS: " ${CUDA_NVCC_FLAGS})


add_subdirectory(Buffers)
add_subdirectory(cudaSirecon)
add_subdirectory(otf)
