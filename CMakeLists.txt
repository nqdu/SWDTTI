#  CMakeLists.txt

cmake_minimum_required(VERSION 3.10.0)

################## compiler options ##################
# compiler
set(CMAKE_CXX_COMPILER          ${CXX})
set(CMAKE_Fortran_COMPILER          ${FC})

# compiler flags
set(CMAKE_CXX_FLAGS             ${CXXFLAGS}) 

################## dependencies ##################
# preferred installation prefix of dependencies
set(EIGEN_INC                 ${EIGEN_INC})
set(USE_LAPACK             ${USE_LAPACK}) 

# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
set(${var} ${val})
endif()
endmacro()
setDefault(EIGEN_INC   "/Users/nqdu/software/eigen-3.4.0")
setDefault(USE_LAPACK FALSE)

# project
project(SPMTOMO LANGUAGES CXX Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/build/obj)

################# cmake setup #################
set(CMAKE_CXX_FLAGS  " -g ")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") 
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native -ffast-math")
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -march=native")
  #set(CMAKE_CXX_FLAGS_RELEASE  "-O2")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -xHost -DEIGEN_STRONG_INLINE=inline")
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -xHost")
endif()

set(CMAKE_VERBOSE_MAKEFILE off)

################# find packages #################
# path of *.cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/include/cmake)

# eigen
include_directories(${EIGEN_INC})


# lapack
if(${USE_LAPACK})
  find_package(LAPACK)
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DEIGEN_USE_LAPACKE")
endif()

################# include paths #################
include_directories(${PROJECT_SOURCE_DIR}/src)

# subs
add_subdirectory(${PROJECT_SOURCE_DIR}/src/)

# exe
add_executable(surf24 ${PROJECT_SOURCE_DIR}/src/main_vti.cpp)
add_executable(surftti ${PROJECT_SOURCE_DIR}/src/main_tti.cpp)
target_link_libraries(surf24 swd ${LAPACK_LIBRARIES})
target_link_libraries(surftti swd ${LAPACK_LIBRARIES})