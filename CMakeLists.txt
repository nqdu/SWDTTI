#  CMakeLists.txt

cmake_minimum_required(VERSION 3.13.0)

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
set(BUILD_PYLIB             ${BUILD_PYLIB})


# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
set(${var} ${val})
endif()
endmacro()
setDefault(EIGEN_INC   "~/software/eigen-3.4.0")
setDefault(USE_LAPACK FALSE)
setDefault(BUILD_PYLIB FALSE)

# project
project(SWDTTI LANGUAGES CXX Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/build/obj)

################# cmake setup #################
set(CMAKE_CXX_FLAGS  "-g")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") 
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native -ffast-math")
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -march=native -ffixed-line-length-none")
  #set(CMAKE_CXX_FLAGS_RELEASE  "-O2")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -xHost -DEIGEN_STRONG_INLINE=inline")
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -xHost -extend-source")
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

# check build pylib
if(${BUILD_PYLIB}) 
  #find_package(pybind11 REQUIRED)
  find_package(Python COMPONENTS Interpreter Development)
  find_package(pybind11 CONFIG)
  add_subdirectory(src/cps330)
  add_subdirectory(src/libswd)
  target_include_directories(libswd  PUBLIC ${PROJECT_SOURCE_DIR}/src)
endif()

# exe
add_executable(surfvti ${PROJECT_SOURCE_DIR}/src/main_vti.cpp)
add_executable(surftti ${PROJECT_SOURCE_DIR}/src/main_tti.cpp)
target_link_libraries(surfvti swd ${LAPACK_LIBRARIES})
target_link_libraries(surftti swd ${LAPACK_LIBRARIES})