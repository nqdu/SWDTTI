file(GLOB SOURCE "*.cpp" "*.f90")
list(FILTER SOURCE EXCLUDE REGEX "main_tti.cpp")

add_library(libtti ${SOURCE})


add_executable(surftti ${PROJECT_SOURCE_DIR}/src/tti/main_tti.cpp)
target_link_libraries(surftti libtti shared ${LAPACK_LIBRARIES})