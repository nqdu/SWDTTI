file(GLOB SOURCE "*.cpp" )
list(FILTER SOURCE EXCLUDE REGEX "main.cpp")

add_library(libmultiphy ${SOURCE})


add_executable(surfvti_ac ${PROJECT_SOURCE_DIR}/src/multiphysics/main.cpp)
target_link_libraries(surfvti_ac libmultiphy shared ${LAPACK_LIBRARIES})