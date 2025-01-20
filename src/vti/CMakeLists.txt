file(GLOB SOURCE "*.cpp")
list(FILTER SOURCE EXCLUDE REGEX "main_vti.cpp")

add_library(libvti ${SOURCE})

# exe
add_executable(surfvti ${PROJECT_SOURCE_DIR}/src/vti/main_vti.cpp)
target_link_libraries(surfvti libvti shared ${LAPACK_LIBRARIES})