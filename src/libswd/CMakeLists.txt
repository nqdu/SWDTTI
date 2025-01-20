#file(GLOB SOURCE "*.f")
#list(FILTER SOURCE EXCLUDE REGEX "main.cpp")

#message("myfile = ${SOURCE}")
pybind11_add_module(libswd main.cpp utils.cpp)
target_link_libraries(libswd PRIVATE libvti libtti shared )
set_target_properties(libvti libtti shared  PROPERTIES POSITION_INDEPENDENT_CODE ON)
#target_include_directories(libswd  PUBLIC ${PROJECT_SOURCE_DIR}/src)
install(TARGETS libswd
  COMPONENT python
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
)
# add_library(libcps330 ${SOURCE})