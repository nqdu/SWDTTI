#file(GLOB SOURCE "*.f")
#list(FILTER SOURCE EXCLUDE REGEX "main.cpp")

#message("myfile = ${SOURCE}")
pybind11_add_module(libswd main.cpp utils.cpp)
target_link_libraries(libswd PRIVATE swd)
set_target_properties(swd PROPERTIES POSITION_INDEPENDENT_CODE ON)
install(TARGETS libswd
  COMPONENT python
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
)
# add_library(libcps330 ${SOURCE})