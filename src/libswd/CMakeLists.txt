#file(GLOB SOURCE "*.f")
#list(FILTER SOURCE EXCLUDE REGEX "main.cpp")

#message("myfile = ${SOURCE}")
pybind11_add_module(libswd main.cpp utils.cpp)
target_link_libraries(libswd PRIVATE swd)
install(TARGETS libswd
  COMPONENT python
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
)
# add_library(libcps330 ${SOURCE})