#file(GLOB SOURCE "*.f")
#list(FILTER SOURCE EXCLUDE REGEX "main.cpp")

#message("myfile = ${SOURCE}")
pybind11_add_module(cps330 main.cpp surfdisp96.f)
install(TARGETS cps330
  COMPONENT python
  LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib
)
# add_library(libcps330 ${SOURCE})