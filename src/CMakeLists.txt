
file(GLOB SOURCE "*.cpp" "*.f90")
list(FILTER SOURCE EXCLUDE REGEX "main_*")

#message("myfile = ${SOURCE}")

add_library(swd ${SOURCE})