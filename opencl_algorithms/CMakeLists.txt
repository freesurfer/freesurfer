project(opencl_algorithms)

set(SOURCES oclDijkstraKernel.cpp)

if(OpenCL_FOUND)
  set(SOURCES ${SOURCES} oclCommon.cpp)
endif()

add_library(opencl_algorithms STATIC ${SOURCES})

if(OpenCL_FOUND)
  target_link_libraries(opencl_algorithms ${OpenCL_LIBRARY})
endif()
