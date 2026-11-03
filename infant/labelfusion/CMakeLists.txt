project(labelfusion)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

pybind11_add_module(labelfusion
  numpy.cpp
  python.cpp
  cmf3d_ml.cpp
  fib.cpp
  perform_front_propagation_3d.cpp
)

if(HOST_OS MATCHES "Rocky9")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized")
endif()

# make sure the labelfusion library gets built into the repository even in out-of-source builds
set_target_properties(labelfusion PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/fsbindings)
