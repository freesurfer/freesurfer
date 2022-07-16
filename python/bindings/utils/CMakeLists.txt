include_directories(${CMAKE_CURRENT_SOURCE_DIR})

pybind11_add_module(fsbindings
  module.cpp
  bindings_numpy.cpp
  bindings_surface.cpp
  bindings_mri.cpp
  bindings_transform.cpp 
)

# link utilities
target_link_libraries(fsbindings PRIVATE utils)

# make sure the bindings library gets built into the repository even in out-of-source builds
set_target_properties(fsbindings PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/fsbindings)
