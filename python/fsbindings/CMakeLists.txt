include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${FS_INCLUDE_DIRS})

# temporary fix so that -g doesn't produce linker errors when binding cxx/python code
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

pybind11_add_module(fsbindings
  module.cpp
  bindings_numpy.cpp
  bindings_surface.cpp
  bindings_mri.cpp
  bindings_transform.cpp 
)

# link utilities
target_link_libraries(fsbindings PRIVATE utils nifti ${ITK_LIBRARIES})

# make sure the bindings library gets built into the repository even in out-of-source builds
set_target_properties(fsbindings PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/fsbindings)
