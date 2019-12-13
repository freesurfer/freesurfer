include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${FS_INCLUDE_DIRS})

# temporary fix so that -g doesn't produce linker errors
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

pybind11_add_module(bindings
  python.cpp
  numpy.cpp
  volume.cpp
  surface.cpp
  xform.cpp
)

# link utilities
target_link_libraries(bindings PRIVATE utils)

# make sure the bindings library gets built into the repository even in out-of-source builds
set_target_properties(bindings PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/freesurfer)
