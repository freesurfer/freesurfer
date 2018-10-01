foreach(filename ${SOURCES})
  list(APPEND GEMS_SOURCES "../${filename}")
endforeach()

# temporary fix so that -g doesn't cause linker errors
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

pybind11_add_module(gems_python
  python.cxx
  pyKvlCalculator.cxx
  pyKvlImage.cxx
  pyKvlMesh.cxx
  pyKvlOptimizer.cxx
  pyKvlTransform.cxx
  ${GEMS_SOURCES}
)

target_link_libraries(gems_python PRIVATE ${ITK_LIBRARIES} ${ZLIB_LIBRARIES})
set_target_properties(gems_python PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/freesurfer/gems)
