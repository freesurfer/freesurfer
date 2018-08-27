foreach(filename ${SOURCES})
  LIST(APPEND GEMS_SOURCES "../${filename}")
endforeach()

pybind11_add_module(gems
  python.cxx
  pyKvlCalculator.cxx
  pyKvlImage.cxx
  pyKvlMesh.cxx
  pyKvlOptimizer.cxx
  pyKvlTransform.cxx
  ${GEMS_SOURCES}
)

target_link_libraries(gems PRIVATE ${ITK_LIBRARIES} ${ZLIB_LIBRARIES})
set_target_properties(gems PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/freesurfer/gems)
