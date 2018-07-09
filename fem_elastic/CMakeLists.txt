project(fem_elastic)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${PETSC_INCLUDE_DIR}
  ${TETGEN_INCLUDE_DIR}
  ${ITK_INCLUDE_DIRS}
)

set(COMMON3D_SOURCES
  fem_3d.cpp
  small_matrix.cpp
  morph.cpp
  morph_utils.cpp
  misc_maths.cpp
  transformUtils.cpp
  surf_utils.cpp
  ZlibStringCompressor.cpp 
)

if(Boost_FOUND)

  include_directories(${Boost_INCLUDE_DIRS})

  # createMorph
  add_executable(createMorph fcreateMorph.cpp ${COMMON3D_SOURCES})
  target_link_libraries(createMorph utils ${TETGEN_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  add_external_dependencies(createMorph tetgen)
  install(TARGETS createMorph DESTINATION bin)

  # applyMorph
  add_executable(applyMorph fapplyMorph.cpp ${COMMON3D_SOURCES})
  target_link_libraries(applyMorph utils ${TETGEN_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  add_external_dependencies(applyMorph tetgen)
  install(TARGETS applyMorph DESTINATION bin)

  # exportGcam
  add_executable(exportGcam exportGcam.cpp ${COMMON3D_SOURCES})
  target_link_libraries(exportGcam utils ${TETGEN_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  add_external_dependencies(exportGcam tetgen)
  install(TARGETS exportGcam DESTINATION bin)

endif()

# surf2vol
add_executable(surf2vol
  ${COMMON3D_SOURCES}
  fsurf2vol.cpp
  surf_powell.cpp
  surf_energy.cpp
  misc_maths.cpp
  transformUtils.cpp
  surf_utils.cpp
  ZlibStringCompressor.cpp
  pbCluster_mesh_crop.cpp
  untangler.cpp
)
target_link_libraries(surf2vol
  utils
  ${TETGEN_LIBRARY}
  ${PETSC_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${GFORTRAN_LIBRARIES}
)
add_external_dependencies(applyMorph tetgen petsc)
install(TARGETS surf2vol DESTINATION bin)

# fem_elastic library
add_library(fem_elastic STATIC fcreateMorph.cpp fapplyMorph.cpp ${COMMON3D_SOURCES})
add_external_dependencies(fem_elastic tetgen)
target_link_libraries(fem_elastic utils ${TETGEN_LIBRARY})
