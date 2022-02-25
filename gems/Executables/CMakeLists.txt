option(GEMS_MAKE_SPARSE_INITIAL_MESHES "Make sparse initial meshes" OFF)

if(MAKE_SPARSE_INITIAL_MESHES)
  message(WARNING "Since you're opting to make sparse (i.e., irregular) initial meshes, the TetGen library will be included, which has a restrictive license.")
  if(NOT Tetgen_FOUND)
    find_package(Tetgen REQUIRED)
  endif()
  add_definitions(-DUSE_TETGEN)
  include_directories(${Tetgen_INCLUDE_DIR})
endif()

add_executable(kvlBuildAtlasMesh
  kvlBuildAtlasMesh.cxx
  kvlAtlasMeshBuilder.cxx
  kvlAtlasMeshCollectionFastReferencePositionCost.cxx
  kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx
  kvlMultiResolutionAtlasMesher.cxx
)
target_link_libraries(kvlBuildAtlasMesh kvlGEMSCommon_dynmesh)
if(MAKE_SPARSE_INITIAL_MESHES)
  target_link_libraries(kvlBuildAtlasMesh ${Tetgen_LIBRARIES})
endif()
