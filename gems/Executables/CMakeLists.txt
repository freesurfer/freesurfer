option(GEMS_MAKE_SPARSE_INITIAL_MESHES "Make sparse initial meshes" ON)

include_directories(${FS_INCLUDE_DIRS}) # added for linking to FS

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
add_executable(kvlBuildAtlasMesh.tetgen
  kvlBuildAtlasMesh.cxx
  kvlAtlasMeshBuilder.cxx
  kvlAtlasMeshCollectionFastReferencePositionCost.cxx
  kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx
  kvlMultiResolutionAtlasMesher.cxx
)

add_executable(samseg-atlas
  samseg-atlas.cxx
  kvlAtlasMeshBuilder.cxx
  kvlAtlasMeshCollectionFastReferencePositionCost.cxx
  kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx
  kvlMultiResolutionAtlasMesher.cxx
)
add_executable(samseg-atlas.tetgen
  samseg-atlas.cxx
  kvlAtlasMeshBuilder.cxx
  kvlAtlasMeshCollectionFastReferencePositionCost.cxx
  kvlAtlasMeshCollectionModelLikelihoodCalculator.cxx
  kvlMultiResolutionAtlasMesher.cxx
)
add_help(samseg-atlas samseg-atlas.help.xml)

add_definitions(-DUSE_DYNAMIC_MESH)
target_link_libraries(kvlBuildAtlasMesh kvlGEMSCommon_dynmesh utils)
target_link_libraries(kvlBuildAtlasMesh.tetgen kvlGEMSCommon_dynmesh utils)
target_link_libraries(samseg-atlas kvlGEMSCommon_dynmesh utils)
target_link_libraries(samseg-atlas.tetgen kvlGEMSCommon_dynmesh utils)
if(MAKE_SPARSE_INITIAL_MESHES)
  target_link_libraries(kvlBuildAtlasMesh.tetgen samseg-atlas.tetgen ${Tetgen_LIBRARIES})
endif()

