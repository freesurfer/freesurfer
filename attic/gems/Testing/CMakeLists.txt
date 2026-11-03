set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.59.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

# add an executable
add_executable(kvlAtlasMeshRasterizorTestGPU kvlAtlasMeshRasterizorTestGPU.cxx)
target_link_libraries(kvlAtlasMeshRasterizorTestGPU kvlGEMSCommon)

# testCompressionLookupTable
add_executable(testCompressionLookupTable testCompressionLookupTable.cxx)
target_link_libraries(testCompressionLookupTable kvlGEMSCommon)

add_executable(testdeformMeshPython testdeformMeshPython.cxx)
target_link_libraries(testdeformMeshPython kvlGEMSCommon)

#
add_executable(testdeformMeshPython.dynmesh testdeformMeshPython.cxx)
target_compile_definitions(testdeformMeshPython.dynmesh PRIVATE -DUSE_DYNAMIC_MESH)
target_link_libraries(testdeformMeshPython.dynmesh kvlGEMSCommon_dynmesh)

# boost tests
if(Boost_FOUND)
  set(GEMS2libs kvlGEMSCommon )
  set(testsrcs boosttests.cpp)
  list(APPEND testsrcs testatlasmeshrasterizorbasic.cpp)
  list(APPEND testsrcs atlasmeshvisitcountercpuwrapper.cpp)
  list(APPEND testsrcs atlasmeshalphadrawercpuwrapper.cpp)
  list(APPEND testsrcs testatlasmeshvisitcounter.cpp)
  list(APPEND testsrcs testatlasmeshalphadrawer.cpp)
  list(APPEND testsrcs teststopwatch.cpp)

  list(APPEND testsrcs imageutils.cpp)

  if(CUDA_FOUND)
    list(APPEND testsrcs testcudaimage.cpp)
    list(APPEND testsrcs cudaglobalfixture.cpp)
    list(APPEND testsrcs cudaimagetests.cu)
    list(APPEND testsrcs testcudatetrahedralmesh.cpp)
    list(APPEND testsrcs testdimensioncuda.cpp)
    list(APPEND GEMS2libs kvlGEMSCUDA)
    cuda_add_executable( TestGEMS2 ${testsrcs} ${cudatestsrcs})
  else()
    add_executable(TestGEMS2 ${testsrcs})
  endif()

  target_link_libraries(TestGEMS2 ${Boost_LIBRARIES} ${GEMS2libs})
endif()
               
# copy some files
configure_file(test.nii ${GEMS_RUNTIME_PATH}/test.nii COPYONLY)
configure_file(test.txt.gz ${GEMS_RUNTIME_PATH}/test.txt.gz COPYONLY)

configure_file(inp_image_deformmesh.mgz ${GEMS_RUNTIME_PATH}/inp_image_deformmesh.mgz COPYONLY)
configure_file(atlas_level1_deformmesh.txt.gz ${GEMS_RUNTIME_PATH}/atlas_level1_deformmesh.txt.gz COPYONLY)
