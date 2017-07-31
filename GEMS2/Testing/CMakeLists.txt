# Add an executable
ADD_EXECUTABLE( kvlAtlasMeshRasterizorTestGPU
                kvlAtlasMeshRasterizorTestGPU.cxx 
                )
TARGET_LINK_LIBRARIES( kvlAtlasMeshRasterizorTestGPU
                       kvlGEMSCommon
                      )


# Boost tests
IF(Boost_FOUND)
  SET(GEMS2libs kvlGEMSCommon )
  SET(testsrcs boosttests.cpp)
  LIST(APPEND testsrcs testatlasmeshrasterizorbasic.cpp)
  LIST(APPEND testsrcs atlasmeshvisitcountercpuwrapper.cpp)
  LIST(APPEND testsrcs atlasmeshalphadrawercpuwrapper.cpp)
  LIST(APPEND testsrcs testatlasmeshvisitcounter.cpp)
  LIST(APPEND testsrcs testatlasmeshalphadrawer.cpp)
  LIST(APPEND testsrcs testdimensioncuda.cpp)
  LIST(APPEND testsrcs teststopwatch.cpp)

  IF( CUDA_FOUND )
    LIST(APPEND testsrcs testcudaimage.cpp)
    LIST(APPEND testsrcs cudaglobalfixture.cpp)
    LIST(APPEND testsrcs cudaimagetests.cu)
    LIST(APPEND testsrcs testcudatetrahedralmesh.cpp)

    LIST(APPEND GEMS2libs kvlGEMSCUDA)
  ENDIF( CUDA_FOUND )

  IF( CUDA_FOUND )
    CUDA_ADD_EXECUTABLE( TestGEMS2 ${testsrcs} ${cudatestsrcs} )
  ELSE( CUDA_FOUND )
    ADD_EXECUTABLE( TestGEMS2 ${testsrcs} )
  ENDIF(CUDA_FOUND)

  TARGET_LINK_LIBRARIES( TestGEMS2 ${Boost_LIBRARIES} ${GEMS2libs} )
ENDIF(Boost_FOUND)

                  
# Copy some files
CONFIGURE_FILE( test.nii ${GEMS_LIBRARY_PATH}/test.nii COPYONLY )
CONFIGURE_FILE( test.txt.gz ${GEMS_LIBRARY_PATH}/test.txt.gz COPYONLY )
