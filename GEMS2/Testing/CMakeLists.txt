# Add an executable
ADD_EXECUTABLE( kvlAtlasMeshRasterizorTestGPU
                kvlAtlasMeshRasterizorTestGPU.cxx 
                )
TARGET_LINK_LIBRARIES( kvlAtlasMeshRasterizorTestGPU
                       kvlGEMSCommon
                      )


# Boost tests
IF(Boost_FOUND)
  SET(testsrcs boosttests.cpp)
  LIST(APPEND testsrcs testatlasmeshrasterizorbasic.cpp)

  ADD_EXECUTABLE( TestGEMS2 ${testsrcs} )

  TARGET_LINK_LIBRARIES( TestGEMS2 ${Boost_LIBRARIES} kvlGEMSCommon )
ENDIF(Boost_FOUND)

                  
# Copy some files
CONFIGURE_FILE( test.nii ${GEMS_LIBRARY_PATH}/test.nii COPYONLY )
CONFIGURE_FILE( test.txt.gz ${GEMS_LIBRARY_PATH}/test.txt.gz COPYONLY )
