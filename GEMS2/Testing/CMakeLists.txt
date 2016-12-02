# Add an executable
ADD_EXECUTABLE( kvlAtlasMeshRasterizorTestGPU
                kvlAtlasMeshRasterizorTestGPU.cxx 
                )
TARGET_LINK_LIBRARIES( kvlAtlasMeshRasterizorTestGPU
                       kvlGEMSCommon
                      )
                  
# Copy some files
CONFIGURE_FILE( test.nii ${GEMS_LIBRARY_PATH}/test.nii COPYONLY )
CONFIGURE_FILE( test.txt.gz ${GEMS_LIBRARY_PATH}/test.txt.gz COPYONLY )
