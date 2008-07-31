project(dmri_poistats)

SET (DMRI_DATAMODEL_SRCS
	datamodel/PoistatsModel.cxx 
	datamodel/PoistatsModel.h
	)

SET (DMRI_DATA_SRCS
	datamodel/data/PoistatsReplica.cxx 
	datamodel/data/PoistatsReplica.h 
	datamodel/data/PoistatsReplicas.cxx 
	datamodel/data/PoistatsReplicas.h
	)

SET (DMRI_EVENTS_SRCS
  datamodel/events/PoistatsEvents.h 
  datamodel/events/CommandUpdate.cxx 
  datamodel/events/CommandUpdate.h
 )

SET (DMRI_UTILS_SRCS
	datamodel/utils/itkPoistatsFilter.txx 
	datamodel/utils/itkPoistatsFilter.h 
	datamodel/utils/InitializePath.cxx 
	datamodel/utils/InitializePath.h 
	datamodel/utils/EigenVectorInitPathStrategy.cxx 
	datamodel/utils/EigenVectorInitPathStrategy.h 
	datamodel/utils/FieldLineInitPathStrategy.cxx 
	datamodel/utils/FieldLineInitPathStrategy.h	
)

SET (DMRI_IO_SRCS 
	io/AsymmetricTensorReaderStrategy.cxx 
	io/AsymmetricTensorReaderStrategy.h 
	io/AsymmetricTensorVectorReaderStrategy.cxx 
	io/AsymmetricTensorVectorReaderStrategy.h 
	io/SymmetricTensorReaderStrategy.cxx 
	io/SymmetricTensorReaderStrategy.h 
	io/TensorReaderStrategy.cxx 
	io/TensorReaderStrategy.h
)

# poistats
SET(dmri_poistats_SOURCES  
  ${DMRI_DATAMODEL_SRCS}
  ${DMRI_DATA_SRCS}
  ${DMRI_EVENTS_SRCS}
  ${DMRI_UTILS_SRCS}
  ${DMRI_IO_SRCS}
  dmri_poistats.cxx 
  )
  
SET (dmri_poistats_LDADD 
itkutils
ITKBasicFilters 
litkgdcm 
litkpng 
litktiff 
lITKSpatialObject
lITKNumerics
litkvnl_inst
litkvnl_algo
litkv3p_netlib
litkvnl
litkvcl
lITKMetaIO 
litksys 
lITKEXPAT 
lITKniftiio 
lITKznz 
litkzlib
${FS_LIBS}
)
add_executable (dmri_poistats ${dmri_poistats_SOURCES})
target_link_libraries(dmri_poistats ${dmri_poistats_LDADD})
install(TARGETS dmri_poistats DESTINATION bin)	

SET( dmri_initialize_path_SOURCES 
  ${DATAMODEL_SRCS}
  ${UTILS_SRCS}
  dmri_initialize_path.cxx 
  )
  
SET (dmri_initialize_path_LDADD 
itkutils
ITKBasicFilters 
litkgdcm 
litkpng 
litktiff 
lITKSpatialObject
lITKNumerics
litkvnl_inst
litkvnl_algo
litkv3p_netlib
litkvnl
litkvcl
lITKMetaIO 
litksys 
lITKEXPAT 
lITKniftiio 
lITKznz 
litkzlib
${FS_LIBS}
)
add_executable (dmri_initialize_path ${dmri_initialize_path_SOURCES})
target_link_libraries(dmri_initialize_path ${dmri_initialize_path_LDADD})
install(TARGETS dmri_initialize_path DESTINATION bin)	

# mask path
SET (dmri_mask_SOURCES
  dmri_mask.cxx 
)

SET (dmri_mask_LDADD 
	${FS_LIBS}
	itkutils
)

add_executable (dmri_mask ${dmri_mask_SOURCES})
target_link_libraries(dmri_mask ${dmri_mask_LDADD})
install(TARGETS dmri_mask DESTINATION bin)	

