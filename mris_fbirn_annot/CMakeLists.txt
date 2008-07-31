project(mris_fbirn_annot)
include_directories(${mris_fbirn_annot_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_fbirn_annot_SRCS
mris_fbirn_annot.c
)


add_executable(mris_fbirn_annot ${mris_fbirn_annot_SRCS})
target_link_libraries(mris_fbirn_annot ${FS_LIBS})
install(TARGETS mris_fbirn_annot DESTINATION bin)	

