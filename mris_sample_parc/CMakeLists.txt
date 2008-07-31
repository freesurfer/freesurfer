project(mris_sample_parc)
include_directories(${mris_sample_parc_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_sample_parc_SRCS
mris_sample_parc.c
)


add_executable(mris_sample_parc ${mris_sample_parc_SRCS})
target_link_libraries(mris_sample_parc ${FS_LIBS})
install(TARGETS mris_sample_parc DESTINATION bin)	

