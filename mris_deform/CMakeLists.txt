project(mris_deform)
include_directories(${mris_deform_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_deform_SRCS
mris_deform.c
)


add_executable(mris_deform ${mris_deform_SRCS})
target_link_libraries(mris_deform ${FS_LIBS})
install(TARGETS mris_deform DESTINATION bin)	

