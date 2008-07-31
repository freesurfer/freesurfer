project(mris_glm)
include_directories(${mris_glm_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_glm_SRCS
mris_glm.c
)

SET (MRIS_GLM_LIBS
${FS_LIBS}
fsgdf
)

add_executable(mris_glm ${mris_glm_SRCS})
target_link_libraries(mris_glm ${MRIS_GLM_LIBS})
install(TARGETS mris_glm DESTINATION bin)	

