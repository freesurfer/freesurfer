project(mris_pval_fill)
include_directories(${mris_pval_fill_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_pval_fill_SRCS
mris_pval_fill.c
)


add_executable(mris_pval_fill ${mris_pval_fill_SRCS})
target_link_libraries(mris_pval_fill ${FS_LIBS})
install(TARGETS mris_pval_fill DESTINATION bin)	

