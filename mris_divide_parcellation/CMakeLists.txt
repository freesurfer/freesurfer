project(mris_divide_parcellation)
include_directories(${mris_divide_parcellation_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_divide_parcellation_SRCS
mris_divide_parcellation.c
)


add_executable(mris_divide_parcellation ${mris_divide_parcellation_SRCS})
target_link_libraries(mris_divide_parcellation ${FS_LIBS})
install(TARGETS mris_divide_parcellation DESTINATION bin)	

