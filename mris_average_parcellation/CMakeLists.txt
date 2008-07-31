project(mris_average_parcellation)
include_directories(${mris_average_parcellation_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_average_parcellation_SRCS
mris_average_parcellation.c
)


add_executable(mris_average_parcellation ${mris_average_parcellation_SRCS})
target_link_libraries(mris_average_parcellation ${FS_LIBS})
install(TARGETS mris_average_parcellation DESTINATION bin)	

