project(mri_tessellate)
include_directories(${mri_tessellate_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_tessellate_SRCS
ico_supersample.c
mri_tessellate.c
)


add_executable(mri_tessellate ${mri_tessellate_SRCS})
target_link_libraries(mri_tessellate ${FS_LIBS})
install(TARGETS mri_tessellate DESTINATION bin)	

add_subdirectory(test)

