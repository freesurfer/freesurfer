project(mri_surfacemask)
include_directories(${mri_surfacemask_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_surfacemask_SRCS
mri_surfacemask.cpp
)


add_executable(mri_surfacemask ${mri_surfacemask_SRCS})
target_link_libraries(mri_surfacemask ${FS_LIBS})
install(TARGETS mri_surfacemask DESTINATION bin)	

