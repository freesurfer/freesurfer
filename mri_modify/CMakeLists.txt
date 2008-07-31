project(mri_modify)
include_directories(${mri_modify_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_modify_SRCS
mri_copy_params.cpp
mri_modify.cpp
)


add_executable(mri_modify ${mri_modify_SRCS})
target_link_libraries(mri_modify ${FS_LIBS})
install(TARGETS mri_modify DESTINATION bin)	

