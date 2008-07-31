project(mri_parselabel)
include_directories(${mri_parselabel_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_parselabel_SRCS
mri_parselabel.cpp
)


add_executable(mri_parselabel ${mri_parselabel_SRCS})
target_link_libraries(mri_parselabel ${FS_LIBS})
install(TARGETS mri_parselabel DESTINATION bin)	

add_subdirectory(test)

