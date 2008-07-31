project(mri_watershed)
include_directories(${mri_watershed_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_watershed_SRCS
mri_validate_skull_stripped.cpp
mri_watershed.cpp
)


add_executable(mri_watershed ${mri_watershed_SRCS})
target_link_libraries(mri_watershed ${FS_LIBS})
install(TARGETS mri_watershed DESTINATION bin)	

add_subdirectory(brain_volume)
