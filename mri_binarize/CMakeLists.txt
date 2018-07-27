project(mri_binarize)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_binarize mri_binarize.c)
target_link_libraries(mri_binarize utils)

install(TARGETS mri_binarize DESTINATION bin)
