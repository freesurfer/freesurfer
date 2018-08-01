project(mri_cnr)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_cnr mri_cnr.c)
target_link_libraries(mri_cnr utils)

install(TARGETS mri_cnr DESTINATION bin)
