project(mri_volcluster)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_volcluster mri_volcluster.c)
target_link_libraries(mri_volcluster utils)

install(TARGETS mri_volcluster DESTINATION bin)
