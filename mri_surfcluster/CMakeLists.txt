project(mri_surfcluster)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_surfcluster mri_surfcluster.c)
target_link_libraries(mri_surfcluster utils)

install(TARGETS mri_surfcluster DESTINATION bin)
