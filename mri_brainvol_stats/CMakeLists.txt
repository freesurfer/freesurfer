project(mri_brainvol_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_brainvol_stats mri_brainvol_stats.cpp)
target_link_libraries(mri_brainvol_stats utils)
install(TARGETS mri_brainvol_stats DESTINATION bin) 
