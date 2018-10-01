project(fslutils)

include_directories(${FS_INCLUDE_DIRS})

add_executable(fsl_label2voxel fsl_label2voxel.cpp)
target_link_libraries(fsl_label2voxel utils)
install(TARGETS fsl_label2voxel DESTINATION bin)
