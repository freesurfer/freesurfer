project(mri_probe_ima)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_probe_ima mri_probe_ima.cpp)
target_link_libraries(mri_probe_ima utils)

install(TARGETS mri_probe_ima DESTINATION bin)
